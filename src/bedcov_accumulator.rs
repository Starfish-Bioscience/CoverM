use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use memmap2::Mmap;

use rayon::prelude::*;

use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam;

use crate::mosdepth_genome_coverage_estimators::{
    CoverageEstimator, MosdepthGenomeCoverageEstimator,
};

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct BedRegion {
    pub start: u32,
    pub end: u32,
    pub label_idx: usize,
    pub global_idx: usize,
}

#[derive(Debug, Clone)]
pub struct RegionInfo {
    pub chrom_idx: u32, // index into ParsedBed.chrom_table
    pub start: u32,
    pub end: u32,
    pub label_idx: usize,
}

/// Parsed BED file: labels (alphabetical), regions indexed two ways.
pub struct ParsedBed {
    /// chrom → regions sorted by start
    pub regions_by_chrom: HashMap<String, Vec<BedRegion>>,
    /// alphabetical list of distinct labels (real labels first, unlabeled pseudo-label last)
    pub labels: Vec<String>,
    pub label_to_idx: HashMap<String, usize>,
    /// all_regions[global_idx] — original BED order (Some only when --output-bedcov requested)
    pub all_regions: Option<Vec<RegionInfo>>,
    /// length of each region in bp, parallel to all_regions (Some only when --output-bedcov requested)
    pub region_lengths: Option<Vec<u32>>,
    /// ordered list of distinct chrom names; region.chrom_idx indexes into this vec
    pub chrom_table: Vec<String>,
    /// true when --regions-bed-unlabeled was requested
    pub with_unlabeled: bool,
    /// index of the synthetic unlabeled label (= labels.len()-1 when with_unlabeled)
    pub unlabeled_label_idx: Option<usize>,
}

impl ParsedBed {
    pub fn from_file(path: &str, unlabeled_label: Option<&str>, need_bedgraph: bool) -> ParsedBed {
        let file = fs::File::open(path)
            .unwrap_or_else(|e| panic!("Cannot open --regions-bed file '{}': {}", path, e));
        let mmap = unsafe { Mmap::map(&file) }
            .unwrap_or_else(|e| panic!("Cannot mmap --regions-bed file '{}': {}", path, e));
        let text = std::str::from_utf8(&mmap)
            .unwrap_or_else(|e| panic!("--regions-bed file '{}' is not valid UTF-8: {}", path, e));

        // --- first pass: collect raw records (parallel parse) ---
        struct RawRecord {
            chrom: String,
            start: u32,
            end: u32,
            label: String,
        }

        // Collect non-empty, non-comment lines with their original 1-based line numbers.
        let data_lines: Vec<(usize, &str)> = text
            .lines()
            .enumerate()
            .filter(|(_, l)| {
                let t = l.trim();
                !t.is_empty() && !t.starts_with('#')
            })
            .collect();

        // Parse each line in parallel; panics propagate via unwrap on the join.
        let raw: Vec<RawRecord> = data_lines
            .into_par_iter()
            .map(|(lineno, line)| {
                let line = line.trim();
                let cols: Vec<&str> = line.splitn(5, '\t').collect();
                if cols.len() < 4 {
                    panic!(
                        "--regions-bed '{}' line {}: expected at least 4 tab-separated columns, got {}",
                        path,
                        lineno + 1,
                        cols.len()
                    );
                }
                let chrom = cols[0].to_string();
                let start: u32 = cols[1].parse().unwrap_or_else(|_| {
                    panic!(
                        "--regions-bed '{}' line {}: cannot parse start '{}'",
                        path,
                        lineno + 1,
                        cols[1]
                    )
                });
                let end: u32 = cols[2].parse().unwrap_or_else(|_| {
                    panic!(
                        "--regions-bed '{}' line {}: cannot parse end '{}'",
                        path,
                        lineno + 1,
                        cols[2]
                    )
                });
                if end <= start {
                    panic!(
                        "--regions-bed '{}' line {}: end ({}) must be > start ({})",
                        path,
                        lineno + 1,
                        end,
                        start
                    );
                }
                let label = cols[3].to_string();
                RawRecord {
                    chrom,
                    start,
                    end,
                    label,
                }
            })
            .collect();

        // Rebuild label_set from the parsed records (sequential, cheap).
        let mut label_set: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
        for r in &raw {
            label_set.insert(r.label.clone());
        }

        // --- build label index (alphabetical) ---
        let mut labels: Vec<String> = label_set.into_iter().collect();
        let mut label_to_idx: HashMap<String, usize> = labels
            .iter()
            .enumerate()
            .map(|(i, l)| (l.clone(), i))
            .collect();

        // --- optional synthetic unlabeled label (appended after alphabetical labels) ---
        let (with_unlabeled, unlabeled_label_idx) = match unlabeled_label {
            None => (false, None),
            Some(name) => {
                if label_to_idx.contains_key(name) {
                    panic!(
                        "--regions-bed-unlabeled label '{}' conflicts with an existing \
                         label in '{}'; choose a different name",
                        name, path
                    );
                }
                let idx = labels.len();
                labels.push(name.to_string());
                label_to_idx.insert(name.to_string(), idx);
                (true, Some(idx))
            }
        };

        // --- build all_regions, regions_by_chrom ---
        let n_raw = raw.len();
        let mut all_regions: Option<Vec<RegionInfo>> = if need_bedgraph {
            Some(Vec::with_capacity(n_raw))
        } else {
            None
        };
        let mut region_lengths: Option<Vec<u32>> = if need_bedgraph {
            Some(Vec::with_capacity(n_raw))
        } else {
            None
        };
        let mut regions_by_chrom: HashMap<String, Vec<BedRegion>> = HashMap::new();
        // Intern chrom strings: deduplicate O(N) heap allocations down to O(distinct chroms).
        let mut chrom_intern: HashMap<String, u32> = HashMap::new();
        let mut chrom_table: Vec<String> = Vec::new();

        for (global_idx, r) in raw.into_iter().enumerate() {
            let label_idx = label_to_idx[&r.label];
            let chrom_idx = *chrom_intern.entry(r.chrom.clone()).or_insert_with(|| {
                let idx = chrom_table.len() as u32;
                chrom_table.push(r.chrom.clone());
                idx
            });
            if let Some(ref mut ar) = all_regions {
                ar.push(RegionInfo {
                    chrom_idx,
                    start: r.start,
                    end: r.end,
                    label_idx,
                });
            }
            if let Some(ref mut rl) = region_lengths {
                rl.push(r.end - r.start);
            }
            regions_by_chrom
                .entry(r.chrom)
                .or_default()
                .push(BedRegion {
                    start: r.start,
                    end: r.end,
                    label_idx,
                    global_idx,
                });
        }

        // sort per-chrom regions by start
        for v in regions_by_chrom.values_mut() {
            v.sort_unstable_by_key(|r| r.start);
        }

        // Count overlapping region pairs; emit a single summary warn instead of
        // one message per pair (large BED files can produce millions of pairs).
        let mut overlap_contigs: usize = 0;
        let mut overlap_pairs: usize = 0;
        for v in regions_by_chrom.values() {
            let mut contig_had_overlap = false;
            for w in v.windows(2) {
                if w[1].start < w[0].end {
                    overlap_pairs += 1;
                    if !contig_had_overlap {
                        overlap_contigs += 1;
                        contig_had_overlap = true;
                    }
                }
            }
        }
        if overlap_pairs > 0 {
            warn!(
                "Overlapping BED regions detected: {} overlapping pair(s) across {} contig(s); \
                 each region is counted independently for its own label, \
                 and their union is used for the unlabeled complement",
                overlap_pairs, overlap_contigs
            );
        }

        info!(
            "Loaded {} regions across {} distinct label(s) from '{}'",
            n_raw,
            labels.len(),
            path
        );

        ParsedBed {
            regions_by_chrom,
            labels,
            label_to_idx,
            all_regions,
            region_lengths,
            chrom_table,
            with_unlabeled,
            unlabeled_label_idx,
        }
    }
}

// ---------------------------------------------------------------------------
// make_label_estimators — create fresh per-region estimators from templates
// ---------------------------------------------------------------------------

/// Build a set of fresh `CoverageEstimator`s suitable for per-region
/// accumulation, derived from the user-specified method templates.
///
/// All estimators are forced to `contig_end_exclusion = 0` and
/// `min_fraction_covered_bases = 0.0`, because:
/// - `--contig-end-exclusion` is CLI-incompatible with `--regions-bed`
/// - Minimum-fraction filtering does not make sense per-region
///
/// `coverage_histogram`, `anir`, and `strobealign-aemb` are rejected at the
/// CLI level before this function is called.
pub fn make_label_estimators(templates: &[CoverageEstimator]) -> Vec<CoverageEstimator> {
    templates
        .iter()
        .map(|e| match e {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                exclude_mismatches, ..
            } => CoverageEstimator::new_estimator_mean(0.0, 0, *exclude_mismatches),
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator { min, max, .. } => {
                CoverageEstimator::new_estimator_trimmed_mean(*min, *max, 0.0, 0)
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator { .. } => {
                CoverageEstimator::new_estimator_covered_fraction(0.0)
            }
            CoverageEstimator::NumCoveredBasesCoverageEstimator { .. } => {
                CoverageEstimator::new_estimator_covered_bases(0.0)
            }
            CoverageEstimator::VarianceGenomeCoverageEstimator { .. } => {
                CoverageEstimator::new_estimator_variance(0.0, 0)
            }
            CoverageEstimator::RPKMCoverageEstimator { .. } => {
                CoverageEstimator::new_estimator_rpkm(0.0)
            }
            CoverageEstimator::TPMCoverageEstimator { .. } => {
                CoverageEstimator::new_estimator_tpm(0.0)
            }
            CoverageEstimator::ReferenceLengthCalculator { .. } => {
                CoverageEstimator::new_estimator_length()
            }
            CoverageEstimator::ReadCountCalculator { .. } => {
                CoverageEstimator::new_estimator_read_count()
            }
            CoverageEstimator::ReadsPerBaseCalculator { .. } => {
                CoverageEstimator::new_estimator_reads_per_base()
            }
            CoverageEstimator::PileupCountsGenomeCoverageEstimator { .. } => {
                unreachable!("coverage_histogram is incompatible with --regions-bed")
            }
            CoverageEstimator::AverageIdentityEstimator { .. } => {
                unreachable!("anir is incompatible with --regions-bed")
            }
            CoverageEstimator::StrobealignAembEstimator { .. } => {
                unreachable!("strobealign-aemb is incompatible with --regions-bed")
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// BedIndex — per-BAM tid → regions mapping (rebuilt for each BAM)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct IndexedRegion {
    pub start: u32,
    pub end: u32,
    pub label_idx: usize,
    pub global_idx: usize,
}

pub struct BedIndex {
    /// tid → sorted regions on that contig
    pub regions_by_tid: HashMap<u32, Vec<IndexedRegion>>,
    /// tid → index in ParsedBed.chrom_table (only for tids that have BED regions)
    pub tid_to_chrom_idx: HashMap<u32, u32>,
}

/// One unlabeled gap interval recorded for bedGraph output.
struct UnlabeledEntry {
    chrom_idx: u32,
    start: u32,
    end: u32,
    sum: i64,
}

// ---------------------------------------------------------------------------
// BedcovAccumulator — main state
// ---------------------------------------------------------------------------

pub struct BedcovAccumulator {
    pub parsed_bed: ParsedBed,
    pub current_index: BedIndex,
    /// sum of coverage per region (reset at start of each BAM); None when --output-bedcov not requested
    pub region_sums: Option<Vec<i64>>,
    /// Feature 2 — method templates (with end_excl=0, min_frac=0)
    pub label_estimator_templates: Vec<CoverageEstimator>,
    /// Feature 2 — per-genome per-label per-method estimators (lazily created)
    /// Key: genome_name, Value: [label_idx][method_idx]
    pub genome_label_estimators: HashMap<String, Vec<Vec<CoverageEstimator>>>,
    /// prefix-sum scratch (reused, never shrunk)
    pcov_scratch: Vec<i64>,
    /// sub_ups scratch for Feature 2 (reused, never shrunk)
    sub_ups_scratch: Vec<i32>,
    /// merged-interval scratch for unlabeled gap computation (reused, never shrunk)
    merged_scratch: Vec<(usize, usize)>,
    /// Gap intervals accumulated for bedGraph unlabeled track (None when not needed)
    unlabeled_bg_entries: Option<Vec<UnlabeledEntry>>,
}

impl BedcovAccumulator {
    pub fn new(
        parsed_bed: ParsedBed,
        label_estimator_templates: Vec<CoverageEstimator>,
    ) -> BedcovAccumulator {
        let region_sums = parsed_bed
            .all_regions
            .as_ref()
            .map(|ar| vec![0i64; ar.len()]);
        let unlabeled_bg_entries = if parsed_bed.all_regions.is_some() && parsed_bed.with_unlabeled
        {
            Some(Vec::new())
        } else {
            None
        };
        BedcovAccumulator {
            parsed_bed,
            current_index: BedIndex {
                regions_by_tid: HashMap::new(),
                tid_to_chrom_idx: HashMap::new(),
            },
            region_sums,
            label_estimator_templates,
            genome_label_estimators: HashMap::new(),
            pcov_scratch: Vec::new(),
            sub_ups_scratch: Vec::new(),
            merged_scratch: Vec::new(),
            unlabeled_bg_entries,
        }
    }

    /// Rebuild BedIndex from the BAM header (call once per BAM, before processing reads).
    pub fn reinit_for_bam(&mut self, header: &bam::HeaderView) {
        if let Some(ref mut sums) = self.region_sums {
            sums.fill(0);
        }
        if let Some(ref mut entries) = self.unlabeled_bg_entries {
            entries.clear();
        }

        // Build chrom_name → chrom_idx lookup once (used for tid_to_chrom_idx below).
        let need_chrom_idx = self.unlabeled_bg_entries.is_some();
        let chrom_name_to_idx: HashMap<&str, u32> = if need_chrom_idx {
            self.parsed_bed
                .chrom_table
                .iter()
                .enumerate()
                .map(|(i, s)| (s.as_str(), i as u32))
                .collect()
        } else {
            HashMap::new()
        };

        let mut regions_by_tid: HashMap<u32, Vec<IndexedRegion>> = HashMap::new();
        let mut tid_to_chrom_idx: HashMap<u32, u32> = HashMap::new();

        for (tid, name_bytes) in header.target_names().iter().enumerate() {
            let name = match std::str::from_utf8(name_bytes) {
                Ok(s) => s,
                Err(_) => continue,
            };
            if let Some(bed_regions) = self.parsed_bed.regions_by_chrom.get(name) {
                let indexed: Vec<IndexedRegion> = bed_regions
                    .iter()
                    .map(|r| IndexedRegion {
                        start: r.start,
                        end: r.end,
                        label_idx: r.label_idx,
                        global_idx: r.global_idx,
                    })
                    .collect();
                regions_by_tid.insert(tid as u32, indexed);
                if let Some(&cidx) = chrom_name_to_idx.get(name) {
                    tid_to_chrom_idx.insert(tid as u32, cidx);
                }
            }
        }
        self.current_index = BedIndex {
            regions_by_tid,
            tid_to_chrom_idx,
        };
        self.genome_label_estimators.clear();
    }

    /// Accumulate coverage sums for all BED regions on this contig.
    ///
    /// `ups_and_downs` is the mosdepth difference array: coverage at position i =
    /// prefix_sum(ups_and_downs[0..=i]).
    ///
    /// `genome_name` is used to key the per-genome label estimators (Feature 2).
    pub fn process_contig(&mut self, tid: u32, genome_name: &str, ups_and_downs: &[i32]) {
        let has_bed_regions = self.current_index.regions_by_tid.contains_key(&tid);
        // Unlabeled gaps only matter when Feature 2 estimators are active.
        let need_unlabeled =
            self.parsed_bed.with_unlabeled && !self.label_estimator_templates.is_empty();

        if !has_bed_regions && !need_unlabeled {
            return;
        }

        let n = ups_and_downs.len();

        // Build prefix-coverage array: pcov[i+1] - pcov[i] = absolute_coverage[i]
        // pcov[e] - pcov[s] = sum of coverage over [s, e)
        if self.pcov_scratch.len() < n + 1 {
            self.pcov_scratch.resize(n + 1, 0);
        }
        self.pcov_scratch[0] = 0;
        let mut running: i64 = 0;
        for (i, &delta) in ups_and_downs.iter().enumerate() {
            running += delta as i64;
            self.pcov_scratch[i + 1] = self.pcov_scratch[i] + running;
        }

        // Feature 2 — lazily create the genome entry if templates are configured
        if !self.label_estimator_templates.is_empty()
            && !self.genome_label_estimators.contains_key(genome_name)
        {
            let n_labels = self.parsed_bed.labels.len();
            let templates = self.label_estimator_templates.clone();
            let fresh: Vec<Vec<CoverageEstimator>> =
                (0..n_labels).map(|_| templates.clone()).collect();
            self.genome_label_estimators
                .insert(genome_name.to_string(), fresh);
        }

        // Take regions out of the HashMap (zero-cost move).
        // If the contig has no BED regions, use an empty Vec.
        let regions = if has_bed_regions {
            std::mem::take(self.current_index.regions_by_tid.get_mut(&tid).unwrap())
        } else {
            Vec::new()
        };

        // --- Labeled regions (Feature 1 + Feature 2) ---
        for region in &regions {
            let s = (region.start as usize).min(n);
            let e = (region.end as usize).min(n);
            if s >= e {
                if region.end as usize > n {
                    warn!(
                        "BED region [{},{}) on contig tid={} extends beyond contig length {}; \
                         region is skipped",
                        region.start, region.end, tid, n
                    );
                }
                continue;
            }

            // Feature 1 — accumulate region coverage sum (only when bedgraph output requested)
            let delta = self.pcov_scratch[e] - self.pcov_scratch[s];
            if let Some(ref mut sums) = self.region_sums {
                sums[region.global_idx] += delta;
            }

            // Feature 2 — call add_contig on per-label estimators via sub_ups.
            // We need to borrow both `sub_ups_scratch` (immutably, after writing)
            // and `genome_label_estimators` (mutably).  Destructuring `self` into
            // its fields makes the disjointness visible to the borrow checker.
            if !self.label_estimator_templates.is_empty() {
                let len = e - s;
                if self.sub_ups_scratch.len() < len {
                    self.sub_ups_scratch.resize(len, 0);
                }
                // sub_ups[0] = absolute coverage at position s
                self.sub_ups_scratch[0] = (self.pcov_scratch[s + 1] - self.pcov_scratch[s]) as i32;
                // sub_ups[i] = ups_and_downs[s+i] for i >= 1
                // Invariant: prefix_sum(sub_ups)[j] = coverage[s+j]  ✓
                self.sub_ups_scratch[1..len].copy_from_slice(&ups_and_downs[(s + 1)..(s + len)]);
                // Destructure to let the borrow checker see disjoint fields.
                let BedcovAccumulator {
                    sub_ups_scratch,
                    genome_label_estimators,
                    ..
                } = self;
                let sub_slice = &sub_ups_scratch[..len];
                let label_ests = genome_label_estimators.get_mut(genome_name).unwrap();
                for est in label_ests[region.label_idx].iter_mut() {
                    est.add_contig(sub_slice, 0, 0, 0.0);
                }
            }
        }

        // --- Unlabeled gaps (Feature 2 only) ---
        //
        // Build the merged union of all labeled regions on this contig (across all
        // labels, regardless of which label), then feed the complement intervals to
        // the unlabeled estimator.  Merging is required to avoid under-counting gaps
        // when regions from different labels overlap.
        if need_unlabeled {
            let unlabeled_idx = self.parsed_bed.unlabeled_label_idx.unwrap();

            // Merge overlapping/adjacent intervals.  `regions` is already sorted by
            // start (guaranteed by `reinit_for_bam` → `regions_by_chrom` sort).
            // Reuse the scratch buffer to avoid a heap allocation per contig.
            self.merged_scratch.clear();
            for region in &regions {
                let s = (region.start as usize).min(n);
                let e = (region.end as usize).min(n);
                if s >= e {
                    continue;
                }
                if let Some(last) = self.merged_scratch.last_mut() {
                    if s < last.1 {
                        // overlapping or adjacent — extend
                        last.1 = last.1.max(e);
                    } else {
                        self.merged_scratch.push((s, e));
                    }
                } else {
                    self.merged_scratch.push((s, e));
                }
            }

            // Iterate over gaps between merged intervals using an index loop.
            // Reading `self.merged_scratch[i]` as a Copy value (usize, usize) releases
            // the borrow before the destructuring of `self` below — no clone needed.
            // A virtual sentinel at index n_merged represents (n, n) so the final gap
            // [last_me, n) is handled by the same loop body.
            let n_merged = self.merged_scratch.len();
            let mut gap_start = 0usize;
            for i in 0..=n_merged {
                let (ms, me) = if i < n_merged {
                    self.merged_scratch[i]
                } else {
                    (n, n)
                };
                if gap_start < ms {
                    // Gap interval: [gap_start, ms)
                    let gs = gap_start;
                    let len = ms - gs;

                    // Record for bedGraph unlabeled track (before self destructure).
                    // pcov_scratch[ms] - pcov_scratch[gs] is the total coverage sum over [gs, ms).
                    let chrom_idx_opt = self.current_index.tid_to_chrom_idx.get(&tid).copied();
                    if let (Some(ref mut entries), Some(chrom_idx)) =
                        (&mut self.unlabeled_bg_entries, chrom_idx_opt)
                    {
                        entries.push(UnlabeledEntry {
                            chrom_idx,
                            start: gs as u32,
                            end: ms as u32,
                            sum: self.pcov_scratch[ms] - self.pcov_scratch[gs],
                        });
                    }

                    if self.sub_ups_scratch.len() < len {
                        self.sub_ups_scratch.resize(len, 0);
                    }
                    self.sub_ups_scratch[0] =
                        (self.pcov_scratch[gs + 1] - self.pcov_scratch[gs]) as i32;
                    if len > 1 {
                        self.sub_ups_scratch[1..len]
                            .copy_from_slice(&ups_and_downs[(gs + 1)..(gs + len)]);
                    }
                    let BedcovAccumulator {
                        sub_ups_scratch,
                        genome_label_estimators,
                        ..
                    } = self;
                    let sub_slice = &sub_ups_scratch[..len];
                    let label_ests = genome_label_estimators.get_mut(genome_name).unwrap();
                    for est in label_ests[unlabeled_idx].iter_mut() {
                        est.add_contig(sub_slice, 0, 0, 0.0);
                    }
                }
                gap_start = me;
            }
        }

        // Restore regions into the HashMap (zero-cost move, no allocation).
        if has_bed_regions {
            if let Some(slot) = self.current_index.regions_by_tid.get_mut(&tid) {
                *slot = regions;
            }
        }
    }

    /// Retrieve (and reset) per-label coverage values for a genome.
    ///
    /// Returns a `[label_idx][method_idx]` matrix of `f32` coverage values.
    /// Returns an empty vec when `label_estimator_templates` is empty.
    /// Returns all-zero rows when the genome had no BED regions.
    pub fn take_label_coverages(&mut self, genome_name: &str) -> Vec<Vec<f32>> {
        if self.label_estimator_templates.is_empty() {
            return vec![];
        }
        match self.genome_label_estimators.remove(genome_name) {
            None => {
                // Genome had no BED regions — return zeros for all labels/methods
                let n_labels = self.parsed_bed.labels.len();
                let n_methods = self.label_estimator_templates.len();
                vec![vec![0.0f32; n_methods]; n_labels]
            }
            Some(mut label_ests) => label_ests
                .iter_mut()
                .map(|method_ests| {
                    method_ests
                        .iter_mut()
                        .map(|est| est.calculate_coverage(&[]))
                        .collect()
                })
                .collect(),
        }
    }

    /// Write a multi-track bedGraph file for the current BAM sample.
    ///
    /// One track per label (alphabetical). Regions with zero coverage are included.
    pub fn finalise_bedcov(&self, sample_name: &str, output_dir: &Path, compress: bool) {
        let safe_name = sanitize_filename(sample_name);
        let filename = if compress {
            format!("{}.bedgraph.gz", safe_name)
        } else {
            format!("{}.bedgraph", safe_name)
        };
        let out_path = output_dir.join(&filename);

        let file = fs::File::create(&out_path).unwrap_or_else(|e| {
            panic!(
                "Cannot create bedGraph file '{}': {}",
                out_path.display(),
                e
            )
        });
        let buf = BufWriter::new(file);

        if compress {
            let mut gz = GzEncoder::new(buf, Compression::default());
            self.write_bedgraph(&mut gz, sample_name);
            gz.finish()
                .unwrap_or_else(|e| panic!("Error finishing gzip stream '{}': {}", filename, e));
        } else {
            let mut w = buf;
            self.write_bedgraph(&mut w, sample_name);
            w.flush()
                .unwrap_or_else(|e| panic!("Error flushing '{}': {}", filename, e));
        }

        info!("Wrote bedGraph to '{}'", out_path.display());
    }

    fn write_bedgraph<W: Write>(&self, w: &mut W, _sample_name: &str) {
        let bed = &self.parsed_bed;
        let all_regions = bed
            .all_regions
            .as_ref()
            .expect("write_bedgraph called without all_regions (need_bedgraph was false)");
        let region_lengths = bed
            .region_lengths
            .as_ref()
            .expect("write_bedgraph called without region_lengths");
        let region_sums = self
            .region_sums
            .as_ref()
            .expect("write_bedgraph called without region_sums");

        for (label_idx, label) in bed.labels.iter().enumerate() {
            writeln!(w, "track type=bedGraph name=\"{}\"", label)
                .expect("Error writing bedGraph track header");

            // The unlabeled pseudo-label has no entries in all_regions.
            // Its gap intervals are stored separately in unlabeled_bg_entries.
            if Some(label_idx) == bed.unlabeled_label_idx {
                if let Some(ref entries) = self.unlabeled_bg_entries {
                    for e in entries {
                        let len = (e.end - e.start) as i64;
                        let mean_cov = if len > 0 {
                            e.sum as f64 / len as f64
                        } else {
                            0.0
                        };
                        let chrom = &bed.chrom_table[e.chrom_idx as usize];
                        writeln!(w, "{}\t{}\t{}\t{:.6}", chrom, e.start, e.end, mean_cov)
                            .expect("Error writing bedGraph record");
                    }
                }
                continue;
            }

            for (global_idx, region) in all_regions.iter().enumerate() {
                if region.label_idx != label_idx {
                    continue;
                }
                let len = region_lengths[global_idx] as i64;
                let mean_cov = if len > 0 {
                    region_sums[global_idx] as f64 / len as f64
                } else {
                    0.0
                };
                let chrom = &bed.chrom_table[region.chrom_idx as usize];
                writeln!(
                    w,
                    "{}\t{}\t{}\t{:.6}",
                    chrom, region.start, region.end, mean_cov
                )
                .expect("Error writing bedGraph record");
            }
        }
    }
}

/// Replace characters that are unsafe in POSIX filenames.
/// `/` → `__`, spaces and other non-alphanumeric/non-dot/non-hyphen → `_`.
pub fn sanitize_filename(name: &str) -> String {
    let mut out = String::with_capacity(name.len());
    let mut prev_slash = false;
    for c in name.chars() {
        if c == '/' {
            if !prev_slash {
                out.push_str("__");
            }
            prev_slash = true;
        } else {
            prev_slash = false;
            if c.is_alphanumeric() || c == '.' || c == '-' || c == '_' {
                out.push(c);
            } else {
                out.push('_');
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_bed_file(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    fn make_accumulator(bed_content: &str) -> BedcovAccumulator {
        let f = make_bed_file(bed_content);
        // need_bedgraph=true so Feature-1 tests can inspect region_sums directly.
        let parsed = ParsedBed::from_file(f.path().to_str().unwrap(), None, true);
        BedcovAccumulator::new(parsed, vec![])
    }

    fn make_accumulator_unlabeled(
        bed_content: &str,
        templates: Vec<CoverageEstimator>,
        unlabeled_name: &str,
    ) -> BedcovAccumulator {
        let f = make_bed_file(bed_content);
        let parsed = ParsedBed::from_file(f.path().to_str().unwrap(), Some(unlabeled_name), false);
        BedcovAccumulator::new(parsed, templates)
    }

    /// Build a fake BedIndex directly (no real BAM header needed).
    fn fake_index_tid0(acc: &mut BedcovAccumulator) {
        let mut regions_by_tid: HashMap<u32, Vec<IndexedRegion>> = HashMap::new();
        let bed = &acc.parsed_bed;
        let indexed: Vec<IndexedRegion> = bed
            .regions_by_chrom
            .values()
            .flatten()
            .map(|r| IndexedRegion {
                start: r.start,
                end: r.end,
                label_idx: r.label_idx,
                global_idx: r.global_idx,
            })
            .collect();
        regions_by_tid.insert(0, indexed);
        acc.current_index = BedIndex {
            regions_by_tid,
            tid_to_chrom_idx: HashMap::new(),
        };
    }

    #[test]
    fn test_process_contig_basic() {
        // Region [2, 7) on a 10-bp contig.
        // ups_and_downs = [0,0,1,0,0,0,-1,0,0,0]
        // cov[i] = prefix_sum(ups)[i]: 0 0 1 1 1 1 0 0 0 0
        // sum over [2,7) = cov[2]+cov[3]+cov[4]+cov[5]+cov[6] = 1+1+1+1+0 = 4
        let mut acc = make_accumulator("chr1\t2\t7\tlabel_a\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![0, 0, 1, 0, 0, 0, -1, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);
        assert_eq!(acc.region_sums.as_ref().unwrap()[0], 4);
    }

    #[test]
    fn test_region_clamp() {
        // Region [8, 15) on a 10-bp contig — should clamp to [8, 10), no panic.
        let mut acc = make_accumulator("chr1\t8\t15\tlabel_a\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![0, 0, 1, 0, 0, 0, -1, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);
        // coverage at [8,10) = 0, so sum = 0
        assert_eq!(acc.region_sums.as_ref().unwrap()[0], 0);
    }

    #[test]
    fn test_zero_coverage() {
        let mut acc = make_accumulator("chr1\t0\t10\tlabel_a\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![0; 10];
        acc.process_contig(0, "g1", &ups);
        assert_eq!(acc.region_sums.as_ref().unwrap()[0], 0);
    }

    #[test]
    fn test_multiple_contigs_same_genome() {
        // Two regions on the same contig (tid=0), processed in two calls
        // (simulating two different BAM contigs mapped to same BED chrom is not
        // possible with tid=0 twice, so here we test two regions in one call).
        let mut acc = make_accumulator("chr1\t0\t3\tlabel_a\nchr1\t5\t8\tlabel_a\n");
        fake_index_tid0(&mut acc);
        // coverage = [1,1,1,0,0,2,2,2,0,0]
        let ups: Vec<i32> = vec![1, 0, 0, -1, 0, 2, 0, 0, -2, 0];
        acc.process_contig(0, "g1", &ups);
        // region0=[0,3): sum=3, region1=[5,8): sum=6
        assert_eq!(acc.region_sums.as_ref().unwrap()[0], 3);
        assert_eq!(acc.region_sums.as_ref().unwrap()[1], 6);
    }

    #[test]
    fn test_overlapping_regions() {
        // [0,5) and [3,8) — overlap, each counts independently
        let mut acc = make_accumulator("chr1\t0\t5\tlabel_a\nchr1\t3\t8\tlabel_b\n");
        fake_index_tid0(&mut acc);
        // uniform coverage = 1 everywhere
        let ups: Vec<i32> = vec![1, 0, 0, 0, 0, 0, 0, 0, -1, 0];
        acc.process_contig(0, "g1", &ups);
        assert_eq!(acc.region_sums.as_ref().unwrap()[0], 5); // [0,5)
        assert_eq!(acc.region_sums.as_ref().unwrap()[1], 5); // [3,8)
    }

    #[test]
    fn test_sanitize_filename() {
        assert_eq!(sanitize_filename("sample/ref"), "sample__ref");
        assert_eq!(sanitize_filename("my sample"), "my_sample");
        assert_eq!(sanitize_filename("a//b"), "a__b");
        assert_eq!(sanitize_filename("ok-name_1.2"), "ok-name_1.2");
    }

    #[test]
    fn test_label_absent_gives_zero_sum() {
        // Two labels; only label_a has coverage
        let mut acc = make_accumulator("chr1\t0\t5\tlabel_a\nchr1\t5\t10\tlabel_b\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![1, 0, 0, 0, 0, -1, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);
        assert_eq!(acc.region_sums.as_ref().unwrap()[0], 5); // label_a
        assert_eq!(acc.region_sums.as_ref().unwrap()[1], 0); // label_b (no coverage)
    }

    #[test]
    fn test_finalise_bedgraph_content() {
        let mut acc = make_accumulator("chr1\t0\t5\tlabel_a\nchr1\t5\t10\tlabel_b\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![2, 0, 0, 0, 0, -2, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let dir = tempfile::tempdir().unwrap();
        acc.finalise_bedcov("mysample", dir.path(), false);

        let out_path = dir.path().join("mysample.bedgraph");
        let content = std::fs::read_to_string(&out_path).unwrap();
        assert!(content.contains("track type=bedGraph name=\"label_a\""));
        assert!(content.contains("track type=bedGraph name=\"label_b\""));
        // label_a: mean = 10/5 = 2.0
        assert!(content.contains("chr1\t0\t5\t2.000000"));
        // label_b: mean = 0/5 = 0.0
        assert!(content.contains("chr1\t5\t10\t0.000000"));
    }

    // -----------------------------------------------------------------------
    // Feature 2 integration tests
    // -----------------------------------------------------------------------

    fn make_accumulator_with_estimators(
        bed_content: &str,
        templates: Vec<CoverageEstimator>,
    ) -> BedcovAccumulator {
        let f = make_bed_file(bed_content);
        let parsed = ParsedBed::from_file(f.path().to_str().unwrap(), None, false);
        BedcovAccumulator::new(parsed, templates)
    }

    /// Feature 2: basic take_label_coverages with MeanGenomeCoverageEstimator.
    /// BED: two labels, each one region.  Coverage is known from ups_and_downs.
    #[test]
    fn test_feature2_take_label_coverages_mean() {
        // label_a: [0,5) — cov = [2,2,2,2,2], mean = 2.0
        // label_b: [5,10) — cov = [0,0,0,0,0], mean = 0.0
        let mut acc = make_accumulator_with_estimators(
            "chr1\t0\t5\tlabel_a\nchr1\t5\t10\tlabel_b\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
        fake_index_tid0(&mut acc);

        let ups: Vec<i32> = vec![2, 0, 0, 0, 0, -2, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        // covs[label_idx][method_idx]
        // labels sorted alphabetically: label_a=0, label_b=1
        assert_eq!(covs.len(), 2, "expected 2 labels");
        assert_eq!(covs[0].len(), 1, "expected 1 method");
        assert!(
            (covs[0][0] - 2.0).abs() < 1e-5,
            "label_a mean: expected 2.0, got {}",
            covs[0][0]
        );
        assert!(
            covs[1][0].abs() < 1e-5,
            "label_b mean: expected 0.0, got {}",
            covs[1][0]
        );
    }

    /// Feature 2: covered_fraction for a partial-coverage region.
    #[test]
    fn test_feature2_covered_fraction() {
        // label_a: [0,10) — cov = [1,1,1,1,1,0,0,0,0,0], covered 5/10 = 0.5
        let mut acc = make_accumulator_with_estimators(
            "chr1\t0\t10\tlabel_a\n",
            vec![CoverageEstimator::new_estimator_covered_fraction(0.0)],
        );
        fake_index_tid0(&mut acc);

        let ups: Vec<i32> = vec![1, 0, 0, 0, 0, -1, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        assert!(
            (covs[0][0] - 0.5).abs() < 1e-5,
            "covered_fraction: expected 0.5, got {}",
            covs[0][0]
        );
    }

    /// Feature 2: two methods, two labels — verify matrix dimensions and values.
    #[test]
    fn test_feature2_two_methods_two_labels() {
        // label_a: [0,4) cov=[3,3,3,3]  mean=3.0, fraction=1.0
        // label_b: [4,8) cov=[0,0,0,0]  mean=0.0, fraction=0.0
        let mut acc = make_accumulator_with_estimators(
            "chr1\t0\t4\tlabel_a\nchr1\t4\t8\tlabel_b\n",
            vec![
                CoverageEstimator::new_estimator_mean(0.0, 0, false),
                CoverageEstimator::new_estimator_covered_fraction(0.0),
            ],
        );
        fake_index_tid0(&mut acc);

        let ups: Vec<i32> = vec![3, 0, 0, 0, -3, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        assert_eq!(covs.len(), 2, "2 labels");
        assert_eq!(covs[0].len(), 2, "2 methods");
        // label_a (idx 0): mean=3.0, frac=1.0
        assert!((covs[0][0] - 3.0).abs() < 1e-5, "label_a mean");
        assert!((covs[0][1] - 1.0).abs() < 1e-5, "label_a fraction");
        // label_b (idx 1): mean=0.0, frac=0.0
        assert!(covs[1][0].abs() < 1e-5, "label_b mean");
        assert!(covs[1][1].abs() < 1e-5, "label_b fraction");
    }

    /// Feature 2: multiple contigs same genome accumulate correctly.
    #[test]
    fn test_feature2_multi_contig_accumulation() {
        // label_a spans two contigs (tid=0 and tid=1).
        // After two process_contig calls, take_label_coverages returns the sum.
        // tid=0: [0,5)  ups=[1,0,0,0,-1,…] cov=[1,1,1,1,0], sum_cov=4
        // tid=1: [0,5)  ups=[2,0,0,0,-2,…] cov=[2,2,2,2,0], sum_cov=8
        // Total observed bases = 10, total sum = 12 → mean = 1.2
        let f = make_bed_file("chr1\t0\t5\tlabel_a\nchr2\t0\t5\tlabel_a\n");
        let parsed = ParsedBed::from_file(f.path().to_str().unwrap(), None, false);
        let mut acc = BedcovAccumulator::new(
            parsed,
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );

        // Build a two-tid index manually
        let mut regions_by_tid: HashMap<u32, Vec<IndexedRegion>> = HashMap::new();
        let bed = &acc.parsed_bed;
        // tid=0 → chr1 regions, tid=1 → chr2 regions
        for (chrom, tid) in [("chr1", 0u32), ("chr2", 1u32)] {
            if let Some(rs) = bed.regions_by_chrom.get(chrom) {
                regions_by_tid.insert(
                    tid,
                    rs.iter()
                        .map(|r| IndexedRegion {
                            start: r.start,
                            end: r.end,
                            label_idx: r.label_idx,
                            global_idx: r.global_idx,
                        })
                        .collect(),
                );
            }
        }
        acc.current_index = BedIndex {
            regions_by_tid,
            tid_to_chrom_idx: HashMap::new(),
        };

        let ups0: Vec<i32> = vec![1, 0, 0, 0, -1, 0, 0, 0, 0, 0];
        let ups1: Vec<i32> = vec![2, 0, 0, 0, -2, 0, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups0);
        acc.process_contig(1, "g1", &ups1);

        let covs = acc.take_label_coverages("g1");
        // label_a mean across both regions: (4+8) / 10 = 1.2
        assert!(
            (covs[0][0] - 1.2).abs() < 1e-5,
            "multi-contig mean: expected 1.2, got {}",
            covs[0][0]
        );
    }

    /// Feature 2: region of length 1 bp — no panic, correct value.
    #[test]
    fn test_feature2_region_length_1() {
        // Region [3,4) — single base.  cov at pos 3 = 5 → mean = 5.0
        let mut acc = make_accumulator_with_estimators(
            "chr1\t3\t4\tlabel_a\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
        fake_index_tid0(&mut acc);

        // ups_and_downs: cov = [0,0,0,5,0,0,0,0,0,0]
        let ups: Vec<i32> = vec![0, 0, 0, 5, -5, 0, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        assert!(
            (covs[0][0] - 5.0).abs() < 1e-5,
            "1-bp region mean: expected 5.0, got {}",
            covs[0][0]
        );
    }

    /// Feature 2: genome with no BED regions on any of its contigs returns zeros
    /// (from the None branch — no entry ever created in genome_label_estimators).
    #[test]
    fn test_feature2_genome_no_regions_returns_zeros() {
        let mut acc = make_accumulator_with_estimators(
            "chr1\t0\t5\tlabel_a\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
        // No call to fake_index_tid0 → current_index is empty → process_contig returns early.
        // We call take_label_coverages for a genome that was never seen.
        let covs = acc.take_label_coverages("never_seen_genome");
        assert_eq!(covs.len(), 1, "1 label");
        assert_eq!(covs[0].len(), 1, "1 method");
        assert!(
            covs[0][0].abs() < 1e-5,
            "fresh estimator should return 0.0, got {}",
            covs[0][0]
        );
        assert!(covs[0][0].is_finite(), "must not be NaN or infinite");
    }

    /// Feature 2: label with no BED regions on any processed contig returns 0.0
    /// (entry created for genome via label_a hit, but label_b estimators are fresh).
    #[test]
    fn test_feature2_label_no_regions_returns_zero() {
        // label_a has a region on chr1 (tid=0), label_b only on chr2 (not indexed).
        let mut acc = make_accumulator_with_estimators(
            "chr1\t0\t5\tlabel_a\nchr2\t0\t5\tlabel_b\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );

        // Only map chr1 to tid=0; chr2 gets no tid → label_b has no indexed regions.
        let mut regions_by_tid: HashMap<u32, Vec<IndexedRegion>> = HashMap::new();
        if let Some(rs) = acc.parsed_bed.regions_by_chrom.get("chr1") {
            regions_by_tid.insert(
                0,
                rs.iter()
                    .map(|r| IndexedRegion {
                        start: r.start,
                        end: r.end,
                        label_idx: r.label_idx,
                        global_idx: r.global_idx,
                    })
                    .collect(),
            );
        }
        acc.current_index = BedIndex {
            regions_by_tid,
            tid_to_chrom_idx: HashMap::new(),
        };

        let ups: Vec<i32> = vec![3, 0, 0, 0, 0, -3, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        // labels sorted alphabetically: label_a=0, label_b=1
        assert!(
            (covs[0][0] - 3.0).abs() < 1e-5,
            "label_a mean: expected 3.0, got {}",
            covs[0][0]
        );
        assert!(
            covs[1][0].abs() < 1e-5,
            "label_b (no regions on tid=0): expected 0.0, got {}",
            covs[1][0]
        );
        assert!(covs[1][0].is_finite(), "label_b must not be NaN");
    }

    /// Feature 2: take_label_coverages resets state — second call returns zeros.
    #[test]
    fn test_feature2_take_is_consume() {
        let mut acc = make_accumulator_with_estimators(
            "chr1\t0\t5\tlabel_a\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
        fake_index_tid0(&mut acc);

        let ups: Vec<i32> = vec![1, 0, 0, 0, 0, -1, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs1 = acc.take_label_coverages("g1");
        assert!((covs1[0][0] - 1.0).abs() < 1e-5, "first call");

        // Second call — entry was removed, should return zeros
        let covs2 = acc.take_label_coverages("g1");
        assert!(covs2[0][0].abs() < 1e-5, "second call should be 0.0");
    }

    /// Gate test for Phase 2 / Feature 2.
    ///
    /// Validates the sub_ups reconstruction technique: for a region [s, e), we
    /// build a synthetic ups_and_downs slice where
    ///   sub_ups[0] = pcov[s+1] - pcov[s]   (absolute coverage at position s)
    ///   sub_ups[i] = ups_and_downs[s+i]      for i ≥ 1
    ///
    /// Claim: prefix_sum(sub_ups)[j] == coverage[s+j] for all j in 0..len.
    ///
    /// If this holds, calling CoverageEstimator::add_contig(&sub_ups, ...) with
    /// contig_end_exclusion=0 gives the same result as processing the original
    /// ups_and_downs restricted to [s, e).
    ///
    /// The test verifies MeanGenomeCoverageEstimator and
    /// CoverageFractionGenomeCoverageEstimator across three non-trivial
    /// sub-regions of an 11-bp contig.
    #[test]
    fn test_sub_ups_correctness() {
        // Contig of 11 bp:
        //   ups = [2, -1, 0, 0, -1, 0, 3, -3, 0, 1, -1]
        //   cov  = [2,  1,  1,  1,  0, 0,  3,  0, 0,  1,  0]
        //   pcov = [0,  2,  3,  4,  5, 5,  5,  8, 8,  8,  9, 9]
        let ups: Vec<i32> = vec![2, -1, 0, 0, -1, 0, 3, -3, 0, 1, -1];
        let n = ups.len();

        // Build prefix-coverage array: pcov[i] = sum of coverage[0..i]
        let mut pcov = vec![0i64; n + 1];
        let mut running: i64 = 0;
        for (i, &delta) in ups.iter().enumerate() {
            running += delta as i64;
            pcov[i + 1] = pcov[i] + running;
        }
        assert_eq!(
            pcov,
            vec![0, 2, 3, 4, 5, 5, 5, 8, 8, 8, 9, 9],
            "pcov sanity check"
        );

        // Helper: build sub_ups for region [s, e)
        let build_sub_ups = |s: usize, e: usize| -> Vec<i32> {
            let len = e - s;
            let mut v = vec![0i32; len];
            v[0] = (pcov[s + 1] - pcov[s]) as i32; // absolute cov at position s
            for i in 1..len {
                v[i] = ups[s + i];
            }
            v
        };

        // -------------------------------------------------------------------
        // Region [0, 5): cov = [2,1,1,1,0], sum=5, mean=1.0, covered=4/5=0.8
        // -------------------------------------------------------------------
        let sub05 = build_sub_ups(0, 5);
        assert_eq!(sub05, vec![2, -1, 0, 0, -1], "sub_ups [0,5)");

        let mut est = CoverageEstimator::new_estimator_mean(0.0, 0, false);
        est.add_contig(&sub05, 0, 0, 0.0);
        let mean = est.calculate_coverage(&[]);
        assert!(
            (mean - 1.0).abs() < 1e-6,
            "Mean [0,5): expected 1.0, got {}",
            mean
        );

        let mut est = CoverageEstimator::new_estimator_covered_fraction(0.0);
        est.add_contig(&sub05, 0, 0, 0.0);
        let frac = est.calculate_coverage(&[]);
        assert!(
            (frac - 0.8).abs() < 1e-6,
            "CovFrac [0,5): expected 0.8, got {}",
            frac
        );

        // -------------------------------------------------------------------
        // Region [4, 8): cov = [0,0,3,0], sum=3, mean=0.75, covered=1/4=0.25
        // -------------------------------------------------------------------
        let sub48 = build_sub_ups(4, 8);
        assert_eq!(sub48, vec![0, 0, 3, -3], "sub_ups [4,8)");

        let mut est = CoverageEstimator::new_estimator_mean(0.0, 0, false);
        est.add_contig(&sub48, 0, 0, 0.0);
        let mean = est.calculate_coverage(&[]);
        assert!(
            (mean - 0.75).abs() < 1e-6,
            "Mean [4,8): expected 0.75, got {}",
            mean
        );

        let mut est = CoverageEstimator::new_estimator_covered_fraction(0.0);
        est.add_contig(&sub48, 0, 0, 0.0);
        let frac = est.calculate_coverage(&[]);
        assert!(
            (frac - 0.25).abs() < 1e-6,
            "CovFrac [4,8): expected 0.25, got {}",
            frac
        );

        // -------------------------------------------------------------------
        // Region [6, 11): cov = [3,0,0,1,0], sum=4, mean=0.8, covered=2/5=0.4
        // -------------------------------------------------------------------
        let sub611 = build_sub_ups(6, 11);
        assert_eq!(sub611, vec![3, -3, 0, 1, -1], "sub_ups [6,11)");

        let mut est = CoverageEstimator::new_estimator_mean(0.0, 0, false);
        est.add_contig(&sub611, 0, 0, 0.0);
        let mean = est.calculate_coverage(&[]);
        assert!(
            (mean - 0.8).abs() < 1e-6,
            "Mean [6,11): expected 0.8, got {}",
            mean
        );

        let mut est = CoverageEstimator::new_estimator_covered_fraction(0.0);
        est.add_contig(&sub611, 0, 0, 0.0);
        let frac = est.calculate_coverage(&[]);
        assert!(
            (frac - 0.4).abs() < 1e-6,
            "CovFrac [6,11): expected 0.4, got {}",
            frac
        );
    }

    // -----------------------------------------------------------------------
    // Unlabeled-gaps tests
    // -----------------------------------------------------------------------

    /// Basic unlabeled: single labeled region [2,7) on a 10-bp contig.
    /// Gaps are [0,2) and [7,10).
    /// Coverage: cov=[0,0,1,1,1,1,0,0,0,0]
    /// label_a mean = (1+1+1+1+0)/5 = 0.8
    /// unlabeled mean over [0,2)∪[7,10):
    ///   [0,2): sum=0, len=2
    ///   [7,10): sum=0, len=3
    ///   total_sum=0, total_len=5 → mean=0.0
    #[test]
    fn test_unlabeled_basic() {
        let mut acc = make_accumulator_unlabeled(
            "chr1\t2\t7\tlabel_a\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            "unlabeled",
        );
        fake_index_tid0(&mut acc);

        // labels: ["label_a", "unlabeled"]
        assert_eq!(acc.parsed_bed.labels, vec!["label_a", "unlabeled"]);
        assert_eq!(acc.parsed_bed.unlabeled_label_idx, Some(1));

        let ups: Vec<i32> = vec![0, 0, 1, 0, 0, 0, -1, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        assert_eq!(covs.len(), 2);
        assert!(
            (covs[0][0] - 0.8).abs() < 1e-5,
            "label_a mean: expected 0.8, got {}",
            covs[0][0]
        );
        assert!(
            covs[1][0].abs() < 1e-5,
            "unlabeled mean: expected 0.0, got {}",
            covs[1][0]
        );
    }

    /// Unlabeled with coverage in the gaps.
    /// Contig 10 bp, coverage = [3,3,3,3,3,3,3,3,3,3] (uniform 3, no drop at end).
    /// Region [3,7) → label_a mean = 3.0.
    /// Gaps [0,3) + [7,10), each bp has cov=3 → unlabeled mean = 3.0.
    #[test]
    fn test_unlabeled_with_gap_coverage() {
        let mut acc = make_accumulator_unlabeled(
            "chr1\t3\t7\tlabel_a\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            "unlabeled",
        );
        fake_index_tid0(&mut acc);

        // No drop at the end so all 10 positions have cov=3.
        let ups: Vec<i32> = vec![3, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        assert!(
            (covs[0][0] - 3.0).abs() < 1e-5,
            "label_a mean: expected 3.0, got {}",
            covs[0][0]
        );
        assert!(
            (covs[1][0] - 3.0).abs() < 1e-5,
            "unlabeled mean: expected 3.0, got {}",
            covs[1][0]
        );
    }

    /// Unlabeled with cross-label overlapping regions.
    /// Contig 10 bp, cov = [1,1,1,1,1,1,1,1,1,1] (uniform 1, no drop at end).
    /// label_a [1,5), label_b [3,8) — overlap at [3,5) is cross-label.
    /// Merged union: [1,8) → gaps: [0,1) and [8,10).
    /// label_a mean over [1,5): 1.0
    /// label_b mean over [3,8): 1.0
    /// unlabeled mean over [0,1)∪[8,10): 3 bp with cov=1 → mean=1.0
    #[test]
    fn test_unlabeled_cross_label_overlap() {
        let mut acc = make_accumulator_unlabeled(
            "chr1\t1\t5\tlabel_a\nchr1\t3\t8\tlabel_b\n",
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
            "unlabeled",
        );
        fake_index_tid0(&mut acc);

        // No drop so all 10 positions have cov=1.
        let ups: Vec<i32> = vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        // labels sorted alphabetically: label_a=0, label_b=1, unlabeled=2
        assert_eq!(covs.len(), 3, "expected 3 labels (a, b, unlabeled)");
        assert!(
            (covs[0][0] - 1.0).abs() < 1e-5,
            "label_a mean: expected 1.0, got {}",
            covs[0][0]
        );
        assert!(
            (covs[1][0] - 1.0).abs() < 1e-5,
            "label_b mean: expected 1.0, got {}",
            covs[1][0]
        );
        assert!(
            (covs[2][0] - 1.0).abs() < 1e-5,
            "unlabeled mean: expected 1.0, got {}",
            covs[2][0]
        );
    }

    /// Contig with NO BED regions at all — the entire contig is unlabeled.
    /// Contig 5 bp, cov = [2,2,2,2,2].
    /// process_contig is called; the entire contig is fed to the unlabeled estimator.
    #[test]
    fn test_unlabeled_no_bed_regions_on_contig() {
        // BED has a region on chr2 but we process chr1 (tid=0 mapped to empty slot).
        let f = make_bed_file("chr2\t0\t5\tlabel_a\n");
        let parsed = ParsedBed::from_file(f.path().to_str().unwrap(), Some("unlabeled"), false);
        let mut acc = BedcovAccumulator::new(
            parsed,
            vec![CoverageEstimator::new_estimator_mean(0.0, 0, false)],
        );
        // tid=0 has NO entry in current_index (chr1 not in BED) — leave index empty.
        // process_contig should still run and populate the unlabeled estimator.
        // No drop so all 5 positions have cov=2.
        let ups: Vec<i32> = vec![2, 0, 0, 0, 0];
        acc.process_contig(0, "g1", &ups);

        let covs = acc.take_label_coverages("g1");
        assert_eq!(covs.len(), 2, "expected 2 labels (label_a, unlabeled)");
        assert!(
            covs[0][0].abs() < 1e-5,
            "label_a mean: expected 0.0 (no regions on this contig), got {}",
            covs[0][0]
        );
        assert!(
            (covs[1][0] - 2.0).abs() < 1e-5,
            "unlabeled mean: expected 2.0, got {}",
            covs[1][0]
        );
    }

    /// Confirm that need_bedgraph=false suppresses all three optional allocations.
    #[test]
    fn test_no_alloc_without_bedgraph() {
        let f = make_bed_file("chr1\t0\t10\tlabel_a\n");
        let parsed = ParsedBed::from_file(f.path().to_str().unwrap(), None, false);
        assert!(
            parsed.all_regions.is_none(),
            "all_regions should be None when need_bedgraph=false"
        );
        assert!(
            parsed.region_lengths.is_none(),
            "region_lengths should be None when need_bedgraph=false"
        );
        let acc = BedcovAccumulator::new(parsed, vec![]);
        assert!(
            acc.region_sums.is_none(),
            "region_sums should be None when need_bedgraph=false"
        );
    }
}
