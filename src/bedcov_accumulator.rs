use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam;

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
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub label_idx: usize,
}

/// Parsed BED file: labels (alphabetical), regions indexed two ways.
pub struct ParsedBed {
    /// chrom → regions sorted by start
    pub regions_by_chrom: HashMap<String, Vec<BedRegion>>,
    /// alphabetical list of distinct labels
    pub labels: Vec<String>,
    pub label_to_idx: HashMap<String, usize>,
    /// all_regions[global_idx] — original BED order
    pub all_regions: Vec<RegionInfo>,
    /// length of each region in bp, parallel to all_regions
    pub region_lengths: Vec<u32>,
}

impl ParsedBed {
    pub fn from_file(path: &str) -> ParsedBed {
        let file = fs::File::open(path)
            .unwrap_or_else(|e| panic!("Cannot open --regions-bed file '{}': {}", path, e));
        let reader = BufReader::new(file);

        // --- first pass: collect raw records ---
        struct RawRecord {
            chrom: String,
            start: u32,
            end: u32,
            label: String,
        }
        let mut raw: Vec<RawRecord> = Vec::new();
        let mut label_set: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();

        for (lineno, line_res) in reader.lines().enumerate() {
            let line = line_res.unwrap_or_else(|e| {
                panic!(
                    "Error reading --regions-bed '{}' line {}: {}",
                    path,
                    lineno + 1,
                    e
                )
            });
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
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
            label_set.insert(label.clone());
            raw.push(RawRecord {
                chrom,
                start,
                end,
                label,
            });
        }

        // --- build label index (alphabetical) ---
        let labels: Vec<String> = label_set.into_iter().collect();
        let label_to_idx: HashMap<String, usize> = labels
            .iter()
            .enumerate()
            .map(|(i, l)| (l.clone(), i))
            .collect();

        // --- build all_regions, regions_by_chrom ---
        let mut all_regions: Vec<RegionInfo> = Vec::with_capacity(raw.len());
        let mut region_lengths: Vec<u32> = Vec::with_capacity(raw.len());
        let mut regions_by_chrom: HashMap<String, Vec<BedRegion>> = HashMap::new();

        for r in raw {
            let global_idx = all_regions.len();
            let label_idx = label_to_idx[&r.label];
            all_regions.push(RegionInfo {
                chrom: r.chrom.clone(),
                start: r.start,
                end: r.end,
                label_idx,
            });
            region_lengths.push(r.end - r.start);
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

        // warn on overlapping regions within the same chrom
        for (chrom, v) in &regions_by_chrom {
            for w in v.windows(2) {
                if w[1].start < w[0].end {
                    warn!(
                        "Overlapping BED regions detected on contig '{}' \
                         ([{},{}] and [{},{}]); each region is counted independently",
                        chrom, w[0].start, w[0].end, w[1].start, w[1].end
                    );
                }
            }
        }

        info!(
            "Loaded {} regions across {} distinct label(s) from '{}'",
            all_regions.len(),
            labels.len(),
            path
        );

        ParsedBed {
            regions_by_chrom,
            labels,
            label_to_idx,
            all_regions,
            region_lengths,
        }
    }
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
}

// ---------------------------------------------------------------------------
// BedcovAccumulator — main state
// ---------------------------------------------------------------------------

pub struct BedcovAccumulator {
    pub parsed_bed: ParsedBed,
    pub current_index: BedIndex,
    /// sum of coverage per region (reset at start of each BAM)
    pub region_sums: Vec<i64>,
    /// prefix-sum scratch (reused, never shrunk)
    pcov_scratch: Vec<i64>,
}

impl BedcovAccumulator {
    pub fn new(parsed_bed: ParsedBed) -> BedcovAccumulator {
        let n = parsed_bed.all_regions.len();
        BedcovAccumulator {
            parsed_bed,
            current_index: BedIndex {
                regions_by_tid: HashMap::new(),
            },
            region_sums: vec![0i64; n],
            pcov_scratch: Vec::new(),
        }
    }

    /// Rebuild BedIndex from the BAM header (call once per BAM, before processing reads).
    pub fn reinit_for_bam(&mut self, header: &bam::HeaderView) {
        self.region_sums.fill(0);
        let mut regions_by_tid: HashMap<u32, Vec<IndexedRegion>> = HashMap::new();

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
            }
        }
        self.current_index = BedIndex { regions_by_tid };
    }

    /// Accumulate coverage sums for all BED regions on this contig.
    ///
    /// `ups_and_downs` is the mosdepth difference array: coverage at position i =
    /// prefix_sum(ups_and_downs[0..=i]).
    pub fn process_contig(&mut self, tid: u32, ups_and_downs: &[i32]) {
        let regions = match self.current_index.regions_by_tid.get(&tid) {
            Some(r) => r,
            None => return,
        };

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

        for region in regions {
            let s = (region.start as usize).min(n);
            let e = (region.end as usize).min(n);
            if s >= e {
                continue;
            }
            self.region_sums[region.global_idx] += self.pcov_scratch[e] - self.pcov_scratch[s];
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

        for (label_idx, label) in bed.labels.iter().enumerate() {
            writeln!(w, "track type=bedGraph name=\"{}\"", label)
                .expect("Error writing bedGraph track header");

            for (global_idx, region) in bed.all_regions.iter().enumerate() {
                if region.label_idx != label_idx {
                    continue;
                }
                let len = bed.region_lengths[global_idx] as i64;
                let mean_cov = if len > 0 {
                    self.region_sums[global_idx] as f64 / len as f64
                } else {
                    0.0
                };
                writeln!(
                    w,
                    "{}\t{}\t{}\t{:.6}",
                    region.chrom, region.start, region.end, mean_cov
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
        let parsed = ParsedBed::from_file(f.path().to_str().unwrap());
        BedcovAccumulator::new(parsed)
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
        acc.current_index = BedIndex { regions_by_tid };
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
        acc.process_contig(0, &ups);
        assert_eq!(acc.region_sums[0], 4);
    }

    #[test]
    fn test_region_clamp() {
        // Region [8, 15) on a 10-bp contig — should clamp to [8, 10), no panic.
        let mut acc = make_accumulator("chr1\t8\t15\tlabel_a\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![0, 0, 1, 0, 0, 0, -1, 0, 0, 0];
        acc.process_contig(0, &ups);
        // coverage at [8,10) = 0, so sum = 0
        assert_eq!(acc.region_sums[0], 0);
    }

    #[test]
    fn test_zero_coverage() {
        let mut acc = make_accumulator("chr1\t0\t10\tlabel_a\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![0; 10];
        acc.process_contig(0, &ups);
        assert_eq!(acc.region_sums[0], 0);
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
        acc.process_contig(0, &ups);
        // region0=[0,3): sum=3, region1=[5,8): sum=6
        assert_eq!(acc.region_sums[0], 3);
        assert_eq!(acc.region_sums[1], 6);
    }

    #[test]
    fn test_overlapping_regions() {
        // [0,5) and [3,8) — overlap, each counts independently
        let mut acc = make_accumulator("chr1\t0\t5\tlabel_a\nchr1\t3\t8\tlabel_b\n");
        fake_index_tid0(&mut acc);
        // uniform coverage = 1 everywhere
        let ups: Vec<i32> = vec![1, 0, 0, 0, 0, 0, 0, 0, -1, 0];
        acc.process_contig(0, &ups);
        assert_eq!(acc.region_sums[0], 5); // [0,5)
        assert_eq!(acc.region_sums[1], 5); // [3,8)
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
        acc.process_contig(0, &ups);
        assert_eq!(acc.region_sums[0], 5); // label_a
        assert_eq!(acc.region_sums[1], 0); // label_b (no coverage)
    }

    #[test]
    fn test_finalise_bedgraph_content() {
        let mut acc = make_accumulator("chr1\t0\t5\tlabel_a\nchr1\t5\t10\tlabel_b\n");
        fake_index_tid0(&mut acc);
        let ups: Vec<i32> = vec![2, 0, 0, 0, 0, -2, 0, 0, 0, 0];
        acc.process_contig(0, &ups);

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
}
