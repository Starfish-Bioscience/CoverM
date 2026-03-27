use coverage_takers::CoverageTaker;

#[derive(Clone, Debug)]
pub enum CoverageEstimator {
    MeanGenomeCoverageEstimator {
        total_count: u64,
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        total_mismatches: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        exclude_mismatches: bool,
    },
    TrimmedMeanGenomeCoverageEstimator {
        counts: Vec<u64>,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min: f32,
        max: f32,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    },
    PileupCountsGenomeCoverageEstimator {
        counts: Vec<u64>,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    },
    CoverageFractionGenomeCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    NumCoveredBasesCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    RPKMCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    TPMCoverageEstimator {
        total_bases: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
    },
    VarianceGenomeCoverageEstimator {
        counts: Vec<u64>,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    },
    ReferenceLengthCalculator {
        observed_contig_length: u64,
        num_mapped_reads: u64,
    },
    ReadCountCalculator {
        num_mapped_reads: u64,
    },
    ReadsPerBaseCalculator {
        observed_contig_length: u64,
        num_mapped_reads: u64,
    },
    AverageIdentityEstimator {
        sum_identity: f64,
        num_reads: u64,
    },
    StrobealignAembEstimator {},
    IslandsPerMbpEstimator {
        // Metric accumulators
        total_islands: u64,
        covered_contigs_length: u64, // denominator: contigs with ≥1 valid island
        // Gate accumulators (raw coverage, before min_island_length filtering)
        observed_contig_length: u64, // all contigs seen (with or without signal)
        num_covered_bases: u64,      // raw bases at depth>0
        num_mapped_reads: u64,
        // Config (immutable after construction)
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        min_island_length: u64,
    },
    MaxGapEstimator {
        max_gap: u64,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        min_island_length: u64,
    },
    GapFractionEstimator {
        total_internal_gap_bases: u64,
        total_internal_span: u64,
        observed_contig_length: u64,
        num_covered_bases: u64,
        num_mapped_reads: u64,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        min_island_length: u64,
    },
}

impl CoverageEstimator {
    pub fn column_headers(&self) -> Vec<&str> {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator { .. } => vec!["Mean"],
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator { .. } => vec!["Trimmed Mean"],
            CoverageEstimator::PileupCountsGenomeCoverageEstimator { .. } => {
                vec!["Coverage", "Bases"]
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator { .. } => {
                vec!["Covered Fraction"]
            }
            CoverageEstimator::NumCoveredBasesCoverageEstimator { .. } => vec!["Covered Bases"],
            CoverageEstimator::RPKMCoverageEstimator { .. } => vec!["RPKM"],
            CoverageEstimator::TPMCoverageEstimator { .. } => vec!["TPM"],
            CoverageEstimator::VarianceGenomeCoverageEstimator { .. } => vec!["Variance"],
            CoverageEstimator::ReferenceLengthCalculator { .. } => vec!["Length"],
            CoverageEstimator::ReadCountCalculator { .. } => vec!["Read Count"],
            CoverageEstimator::ReadsPerBaseCalculator { .. } => vec!["Reads per base"],
            CoverageEstimator::AverageIdentityEstimator { .. } => vec!["ANIr"],
            CoverageEstimator::StrobealignAembEstimator { .. } => vec!["Strobealign aemb"],
            CoverageEstimator::IslandsPerMbpEstimator { .. } => vec!["Islands per Mbp"],
            CoverageEstimator::MaxGapEstimator { .. } => vec!["Max Gap"],
            CoverageEstimator::GapFractionEstimator { .. } => vec!["Gap Fraction"],
        }
    }
}

impl CoverageEstimator {
    pub fn new_estimator_mean(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        exclude_mismatches: bool,
    ) -> CoverageEstimator {
        CoverageEstimator::MeanGenomeCoverageEstimator {
            total_count: 0,
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            total_mismatches: 0,
            min_fraction_covered_bases,
            contig_end_exclusion,
            exclude_mismatches,
        }
    }
    pub fn new_estimator_trimmed_mean(
        min: f32,
        max: f32,
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
            counts: vec![],
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
            min,
            max,
            contig_end_exclusion,
        }
    }
    pub fn new_estimator_pileup_counts(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::PileupCountsGenomeCoverageEstimator {
            counts: vec![],
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
            contig_end_exclusion,
        }
    }
    pub fn new_estimator_covered_fraction(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_rpkm(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::RPKMCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_tpm(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::TPMCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_covered_bases(min_fraction_covered_bases: f32) -> CoverageEstimator {
        CoverageEstimator::NumCoveredBasesCoverageEstimator {
            total_bases: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
        }
    }
    pub fn new_estimator_variance(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::VarianceGenomeCoverageEstimator {
            counts: vec![],
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
            contig_end_exclusion,
        }
    }
    pub fn new_estimator_length() -> CoverageEstimator {
        CoverageEstimator::ReferenceLengthCalculator {
            observed_contig_length: 0,
            num_mapped_reads: 0,
        }
    }
    pub fn new_estimator_read_count() -> CoverageEstimator {
        CoverageEstimator::ReadCountCalculator {
            num_mapped_reads: 0,
        }
    }
    pub fn new_estimator_reads_per_base() -> CoverageEstimator {
        CoverageEstimator::ReadsPerBaseCalculator {
            observed_contig_length: 0,
            num_mapped_reads: 0,
        }
    }
    pub fn new_estimator_anir() -> CoverageEstimator {
        CoverageEstimator::AverageIdentityEstimator {
            sum_identity: 0.0,
            num_reads: 0,
        }
    }
    pub fn new_estimator_strobealign_aemb() -> CoverageEstimator {
        CoverageEstimator::StrobealignAembEstimator {}
    }

    pub fn new_estimator_islands_per_mbp(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        min_island_length: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::IslandsPerMbpEstimator {
            total_islands: 0,
            covered_contigs_length: 0,
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
            contig_end_exclusion,
            min_island_length,
        }
    }

    pub fn new_estimator_max_gap(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        min_island_length: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::MaxGapEstimator {
            max_gap: 0,
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
            contig_end_exclusion,
            min_island_length,
        }
    }

    pub fn new_estimator_gap_fraction(
        min_fraction_covered_bases: f32,
        contig_end_exclusion: u64,
        min_island_length: u64,
    ) -> CoverageEstimator {
        CoverageEstimator::GapFractionEstimator {
            total_internal_gap_bases: 0,
            total_internal_span: 0,
            observed_contig_length: 0,
            num_covered_bases: 0,
            num_mapped_reads: 0,
            min_fraction_covered_bases,
            contig_end_exclusion,
            min_island_length,
        }
    }

    fn calculate_unobserved_bases(
        unobserved_contig_lengths: &[u64],
        contig_end_exclusion: u64,
    ) -> u64 {
        let unobserved_not_excluded = unobserved_contig_lengths
            .iter()
            .map(|l| {
                let e = &(2 * contig_end_exclusion);
                if l < e {
                    *l
                } else {
                    l - e
                }
            })
            .sum();
        unobserved_not_excluded
    }
}

/// Result of a spatial scan of one contig's coverage profile.
#[derive(Debug, Clone)]
pub struct SpatialScanResult {
    /// Number of valid islands (covered segments ≥ min_island_length)
    pub n_islands: u64,
    /// Total bases at depth 0 between two valid islands (not prefixes/suffixes)
    pub total_internal_gap_bases: u64,
    /// Largest internal gap in this contig
    pub max_internal_gap: u64,
    /// Contig length after contig_end_exclusion
    pub analysed_contig_length: u64,
    /// Positions between first and last covered base, inclusive (last - first + 1).
    /// This is the region where internal gaps can exist.
    pub internal_span: u64,
    /// Raw number of covered bases (depth > 0), before min_island_length filtering.
    /// Used for the min_fraction_covered_bases gate.
    pub num_covered_bases: u64,
}

/// Scan a single contig's pileup for spatial coverage structure.
///
/// Streaming algorithm: O(contig_length) time, O(1) extra memory (no per-base vector).
///
/// Returns None if:
/// - contig too short (< 2 × contig_end_exclusion)
/// - no base has depth > 0
/// - all covered segments are shorter than min_island_length (no valid island)
pub fn spatial_scan(
    ups_and_downs: &[i32],
    contig_end_exclusion: u64,
    min_island_length: u64,
) -> Option<SpatialScanResult> {
    let len = ups_and_downs.len();
    if contig_end_exclusion * 2 >= len as u64 {
        return None; // contig too short
    }

    let start_from = contig_end_exclusion as usize;
    let end_at = len - contig_end_exclusion as usize - 1;
    let analysed_contig_length = (end_at - start_from + 1) as u64;

    // State variables for streaming run detection
    let mut cumulative_sum: i32 = 0;
    let mut num_covered_bases: u64 = 0;

    // Island/gap tracking
    let mut n_islands: u64 = 0;
    let mut total_internal_gap_bases: u64 = 0;
    let mut max_internal_gap: u64 = 0;
    let mut first_island_seen = false;
    let mut first_covered_pos: u64 = 0;
    let mut last_covered_pos: u64 = 0;

    // Current run tracking
    let mut in_covered_run = false;
    let mut current_run_length: u64 = 0;

    // pending_gap_length accumulates gap bases + reclassified micro-islands,
    // waiting to be confirmed as an internal gap when the next valid island is found.
    let mut pending_gap_length: u64 = 0;

    for (i, current) in ups_and_downs.iter().enumerate() {
        cumulative_sum += current;

        if i < start_from || i > end_at {
            continue;
        }

        let is_covered = cumulative_sum > 0;

        if is_covered {
            num_covered_bases += 1;
        }

        if is_covered && !in_covered_run {
            // Transition: uncovered → covered — start a new covered run
            if in_covered_run || current_run_length > 0 {
                // We were in an uncovered run — store its length in pending
                pending_gap_length += current_run_length;
            }
            in_covered_run = true;
            current_run_length = 1;
        } else if !is_covered && in_covered_run {
            // Transition: covered → uncovered — end of a covered run
            if current_run_length < min_island_length {
                // Micro-island: reclassify as uncovered
                pending_gap_length += current_run_length;
            } else {
                // Valid island
                if first_island_seen {
                    // pending_gap_length is now a confirmed internal gap
                    total_internal_gap_bases += pending_gap_length;
                    if pending_gap_length > max_internal_gap {
                        max_internal_gap = pending_gap_length;
                    }
                } else {
                    first_covered_pos = (i as u64) - current_run_length;
                    first_island_seen = true;
                }
                last_covered_pos = (i as u64) - 1; // last position of the island
                n_islands += 1;
                pending_gap_length = 0;
            }
            in_covered_run = false;
            current_run_length = 1;
        } else {
            current_run_length += 1;
        }
    }

    // End of contig: close the last run
    if in_covered_run && current_run_length >= min_island_length {
        // Last run is a valid island
        if first_island_seen {
            total_internal_gap_bases += pending_gap_length;
            if pending_gap_length > max_internal_gap {
                max_internal_gap = pending_gap_length;
            }
        } else {
            first_covered_pos = (end_at as u64 + 1) - current_run_length;
            first_island_seen = true;
        }
        last_covered_pos = end_at as u64;
        n_islands += 1;
    }
    // If last run is a micro-island or uncovered, it becomes suffix → ignored
    // If last run is uncovered, it's a suffix → ignored (pending_gap_length discarded)

    if !first_island_seen || n_islands == 0 {
        return None;
    }

    // internal_span: positions between first and last covered base, inclusive
    // Convention: last - first + 1 (bounds inclusive)
    let internal_span = last_covered_pos - first_covered_pos + 1;

    Some(SpatialScanResult {
        n_islands,
        total_internal_gap_bases,
        max_internal_gap,
        analysed_contig_length,
        internal_span,
        num_covered_bases,
    })
}

pub trait MosdepthGenomeCoverageEstimator {
    fn setup(&mut self);

    fn add_contig(
        &mut self,
        ups_and_downs: &[i32],
        num_mapped_reads: u64,
        total_mismatches: u64,
        sum_identity: f64,
    );

    fn calculate_coverage(&mut self, unobserved_contig_lengths: &[u64]) -> f32;

    fn print_coverage<T: CoverageTaker>(&self, coverage: f32, coverage_taker: &mut T);

    fn print_zero_coverage<T: CoverageTaker>(&self, coverage_taker: &mut T, entry_length: u64);

    fn copy(&self) -> CoverageEstimator;

    fn num_mapped_reads(&self) -> u64;
}

impl MosdepthGenomeCoverageEstimator for CoverageEstimator {
    fn setup(&mut self) {
        debug!("Running setup..");
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                ref mut total_count,
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ref mut total_mismatches,
                ..
            } => {
                *total_count = 0;
                *total_bases = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
                *total_mismatches = 0;
            }
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::VarianceGenomeCoverageEstimator {
                ref mut observed_contig_length,
                ref mut counts,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *counts = vec![];
                *observed_contig_length = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::RPKMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::TPMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *total_bases = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::ReferenceLengthCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            }
            | CoverageEstimator::ReadsPerBaseCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            } => {
                *observed_contig_length = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::ReadCountCalculator {
                ref mut num_mapped_reads,
            } => {
                *num_mapped_reads = 0;
            }
            CoverageEstimator::AverageIdentityEstimator {
                sum_identity,
                num_reads,
            } => {
                *sum_identity = 0.0;
                *num_reads = 0;
            }
            CoverageEstimator::StrobealignAembEstimator { .. } => panic!("Programming error"),
            CoverageEstimator::IslandsPerMbpEstimator {
                ref mut total_islands,
                ref mut covered_contigs_length,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *total_islands = 0;
                *covered_contigs_length = 0;
                *observed_contig_length = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::MaxGapEstimator {
                ref mut max_gap,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *max_gap = 0;
                *observed_contig_length = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
            CoverageEstimator::GapFractionEstimator {
                ref mut total_internal_gap_bases,
                ref mut total_internal_span,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *total_internal_gap_bases = 0;
                *total_internal_span = 0;
                *observed_contig_length = 0;
                *num_covered_bases = 0;
                *num_mapped_reads = 0;
            }
        }
    }

    fn add_contig(
        &mut self,
        ups_and_downs: &[i32],
        num_mapped_reads_in_contig: u64,
        total_mismatches_in_contig: u64,
        sum_identity_in_contig: f64,
    ) {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                ref mut total_count,
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ref mut total_mismatches,
                contig_end_exclusion,
                ..
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
                *total_mismatches += total_mismatches_in_contig;
                let len = ups_and_downs.len();
                match *contig_end_exclusion * 2 < len as u64 {
                    true => *total_bases += len as u64 - 2 * *contig_end_exclusion,
                    false => {
                        debug!("Contig too short - less than twice the contig-end-exclusion");
                        return; //contig is all ends, too short
                    }
                }
                let mut cumulative_sum: i32 = 0;
                let start_from = *contig_end_exclusion as usize;
                let end_at = len - *contig_end_exclusion as usize - 1;
                for (i, current) in ups_and_downs.iter().enumerate() {
                    cumulative_sum += current;
                    if i >= start_from && i <= end_at {
                        if cumulative_sum > 0 {
                            *num_covered_bases += 1
                        }
                        *total_count += cumulative_sum as u64;
                    }
                }
                debug!(
                    "After adding contig, have total_count {total_count}, total_bases {total_bases}, \
                        num_covered_bases {num_covered_bases}, mismatches {total_mismatches_in_contig}"
                );
            }
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                ..
            }
            | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                ..
            }
            | CoverageEstimator::VarianceGenomeCoverageEstimator {
                ref mut counts,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                ..
            } => {
                *num_mapped_reads = num_mapped_reads_in_contig;
                let len1 = ups_and_downs.len();
                match *contig_end_exclusion * 2 < len1 as u64 {
                    true => {
                        debug!("Adding len1 {len1}");
                        *observed_contig_length += len1 as u64 - 2 * *contig_end_exclusion
                    }
                    false => {
                        debug!("Contig too short - less than twice the contig-end-exclusion");
                        return; //contig is all ends, too short
                    }
                }
                debug!("Total observed length now {}", *observed_contig_length);
                let mut cumulative_sum: i32 = 0;
                let start_from = *contig_end_exclusion as usize;
                let end_at = len1 - *contig_end_exclusion as usize - 1;
                debug!("ups and downs {ups_and_downs:?}");
                for (i, current) in ups_and_downs.iter().enumerate() {
                    if *current != 0 {
                        debug!("At i {i}, cumulative sum {cumulative_sum} and current {current}");
                    }
                    cumulative_sum += current;
                    if i >= start_from && i <= end_at {
                        if cumulative_sum > 0 {
                            *num_covered_bases += 1
                        }
                        if counts.len() <= cumulative_sum as usize {
                            (*counts).resize(cumulative_sum as usize + 1, 0);
                        }
                        (*counts)[cumulative_sum as usize] += 1
                    }
                }
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::RPKMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            }
            | CoverageEstimator::TPMCoverageEstimator {
                ref mut total_bases,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                ..
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
                let len = ups_and_downs.len();
                *total_bases += len as u64;
                let mut cumulative_sum: i32 = 0;

                for current in ups_and_downs.iter() {
                    cumulative_sum += current;
                    if cumulative_sum > 0 {
                        *num_covered_bases += 1
                    }
                }
            }
            CoverageEstimator::ReferenceLengthCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            }
            | CoverageEstimator::ReadsPerBaseCalculator {
                ref mut observed_contig_length,
                ref mut num_mapped_reads,
            } => {
                *observed_contig_length += ups_and_downs.len() as u64;
                *num_mapped_reads += num_mapped_reads_in_contig;
            }
            CoverageEstimator::ReadCountCalculator {
                ref mut num_mapped_reads,
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
            }
            CoverageEstimator::AverageIdentityEstimator {
                sum_identity,
                num_reads,
            } => {
                *num_reads += num_mapped_reads_in_contig;
                *sum_identity += sum_identity_in_contig;
            }
            CoverageEstimator::StrobealignAembEstimator { .. } => unreachable!(),
            CoverageEstimator::IslandsPerMbpEstimator {
                ref mut total_islands,
                ref mut covered_contigs_length,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                min_island_length,
                ..
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
                let len = ups_and_downs.len();
                if *contig_end_exclusion * 2 >= len as u64 {
                    return; // contig too short
                }
                let contig_analysed = len as u64 - 2 * *contig_end_exclusion;
                *observed_contig_length += contig_analysed;
                match spatial_scan(ups_and_downs, *contig_end_exclusion, *min_island_length) {
                    Some(result) => {
                        *num_covered_bases += result.num_covered_bases;
                        *total_islands += result.n_islands;
                        *covered_contigs_length += result.analysed_contig_length;
                    }
                    None => {
                        // Count raw covered bases even if no valid island
                        let start = *contig_end_exclusion as usize;
                        let end = len - *contig_end_exclusion as usize - 1;
                        let mut cs: i32 = 0;
                        for (i, v) in ups_and_downs.iter().enumerate() {
                            cs += v;
                            if i >= start && i <= end && cs > 0 {
                                *num_covered_bases += 1;
                            }
                        }
                    }
                }
            }
            CoverageEstimator::MaxGapEstimator {
                ref mut max_gap,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                min_island_length,
                ..
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
                let len = ups_and_downs.len();
                if *contig_end_exclusion * 2 >= len as u64 {
                    return;
                }
                let contig_analysed = len as u64 - 2 * *contig_end_exclusion;
                *observed_contig_length += contig_analysed;
                match spatial_scan(ups_and_downs, *contig_end_exclusion, *min_island_length) {
                    Some(result) => {
                        *num_covered_bases += result.num_covered_bases;
                        if result.max_internal_gap > *max_gap {
                            *max_gap = result.max_internal_gap;
                        }
                    }
                    None => {
                        let start = *contig_end_exclusion as usize;
                        let end = len - *contig_end_exclusion as usize - 1;
                        let mut cs: i32 = 0;
                        for (i, v) in ups_and_downs.iter().enumerate() {
                            cs += v;
                            if i >= start && i <= end && cs > 0 {
                                *num_covered_bases += 1;
                            }
                        }
                    }
                }
            }
            CoverageEstimator::GapFractionEstimator {
                ref mut total_internal_gap_bases,
                ref mut total_internal_span,
                ref mut observed_contig_length,
                ref mut num_covered_bases,
                ref mut num_mapped_reads,
                contig_end_exclusion,
                min_island_length,
                ..
            } => {
                *num_mapped_reads += num_mapped_reads_in_contig;
                let len = ups_and_downs.len();
                if *contig_end_exclusion * 2 >= len as u64 {
                    return;
                }
                let contig_analysed = len as u64 - 2 * *contig_end_exclusion;
                *observed_contig_length += contig_analysed;
                match spatial_scan(ups_and_downs, *contig_end_exclusion, *min_island_length) {
                    Some(result) => {
                        *num_covered_bases += result.num_covered_bases;
                        *total_internal_gap_bases += result.total_internal_gap_bases;
                        *total_internal_span += result.internal_span;
                    }
                    None => {
                        let start = *contig_end_exclusion as usize;
                        let end = len - *contig_end_exclusion as usize - 1;
                        let mut cs: i32 = 0;
                        for (i, v) in ups_and_downs.iter().enumerate() {
                            cs += v;
                            if i >= start && i <= end && cs > 0 {
                                *num_covered_bases += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    fn calculate_coverage(&mut self, unobserved_contig_lengths: &[u64]) -> f32 {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count,
                total_bases,
                num_covered_bases,
                num_mapped_reads: _,
                total_mismatches,
                contig_end_exclusion,
                min_fraction_covered_bases,
                exclude_mismatches,
            } => {
                let final_total_bases = *total_bases
                    + CoverageEstimator::calculate_unobserved_bases(
                        unobserved_contig_lengths,
                        *contig_end_exclusion,
                    );
                debug!(
                    "Calculating coverage with unobserved {unobserved_contig_lengths:?}, \
                        total bases {total_bases}, num_covered_bases {num_covered_bases}, total_count {total_count}, \
                        total_mismatches {total_mismatches}, final_total_bases {final_total_bases}",
                );
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    0.0
                } else {
                    let calculated_coverage = match exclude_mismatches {
                        true => (*total_count - *total_mismatches) as f32,
                        false => *total_count as f32,
                    } / final_total_bases as f32;
                    debug!("Found mean coverage {calculated_coverage}");
                    calculated_coverage
                }
            }
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts,
                observed_contig_length,
                num_covered_bases,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
                min,
                max,
            } => {
                let unobserved_contig_length = CoverageEstimator::calculate_unobserved_bases(
                    unobserved_contig_lengths,
                    *contig_end_exclusion,
                );
                let total_bases = *observed_contig_length + unobserved_contig_length;
                debug!("Calculating coverage with num_covered_bases {num_covered_bases}, observed_length {observed_contig_length}, unobserved_length {unobserved_contig_lengths:?} and counts {counts:?}");

                match total_bases {
                    0 => 0.0,
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32)
                            < *min_fraction_covered_bases
                        {
                            0.0
                        } else {
                            let min_index: usize = (*min * total_bases as f32).floor() as usize;
                            let max_index: usize = (*max * total_bases as f32).ceil() as usize;
                            if *num_covered_bases == 0 {
                                return 0.0;
                            }
                            counts[0] += unobserved_contig_length;

                            let mut num_accounted_for: usize = 0;
                            let mut total: usize = 0;
                            let mut started = false;
                            for (i, num_covered) in counts.iter().enumerate() {
                                num_accounted_for += *num_covered as usize;
                                debug!(
                                    "start: i {i}, num_accounted_for {num_accounted_for}, total {total}, min {min_index}, max {max_index}"
                                );
                                if num_accounted_for >= min_index {
                                    debug!("inside");
                                    if started {
                                        if num_accounted_for > max_index {
                                            debug!(
                                                "num_accounted_for {}, *num_covered {}",
                                                num_accounted_for, *num_covered
                                            );
                                            let num_excess =
                                                num_accounted_for - *num_covered as usize;
                                            let num_wanted = match max_index >= num_excess {
                                                true => max_index - num_excess + 1,
                                                false => 0,
                                            };
                                            debug!("num wanted1: {num_wanted}");
                                            total += num_wanted * i;
                                            break;
                                        } else {
                                            total += *num_covered as usize * i;
                                        }
                                    } else if num_accounted_for > max_index {
                                        // all coverages are the same in the trimmed set
                                        total = (max_index - min_index + 1) * i;
                                        started = true
                                    } else if num_accounted_for < min_index {
                                        debug!("too few on first")
                                    } else {
                                        let num_wanted = num_accounted_for - min_index + 1;
                                        debug!("num wanted2: {num_wanted}");
                                        total = num_wanted * i;
                                        started = true;
                                    }
                                }
                                debug!(
                                    "end i {i}, num_accounted_for {num_accounted_for}, total {total}"
                                );
                            }
                            total as f32 / (max_index - min_index) as f32
                        }
                    }
                }
            }
            CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length,
                num_covered_bases,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => {
                // No need to actually calculate any kind of coverage, just return
                // whether any coverage was detected
                match observed_contig_length {
                    0 => 0.0,
                    _ => {
                        let total_bases = *observed_contig_length
                            + CoverageEstimator::calculate_unobserved_bases(
                                unobserved_contig_lengths,
                                *contig_end_exclusion,
                            );
                        if (*num_covered_bases as f32 / total_bases as f32)
                            < *min_fraction_covered_bases
                        {
                            0.0
                        } else {
                            // Hack: Return the number of zero coverage bases as the
                            // coverage, plus 1 so it is definitely non-zero, so that
                            // the print_genome function knows this info.
                            (total_bases - *num_covered_bases + 1) as f32
                        }
                    }
                }
            }
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    0.0
                } else {
                    *num_covered_bases as f32 / final_total_bases as f32
                }
            }
            CoverageEstimator::NumCoveredBasesCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    0.0
                } else {
                    *num_covered_bases as f32
                }
            }
            CoverageEstimator::RPKMCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    0.0
                } else {
                    // Here we do not know the number of mapped reads total.
                    // Instead we divide by that later.
                    debug!("RPKM: {num_mapped_reads} {final_total_bases}");
                    match final_total_bases == 0 {
                        true => 0.0,
                        false => {
                            (*num_mapped_reads * (10u64.pow(9))) as f32 / final_total_bases as f32
                        }
                    }
                }
            }
            CoverageEstimator::TPMCoverageEstimator {
                total_bases,
                num_covered_bases,
                num_mapped_reads,
                min_fraction_covered_bases,
            } => {
                let final_total_bases: u64 =
                    *total_bases + unobserved_contig_lengths.iter().sum::<u64>();
                if final_total_bases == 0
                    || (*num_covered_bases as f32 / final_total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    0.0
                } else {
                    // Here we do not know the number of mapped reads total.
                    // Instead we divide by that later.
                    debug!("TPM: {num_mapped_reads} {final_total_bases}");
                    match final_total_bases == 0 {
                        true => 0.0,
                        false => {
                            ((*num_mapped_reads as f64).ln() - (final_total_bases as f64).ln())
                                .exp() as f32
                        }
                    }
                }
            }
            CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length,
                ref mut counts,
                num_covered_bases,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => {
                let unobserved_contig_length = CoverageEstimator::calculate_unobserved_bases(
                    unobserved_contig_lengths,
                    *contig_end_exclusion,
                );
                let total_bases = *observed_contig_length + unobserved_contig_length;
                debug!("Calculating coverage with observed length {num_covered_bases}, unobserved_length {unobserved_contig_lengths:?} and counts {counts:?}");
                match total_bases {
                    0 => 0.0,
                    _ => {
                        if (*num_covered_bases as f32 / total_bases as f32)
                            < *min_fraction_covered_bases || total_bases < 3 ||
                            // no mapped reads
                            counts.is_empty()
                        {
                            0.0
                        } else {
                            counts[0] += unobserved_contig_length;
                            // Calculate variance using the shifted method
                            let mut k = 0;
                            // Ensure K is within the range of coverages - take the
                            // lowest coverage.
                            while counts[k] == 0 {
                                k += 1;
                            }
                            let mut ex = 0;
                            let mut ex2 = 0;
                            for (x, num_covered) in counts.iter().enumerate() {
                                if *num_covered == 0 {
                                    continue;
                                }
                                let nc = *num_covered as usize;
                                ex += (x - k) * nc;
                                ex2 += (x - k) * (x - k) * nc;
                            }
                            // Return sample variance not population variance since
                            // almost all MAGs are incomplete.
                            (ex2 as f32 - (ex * ex) as f32 / total_bases as f32)
                                / (total_bases - 1) as f32
                        }
                    }
                }
            }
            CoverageEstimator::ReferenceLengthCalculator {
                observed_contig_length,
                ..
            } => (*observed_contig_length + unobserved_contig_lengths.iter().sum::<u64>()) as f32,
            CoverageEstimator::ReadCountCalculator { num_mapped_reads } => *num_mapped_reads as f32,
            CoverageEstimator::ReadsPerBaseCalculator {
                observed_contig_length,
                num_mapped_reads,
            } => {
                *num_mapped_reads as f32
                    / (*observed_contig_length + unobserved_contig_lengths.iter().sum::<u64>())
                        as f32
            }
            CoverageEstimator::AverageIdentityEstimator {
                sum_identity,
                num_reads,
            } => {
                if *num_reads == 0 {
                    0.0
                } else {
                    (*sum_identity / *num_reads as f64) as f32
                }
            }
            CoverageEstimator::StrobealignAembEstimator {} => unreachable!(),
            CoverageEstimator::IslandsPerMbpEstimator {
                total_islands,
                covered_contigs_length,
                observed_contig_length,
                num_covered_bases,
                contig_end_exclusion,
                min_fraction_covered_bases,
                ..
            } => {
                let total_bases = *observed_contig_length
                    + CoverageEstimator::calculate_unobserved_bases(
                        unobserved_contig_lengths,
                        *contig_end_exclusion,
                    );
                if total_bases == 0
                    || (*num_covered_bases as f32 / total_bases as f32)
                        < *min_fraction_covered_bases
                    || *covered_contigs_length == 0
                {
                    0.0
                } else {
                    *total_islands as f32 / (*covered_contigs_length as f32 / 1_000_000.0)
                }
            }
            CoverageEstimator::MaxGapEstimator {
                max_gap,
                observed_contig_length,
                num_covered_bases,
                contig_end_exclusion,
                min_fraction_covered_bases,
                ..
            } => {
                let total_bases = *observed_contig_length
                    + CoverageEstimator::calculate_unobserved_bases(
                        unobserved_contig_lengths,
                        *contig_end_exclusion,
                    );
                if total_bases == 0
                    || (*num_covered_bases as f32 / total_bases as f32)
                        < *min_fraction_covered_bases
                {
                    0.0
                } else {
                    // Note: u64 → f32 conversion loses precision beyond 16,777,216.
                    // Acceptable for MAG contigs (rarely >10M bases).
                    *max_gap as f32
                }
            }
            CoverageEstimator::GapFractionEstimator {
                total_internal_gap_bases,
                total_internal_span,
                observed_contig_length,
                num_covered_bases,
                contig_end_exclusion,
                min_fraction_covered_bases,
                ..
            } => {
                let total_bases = *observed_contig_length
                    + CoverageEstimator::calculate_unobserved_bases(
                        unobserved_contig_lengths,
                        *contig_end_exclusion,
                    );
                if total_bases == 0
                    || (*num_covered_bases as f32 / total_bases as f32)
                        < *min_fraction_covered_bases
                    || *total_internal_span == 0
                {
                    0.0
                } else {
                    *total_internal_gap_bases as f32 / *total_internal_span as f32
                }
            }
        }
    }

    fn copy(&self) -> CoverageEstimator {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count: _,
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                total_mismatches: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
                exclude_mismatches,
            } => CoverageEstimator::new_estimator_mean(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
                *exclude_mismatches,
            ),
            CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
                min,
                max,
            } => CoverageEstimator::new_estimator_trimmed_mean(
                *min,
                *max,
                *min_fraction_covered_bases,
                *contig_end_exclusion,
            ),
            CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_pileup_counts(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
            ),
            CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_covered_fraction(*min_fraction_covered_bases),
            CoverageEstimator::NumCoveredBasesCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_covered_bases(*min_fraction_covered_bases),
            CoverageEstimator::RPKMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_rpkm(*min_fraction_covered_bases),
            CoverageEstimator::TPMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_tpm(*min_fraction_covered_bases),
            CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length: _,
                counts: _,
                num_covered_bases: _,
                num_mapped_reads: _,
                contig_end_exclusion,
                min_fraction_covered_bases,
            } => CoverageEstimator::new_estimator_variance(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
            ),
            CoverageEstimator::ReferenceLengthCalculator { .. } => {
                CoverageEstimator::new_estimator_length()
            }
            CoverageEstimator::ReadCountCalculator { .. } => {
                CoverageEstimator::new_estimator_read_count()
            }
            CoverageEstimator::ReadsPerBaseCalculator { .. } => {
                CoverageEstimator::new_estimator_reads_per_base()
            }
            CoverageEstimator::AverageIdentityEstimator { .. } => {
                CoverageEstimator::new_estimator_anir()
            }
            CoverageEstimator::StrobealignAembEstimator { .. } => {
                CoverageEstimator::new_estimator_strobealign_aemb()
            }
            CoverageEstimator::IslandsPerMbpEstimator {
                min_fraction_covered_bases,
                contig_end_exclusion,
                min_island_length,
                ..
            } => CoverageEstimator::new_estimator_islands_per_mbp(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
                *min_island_length,
            ),
            CoverageEstimator::MaxGapEstimator {
                min_fraction_covered_bases,
                contig_end_exclusion,
                min_island_length,
                ..
            } => CoverageEstimator::new_estimator_max_gap(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
                *min_island_length,
            ),
            CoverageEstimator::GapFractionEstimator {
                min_fraction_covered_bases,
                contig_end_exclusion,
                min_island_length,
                ..
            } => CoverageEstimator::new_estimator_gap_fraction(
                *min_fraction_covered_bases,
                *contig_end_exclusion,
                *min_island_length,
            ),
        }
    }

    fn print_coverage<T: CoverageTaker>(&self, coverage: f32, coverage_taker: &mut T) {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::TrimmedMeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::CoverageFractionGenomeCoverageEstimator { .. }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator { .. }
            | CoverageEstimator::RPKMCoverageEstimator { .. }
            | CoverageEstimator::TPMCoverageEstimator { .. }
            | CoverageEstimator::VarianceGenomeCoverageEstimator { .. }
            | CoverageEstimator::ReferenceLengthCalculator { .. }
            | CoverageEstimator::ReadCountCalculator { .. }
            | CoverageEstimator::ReadsPerBaseCalculator { .. }
            | CoverageEstimator::AverageIdentityEstimator { .. }
            | CoverageEstimator::StrobealignAembEstimator { .. }
            | CoverageEstimator::IslandsPerMbpEstimator { .. }
            | CoverageEstimator::MaxGapEstimator { .. }
            | CoverageEstimator::GapFractionEstimator { .. } => {
                coverage_taker.add_single_coverage(coverage);
            }
            CoverageEstimator::PileupCountsGenomeCoverageEstimator { counts, .. } => {
                debug!("{counts:?}");
                for (i, num_covered) in counts.iter().enumerate() {
                    let cov: u64 = match i {
                        0 => {
                            let c = coverage.floor() as u64;
                            match c {
                                0 => 0,
                                _ => c - 1,
                            }
                        }
                        _ => *num_covered,
                    };
                    coverage_taker.add_coverage_entry(i, cov);
                }
            }
        }
    }

    fn print_zero_coverage<T: CoverageTaker>(&self, coverage_taker: &mut T, entry_length: u64) {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::TrimmedMeanGenomeCoverageEstimator { .. }
            | CoverageEstimator::CoverageFractionGenomeCoverageEstimator { .. }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator { .. }
            | CoverageEstimator::RPKMCoverageEstimator { .. }
            | CoverageEstimator::TPMCoverageEstimator { .. }
            | CoverageEstimator::VarianceGenomeCoverageEstimator { .. }
            | CoverageEstimator::ReadCountCalculator { .. }
            | CoverageEstimator::ReadsPerBaseCalculator { .. }
            | CoverageEstimator::AverageIdentityEstimator { .. }
            | CoverageEstimator::StrobealignAembEstimator { .. }
            | CoverageEstimator::IslandsPerMbpEstimator { .. }
            | CoverageEstimator::MaxGapEstimator { .. }
            | CoverageEstimator::GapFractionEstimator { .. } => {
                coverage_taker.add_single_coverage(0.0);
            }
            CoverageEstimator::PileupCountsGenomeCoverageEstimator { .. } => {}
            CoverageEstimator::ReferenceLengthCalculator { .. } => {
                coverage_taker.add_single_coverage(entry_length as f32);
            }
        }
    }

    fn num_mapped_reads(&self) -> u64 {
        match self {
            CoverageEstimator::MeanGenomeCoverageEstimator {
                total_count: _,
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::TrimmedMeanGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::PileupCountsGenomeCoverageEstimator {
                counts: _,
                observed_contig_length: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::CoverageFractionGenomeCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::NumCoveredBasesCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::RPKMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::TPMCoverageEstimator {
                total_bases: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::VarianceGenomeCoverageEstimator {
                observed_contig_length: _,
                counts: _,
                num_covered_bases: _,
                num_mapped_reads,
                ..
            }
            | CoverageEstimator::ReferenceLengthCalculator {
                observed_contig_length: _,
                num_mapped_reads,
            }
            | CoverageEstimator::ReadCountCalculator { num_mapped_reads }
            | CoverageEstimator::ReadsPerBaseCalculator {
                observed_contig_length: _,
                num_mapped_reads,
            } => *num_mapped_reads,
            CoverageEstimator::AverageIdentityEstimator { num_reads, .. } => *num_reads,
            CoverageEstimator::StrobealignAembEstimator { .. } => {
                panic!("Strobealign AEMB does not calculate number of mapped reads")
            }
            CoverageEstimator::IslandsPerMbpEstimator {
                num_mapped_reads, ..
            }
            | CoverageEstimator::MaxGapEstimator {
                num_mapped_reads, ..
            }
            | CoverageEstimator::GapFractionEstimator {
                num_mapped_reads, ..
            } => *num_mapped_reads,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build an ups_and_downs vector from (start, end, depth_change) triples.
    /// Each triple adds +depth_change at start and -depth_change at end.
    fn make_ups_and_downs(len: usize, regions: &[(usize, usize, i32)]) -> Vec<i32> {
        let mut v = vec![0i32; len];
        for &(start, end, dc) in regions {
            v[start] += dc;
            if end < len {
                v[end] -= dc;
            }
        }
        v
    }

    #[test]
    fn test_uniform_coverage() {
        // 100bp contig, covered at 10x everywhere, no exclusion
        let ud = make_ups_and_downs(100, &[(0, 100, 10)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.total_internal_gap_bases, 0);
        assert_eq!(result.max_internal_gap, 0);
        assert_eq!(result.analysed_contig_length, 100);
        assert_eq!(result.internal_span, 100);
        assert_eq!(result.num_covered_bases, 100);
    }

    #[test]
    fn test_three_islands() {
        // 100bp contig, 3 islands: [10..30), [50..70), [80..95)
        // Gaps: [30..50) = 20bp, [70..80) = 10bp
        let ud = make_ups_and_downs(100, &[(10, 30, 5), (50, 70, 5), (80, 95, 5)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 3);
        assert_eq!(result.total_internal_gap_bases, 30); // 20 + 10
        assert_eq!(result.max_internal_gap, 20);
        // internal_span: from pos 10 to pos 94 = 85
        assert_eq!(result.internal_span, 85);
        assert_eq!(result.num_covered_bases, 55); // 20 + 20 + 15
    }

    #[test]
    fn test_no_coverage() {
        let ud = vec![0i32; 100];
        assert!(spatial_scan(&ud, 0, 1).is_none());
    }

    #[test]
    fn test_contig_too_short() {
        let ud = vec![0i32; 10];
        // exclusion = 6, 2*6 = 12 > 10
        assert!(spatial_scan(&ud, 6, 1).is_none());
    }

    #[test]
    fn test_prefix_suffix_not_gaps() {
        // 100bp contig, coverage only in the middle: [30..70)
        // Prefix [0..30) and suffix [70..100) should not be gaps
        let ud = make_ups_and_downs(100, &[(30, 70, 5)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.total_internal_gap_bases, 0);
        assert_eq!(result.max_internal_gap, 0);
        assert_eq!(result.internal_span, 40); // 30..69 inclusive = 40
    }

    #[test]
    fn test_min_island_length_filters_micro_islands() {
        // 200bp contig, 3 segments: [10..20)=10bp, [50..60)=10bp, [100..180)=80bp
        // With min_island_length=50, first two are micro-islands
        let ud = make_ups_and_downs(200, &[(10, 20, 5), (50, 60, 5), (100, 180, 5)]);
        let result = spatial_scan(&ud, 0, 50).unwrap();
        assert_eq!(result.n_islands, 1); // only the 80bp island survives
        assert_eq!(result.total_internal_gap_bases, 0); // no gap between valid islands
        assert_eq!(result.max_internal_gap, 0);
        assert_eq!(result.internal_span, 80); // just the one island
    }

    #[test]
    fn test_all_micro_islands_returns_none() {
        // 100bp contig, all segments < min_island_length
        let ud = make_ups_and_downs(100, &[(10, 15, 5), (50, 55, 5), (80, 85, 5)]);
        assert!(spatial_scan(&ud, 0, 50).is_none());
    }

    #[test]
    fn test_single_island_with_micro_islands_around() {
        // 200bp, micro-island at [10..15), valid island at [50..150), micro at [170..175)
        let ud = make_ups_and_downs(200, &[(10, 15, 5), (50, 150, 5), (170, 175, 5)]);
        let result = spatial_scan(&ud, 0, 10).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.total_internal_gap_bases, 0);
        assert_eq!(result.gap_fraction(), 0.0);
    }

    #[test]
    fn test_contig_end_exclusion() {
        // 100bp contig, coverage everywhere, exclusion=10
        // Analysed region: [10..89] = 80 positions
        let ud = make_ups_and_downs(100, &[(0, 100, 10)]);
        let result = spatial_scan(&ud, 10, 1).unwrap();
        assert_eq!(result.analysed_contig_length, 80);
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.internal_span, 80);
    }

    #[test]
    fn test_two_islands_with_exclusion() {
        // 100bp contig, exclusion=5. Coverage at [10..30) and [60..80).
        // Analysed region: [5..94]. Both islands are within.
        // Gap: [30..60) = 30bp
        let ud = make_ups_and_downs(100, &[(10, 30, 5), (60, 80, 5)]);
        let result = spatial_scan(&ud, 5, 1).unwrap();
        assert_eq!(result.n_islands, 2);
        assert_eq!(result.total_internal_gap_bases, 30);
        assert_eq!(result.max_internal_gap, 30);
        // span: 10..79 = 70
        assert_eq!(result.internal_span, 70);
    }

    // Helper method for test convenience
    impl SpatialScanResult {
        fn gap_fraction(&self) -> f32 {
            if self.internal_span == 0 {
                0.0
            } else {
                self.total_internal_gap_bases as f32 / self.internal_span as f32
            }
        }
    }

    // =========================================================================
    // Edge cases for spatial_scan()
    // =========================================================================

    #[test]
    fn test_single_base_island() {
        // 50bp contig, single base covered at position 25
        let ud = make_ups_and_downs(50, &[(25, 26, 1)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.internal_span, 1);
        assert_eq!(result.total_internal_gap_bases, 0);
        assert_eq!(result.num_covered_bases, 1);
    }

    #[test]
    fn test_single_base_gap() {
        // 50bp contig, two islands separated by exactly 1 base: [10..20), [21..30)
        let ud = make_ups_and_downs(50, &[(10, 20, 5), (21, 30, 5)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 2);
        assert_eq!(result.total_internal_gap_bases, 1);
        assert_eq!(result.max_internal_gap, 1);
        // span: 10..29 = 20
        assert_eq!(result.internal_span, 20);
    }

    #[test]
    fn test_adjacent_islands_no_gap() {
        // Two islands that are exactly adjacent: [10..20), [20..30)
        let ud = make_ups_and_downs(50, &[(10, 20, 5), (20, 30, 3)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        // They should merge into one island (no gap between them)
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.total_internal_gap_bases, 0);
        assert_eq!(result.max_internal_gap, 0);
        assert_eq!(result.internal_span, 20);
        assert_eq!(result.num_covered_bases, 20);
    }

    #[test]
    fn test_overlapping_reads_varying_depth() {
        // 100bp contig, overlapping reads creating varying depth
        // [10..50) at 3x, [30..70) at additional 2x → total [10..30)=3x, [30..50)=5x, [50..70)=2x
        let ud = make_ups_and_downs(100, &[(10, 50, 3), (30, 70, 2)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        // All covered from 10..70, one island
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.total_internal_gap_bases, 0);
        assert_eq!(result.internal_span, 60);
        assert_eq!(result.num_covered_bases, 60);
    }

    #[test]
    fn test_exclusion_cuts_into_island() {
        // 100bp contig, island [0..100) full coverage, exclusion=30
        // Analysed region: [30..69] = 40bp, still one island
        let ud = make_ups_and_downs(100, &[(0, 100, 5)]);
        let result = spatial_scan(&ud, 30, 1).unwrap();
        assert_eq!(result.analysed_contig_length, 40);
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.internal_span, 40);
    }

    #[test]
    fn test_exclusion_removes_island() {
        // 100bp contig, island only at [5..10), exclusion=10
        // Analysed region: [10..89] — island is outside → None
        let ud = make_ups_and_downs(100, &[(5, 10, 5)]);
        assert!(spatial_scan(&ud, 10, 1).is_none());
    }

    #[test]
    fn test_exclusion_exactly_at_boundary() {
        // 100bp contig, exclusion=10, island starts exactly at exclusion boundary [10..50)
        let ud = make_ups_and_downs(100, &[(10, 50, 5)]);
        let result = spatial_scan(&ud, 10, 1).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.num_covered_bases, 40);
    }

    #[test]
    fn test_micro_island_between_two_valid_islands() {
        // Two valid islands with a micro-island in the gap between them
        // [10..60), micro [70..75), [100..150)
        // min_island_length=20 → micro at [70..75) is reclassified
        // Gap should be 10..60 → gap from 60 to 100 = 40bp (including the 5bp micro)
        let ud = make_ups_and_downs(200, &[(10, 60, 5), (70, 75, 5), (100, 150, 5)]);
        let result = spatial_scan(&ud, 0, 20).unwrap();
        assert_eq!(result.n_islands, 2);
        // Internal gap: [60..100) but micro [70..75) is absorbed → total = 40
        assert_eq!(result.total_internal_gap_bases, 40);
        assert_eq!(result.max_internal_gap, 40);
    }

    #[test]
    fn test_two_micro_islands_between_valid_islands() {
        // [10..60), micro [70..75), micro [80..85), [100..150)
        // min_island_length=20 → both micros reclassified
        let ud = make_ups_and_downs(200, &[(10, 60, 5), (70, 75, 5), (80, 85, 5), (100, 150, 5)]);
        let result = spatial_scan(&ud, 0, 20).unwrap();
        assert_eq!(result.n_islands, 2);
        // Gap: 60 to 100 = 40bp (both micros absorbed into gap)
        assert_eq!(result.total_internal_gap_bases, 40);
    }

    #[test]
    fn test_island_at_very_end_of_contig() {
        // 100bp contig, island at the very end: [90..100)
        let ud = make_ups_and_downs(100, &[(90, 100, 5)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.internal_span, 10);
        assert_eq!(result.total_internal_gap_bases, 0);
    }

    #[test]
    fn test_island_at_very_start_of_contig() {
        // 100bp contig, island at the very start: [0..10)
        let ud = make_ups_and_downs(100, &[(0, 10, 5)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.internal_span, 10);
    }

    #[test]
    fn test_many_small_islands_high_fragmentation() {
        // Simulate false positive: 10 small islands of 5bp each across 1000bp
        let mut regions: Vec<(usize, usize, i32)> = vec![];
        for i in 0..10 {
            let start = i * 100;
            regions.push((start, start + 5, 3));
        }
        let ud = make_ups_and_downs(1000, &regions);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 10);
        assert_eq!(result.num_covered_bases, 50);
        // Internal gaps: 9 gaps of 95bp each = 855
        assert_eq!(result.total_internal_gap_bases, 855);
        assert_eq!(result.max_internal_gap, 95);
        // Span: 0..904 = 905
        assert_eq!(result.internal_span, 905);
        // gap_fraction should be very high
        let gf = result.gap_fraction();
        assert!(gf > 0.9, "gap_fraction should be >0.9, got {}", gf);
    }

    #[test]
    fn test_continuous_coverage_low_fragmentation() {
        // Simulate true MAG: one large island covering most of the contig
        let ud = make_ups_and_downs(1000, &[(10, 980, 8)]);
        let result = spatial_scan(&ud, 0, 1).unwrap();
        assert_eq!(result.n_islands, 1);
        assert_eq!(result.total_internal_gap_bases, 0);
        assert_eq!(result.gap_fraction(), 0.0);
    }

    // =========================================================================
    // Estimator integration tests (add_contig → calculate_coverage)
    // =========================================================================

    #[test]
    fn test_islands_per_mbp_estimator_single_contig() {
        let mut est = CoverageEstimator::new_estimator_islands_per_mbp(0.0, 0, 1);
        // 1000bp contig, 3 islands
        let ud = make_ups_and_downs(1000, &[(10, 100, 5), (200, 400, 5), (600, 900, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        // 3 islands / (1000bp / 1e6) = 3 / 0.001 = 3000
        assert!((coverage - 3000.0).abs() < 1.0, "got {}", coverage);
    }

    #[test]
    fn test_islands_per_mbp_estimator_multi_contig() {
        let mut est = CoverageEstimator::new_estimator_islands_per_mbp(0.0, 0, 1);

        // Contig 1: 1000bp, 2 islands
        let ud1 = make_ups_and_downs(1000, &[(10, 200, 5), (500, 800, 5)]);
        est.add_contig(&ud1, 10, 0, 0.0);

        // Contig 2: 2000bp, 1 island
        let ud2 = make_ups_and_downs(2000, &[(100, 1500, 5)]);
        est.add_contig(&ud2, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        // 3 islands / (3000bp / 1e6) = 3 / 0.003 = 1000
        assert!((coverage - 1000.0).abs() < 1.0, "got {}", coverage);
    }

    #[test]
    fn test_islands_per_mbp_estimator_uncovered_contig_excluded() {
        let mut est = CoverageEstimator::new_estimator_islands_per_mbp(0.0, 0, 1);

        // Contig 1: 1000bp, 2 islands
        let ud1 = make_ups_and_downs(1000, &[(10, 200, 5), (500, 800, 5)]);
        est.add_contig(&ud1, 10, 0, 0.0);

        // Contig 2: 2000bp, no coverage
        let ud2 = vec![0i32; 2000];
        est.add_contig(&ud2, 0, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        // 2 islands / (1000bp / 1e6) = 2000 — uncovered contig not in denominator
        assert!((coverage - 2000.0).abs() < 1.0, "got {}", coverage);
    }

    #[test]
    fn test_max_gap_estimator_multi_contig() {
        let mut est = CoverageEstimator::new_estimator_max_gap(0.0, 0, 1);

        // Contig 1: gap of 200bp
        let ud1 = make_ups_and_downs(500, &[(10, 50, 5), (250, 400, 5)]);
        est.add_contig(&ud1, 10, 0, 0.0);

        // Contig 2: gap of 500bp (larger)
        let ud2 = make_ups_and_downs(1000, &[(10, 50, 5), (550, 800, 5)]);
        est.add_contig(&ud2, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        // max_gap should be 500 (from contig 2)
        assert!((coverage - 500.0).abs() < 1.0, "got {}", coverage);
    }

    #[test]
    fn test_max_gap_estimator_no_gaps() {
        let mut est = CoverageEstimator::new_estimator_max_gap(0.0, 0, 1);
        let ud = make_ups_and_downs(100, &[(0, 100, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);
        let coverage = est.calculate_coverage(&[]);
        assert_eq!(coverage, 0.0);
    }

    #[test]
    fn test_gap_fraction_estimator() {
        let mut est = CoverageEstimator::new_estimator_gap_fraction(0.0, 0, 1);

        // 1000bp contig: islands [10..100), [200..400), [600..900)
        // Internal gaps: [100..200)=100bp, [400..600)=200bp → total=300bp
        // Internal span: 10..899 = 890bp
        let ud = make_ups_and_downs(1000, &[(10, 100, 5), (200, 400, 5), (600, 900, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        let expected = 300.0 / 890.0;
        assert!(
            (coverage - expected).abs() < 0.001,
            "got {}, expected {}",
            coverage,
            expected
        );
    }

    #[test]
    fn test_gap_fraction_estimator_multi_contig() {
        let mut est = CoverageEstimator::new_estimator_gap_fraction(0.0, 0, 1);

        // Contig 1: span 90, gap 50
        let ud1 = make_ups_and_downs(200, &[(10, 50, 5), (100, 150, 5)]);
        est.add_contig(&ud1, 10, 0, 0.0);

        // Contig 2: span 80, gap 30
        let ud2 = make_ups_and_downs(200, &[(10, 60, 5), (90, 140, 5)]);
        est.add_contig(&ud2, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        // total_internal_gap = 50 + 30 = 80
        // total_internal_span = (149-10+1) + (139-10+1) = 140 + 130 = 270
        let expected = 80.0 / 270.0;
        assert!(
            (coverage - expected).abs() < 0.01,
            "got {}, expected {}",
            coverage,
            expected
        );
    }

    #[test]
    fn test_estimator_min_fraction_gate() {
        // Set min_fraction_covered_bases = 0.5
        let mut est = CoverageEstimator::new_estimator_islands_per_mbp(0.5, 0, 1);

        // 1000bp contig with only 100bp covered → 10% < 50% → should be gated
        let ud = make_ups_and_downs(1000, &[(100, 200, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        assert_eq!(coverage, 0.0, "should be 0.0 due to min_fraction gate");
    }

    #[test]
    fn test_estimator_min_fraction_gate_passes() {
        let mut est = CoverageEstimator::new_estimator_max_gap(0.5, 0, 1);

        // 1000bp contig with 800bp covered → 80% > 50% → gate passes
        let ud = make_ups_and_downs(1000, &[(10, 500, 5), (600, 810, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        // Gap of 100bp between islands
        assert!((coverage - 100.0).abs() < 1.0, "got {}", coverage);
    }

    #[test]
    fn test_estimator_with_unobserved_contigs() {
        let mut est = CoverageEstimator::new_estimator_islands_per_mbp(0.5, 0, 1);

        // 1000bp contig, 800bp covered
        let ud = make_ups_and_downs(1000, &[(100, 900, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        // Unobserved contig of 9000bp → total = 10000bp, covered = 800/10000 = 8% < 50%
        let coverage = est.calculate_coverage(&[9000]);
        assert_eq!(
            coverage, 0.0,
            "unobserved contigs should push below min_fraction gate"
        );
    }

    #[test]
    fn test_estimator_setup_resets() {
        let mut est = CoverageEstimator::new_estimator_gap_fraction(0.0, 0, 1);

        // Add some data
        let ud = make_ups_and_downs(100, &[(10, 50, 5), (60, 90, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        // Reset
        est.setup();

        // After reset, should return 0.0
        let coverage = est.calculate_coverage(&[]);
        assert_eq!(coverage, 0.0, "setup should reset all accumulators");
    }

    #[test]
    fn test_estimator_with_contig_end_exclusion() {
        let mut est = CoverageEstimator::new_estimator_max_gap(0.0, 10, 1);

        // 100bp contig, exclusion=10, islands at [5..15) and [85..95)
        // After exclusion [10..89]: first island is cut to [10..15)=5bp,
        // second island is cut to [85..89]=5bp
        // Gap: [15..85) = 70bp
        let ud = make_ups_and_downs(100, &[(5, 15, 5), (85, 95, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        assert!((coverage - 70.0).abs() < 1.0, "got {}", coverage);
    }

    #[test]
    fn test_min_island_length_in_estimator() {
        let mut est = CoverageEstimator::new_estimator_islands_per_mbp(0.0, 0, 50);

        // 1000bp contig: micro [10..15)=5bp, valid [100..800)=700bp, micro [900..905)=5bp
        let ud = make_ups_and_downs(1000, &[(10, 15, 5), (100, 800, 5), (900, 905, 5)]);
        est.add_contig(&ud, 10, 0, 0.0);

        let coverage = est.calculate_coverage(&[]);
        // Only 1 valid island, denominator = 1000bp
        // 1 / (1000/1e6) = 1000
        assert!((coverage - 1000.0).abs() < 1.0, "got {}", coverage);
    }

    // =========================================================================
    // Tests for column_headers(), copy(), num_mapped_reads()
    // =========================================================================

    #[test]
    fn test_column_headers() {
        let est1 = CoverageEstimator::new_estimator_islands_per_mbp(0.0, 0, 1);
        assert_eq!(est1.column_headers(), vec!["Islands per Mbp"]);

        let est2 = CoverageEstimator::new_estimator_max_gap(0.0, 0, 1);
        assert_eq!(est2.column_headers(), vec!["Max Gap"]);

        let est3 = CoverageEstimator::new_estimator_gap_fraction(0.0, 0, 1);
        assert_eq!(est3.column_headers(), vec!["Gap Fraction"]);
    }

    #[test]
    fn test_copy_preserves_config() {
        let est = CoverageEstimator::new_estimator_islands_per_mbp(0.3, 75, 50);
        // Add some data to ensure accumulators are non-zero
        let mut est_mut = est;
        let ud = make_ups_and_downs(1000, &[(10, 500, 5)]);
        est_mut.add_contig(&ud, 10, 0, 0.0);

        // Copy should reset accumulators but preserve config
        let mut copied = est_mut.copy();
        let coverage = copied.calculate_coverage(&[]);
        // Copied estimator has no data → should return 0.0
        assert_eq!(coverage, 0.0, "Copied estimator should have no data");

        // Verify config is preserved by checking it works with the same params
        let mut copied_mut = copied;
        copied_mut.add_contig(&ud, 10, 0, 0.0);
        let coverage = copied_mut.calculate_coverage(&[]);
        assert!(
            coverage > 0.0,
            "Copied estimator should work after adding data"
        );
    }

    #[test]
    fn test_copy_preserves_config_max_gap() {
        let est = CoverageEstimator::new_estimator_max_gap(0.0, 10, 25);
        let copied = est.copy();
        // Verify the copy works (config preserved)
        let mut copied_mut = copied;
        let ud = make_ups_and_downs(200, &[(20, 50, 5), (100, 150, 5)]);
        copied_mut.add_contig(&ud, 10, 0, 0.0);
        let coverage = copied_mut.calculate_coverage(&[]);
        assert!(coverage > 0.0, "Copied max_gap estimator should work");
    }

    #[test]
    fn test_copy_preserves_config_gap_fraction() {
        let est = CoverageEstimator::new_estimator_gap_fraction(0.0, 5, 10);
        let copied = est.copy();
        let mut copied_mut = copied;
        let ud = make_ups_and_downs(200, &[(20, 50, 5), (100, 150, 5)]);
        copied_mut.add_contig(&ud, 10, 0, 0.0);
        let coverage = copied_mut.calculate_coverage(&[]);
        assert!(coverage > 0.0, "Copied gap_fraction estimator should work");
    }

    #[test]
    fn test_num_mapped_reads_accumulates() {
        let mut est = CoverageEstimator::new_estimator_islands_per_mbp(0.0, 0, 1);

        let ud1 = make_ups_and_downs(1000, &[(10, 500, 5)]);
        est.add_contig(&ud1, 100, 0, 0.0);

        let ud2 = make_ups_and_downs(1000, &[(10, 500, 5)]);
        est.add_contig(&ud2, 200, 0, 0.0);

        // Should accumulate: 100 + 200 = 300
        assert_eq!(est.num_mapped_reads(), 300);
    }

    #[test]
    fn test_num_mapped_reads_max_gap() {
        let mut est = CoverageEstimator::new_estimator_max_gap(0.0, 0, 1);
        let ud = make_ups_and_downs(100, &[(10, 90, 5)]);
        est.add_contig(&ud, 42, 0, 0.0);
        assert_eq!(est.num_mapped_reads(), 42);
    }

    #[test]
    fn test_num_mapped_reads_gap_fraction() {
        let mut est = CoverageEstimator::new_estimator_gap_fraction(0.0, 0, 1);
        let ud = make_ups_and_downs(100, &[(10, 90, 5)]);
        est.add_contig(&ud, 55, 0, 0.0);
        est.add_contig(&ud, 45, 0, 0.0);
        assert_eq!(est.num_mapped_reads(), 100);
    }
}
