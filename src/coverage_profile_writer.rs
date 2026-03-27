use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Writes per-base coverage profiles as BigWig files.
///
/// The writer produces a raw, unfiltered coverage profile — no
/// contig_end_exclusion, no min_island_length. It's intentionally decoupled
/// from the spatial metrics estimators: the BigWig is primary data, the
/// metrics are derived statistics.
///
/// Implementation: coverage data is first written as BedGraph to a temporary
/// file, then converted to BigWig using the bigtools crate. This two-step
/// approach is necessary because BigWig requires all data to be written in
/// a specific format with indexes. The temporary BedGraph is very compact
/// (RLE encoded) so the overhead is minimal.
///
/// Note: this writer performs its own cumulative sum over ups_and_downs,
/// independently from the estimators. This is intentional — the bottleneck
/// is BAM I/O, not an in-memory vector scan.
pub struct CoverageProfileWriter {
    /// Temporary BedGraph file (plain text, uncompressed)
    temp_writer: std::io::BufWriter<std::fs::File>,
    temp_path: PathBuf,
    /// Final BigWig output path
    output_path: PathBuf,
    /// Chromosome sizes collected during writing (name → length)
    chrom_sizes: HashMap<String, u32>,
}

impl CoverageProfileWriter {
    /// Create a new BigWig writer at the given path.
    ///
    /// The final output will be a .bw file. A temporary BedGraph file is
    /// used during writing and removed after conversion.
    ///
    /// Panics if the output file already exists (no silent overwrite).
    pub fn new(path: &Path) -> Self {
        if path.exists() {
            panic!(
                "Coverage profile output file already exists: {}. \
                 Refusing to overwrite.",
                path.display()
            );
        }

        // Create a temporary BedGraph file alongside the output
        let temp_path = path.with_extension("bedgraph.tmp");
        let temp_file = std::fs::File::create(&temp_path).unwrap_or_else(|e| {
            panic!(
                "Failed to create temporary BedGraph file at {}: {}",
                temp_path.display(),
                e
            )
        });

        CoverageProfileWriter {
            temp_writer: std::io::BufWriter::new(temp_file),
            temp_path,
            output_path: path.to_path_buf(),
            chrom_sizes: HashMap::new(),
        }
    }

    /// Write the per-base coverage profile of a single contig in BedGraph RLE format.
    ///
    /// The profile is raw (no contig_end_exclusion applied) and run-length encoded:
    /// consecutive positions with the same depth are merged into one segment.
    pub fn write_contig(&mut self, contig_name: &str, ups_and_downs: &[i32]) {
        let len = ups_and_downs.len();
        if len == 0 {
            return;
        }

        // Record chromosome size for BigWig header
        self.chrom_sizes.insert(contig_name.to_string(), len as u32);

        let mut cumulative_sum: i32 = 0;
        let mut segment_start: usize = 0;
        let mut segment_depth: i32 = 0;

        for (i, &change) in ups_and_downs.iter().enumerate() {
            cumulative_sum += change;

            if i == 0 {
                segment_depth = cumulative_sum;
                continue;
            }

            if cumulative_sum != segment_depth {
                writeln!(
                    self.temp_writer,
                    "{}\t{}\t{}\t{}",
                    contig_name, segment_start, i, segment_depth
                )
                .unwrap_or_else(|e| panic!("Failed to write BedGraph: {}", e));
                segment_start = i;
                segment_depth = cumulative_sum;
            }
        }

        // Write the last segment
        writeln!(
            self.temp_writer,
            "{}\t{}\t{}\t{}",
            contig_name, segment_start, len, segment_depth
        )
        .unwrap_or_else(|e| panic!("Failed to write BedGraph: {}", e));
    }

    /// Convert the temporary BedGraph to BigWig and clean up.
    ///
    /// This flushes the BedGraph, converts it to BigWig using the bigtools crate,
    /// then removes the temporary file.
    pub fn finish(self) {
        // Flush and close the temp writer
        drop(self.temp_writer);

        if self.chrom_sizes.is_empty() {
            // No data was written — clean up temp file, don't create BigWig
            let _ = std::fs::remove_file(&self.temp_path);
            return;
        }

        // Convert BedGraph → BigWig using bigtools
        let infile = std::fs::File::open(&self.temp_path)
            .unwrap_or_else(|e| panic!("Failed to open temporary BedGraph for conversion: {}", e));

        let output_path = self.output_path;
        let temp_path = self.temp_path;

        let mut outfile = bigtools::BigWigWrite::create_file(&output_path, self.chrom_sizes)
            .unwrap_or_else(|e| {
                panic!(
                    "Failed to create BigWig file at {}: {}",
                    output_path.display(),
                    e
                )
            });
        // Contigs come in BAM tid order, not lexicographic order.
        // InputSortType::START allows out-of-order chromosomes while still
        // requiring sorted positions within each chromosome (which we guarantee).
        outfile.options.input_sort_type = bigtools::InputSortType::START;

        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .build()
            .expect("Failed to create tokio runtime for BigWig writing");

        // allow_out_of_order_chroms=true matches InputSortType::START
        let data = bigtools::beddata::BedParserStreamingIterator::from_bedgraph_file(infile, true);

        outfile.write(data, runtime).unwrap_or_else(|e| {
            panic!(
                "Failed to write BigWig file {}: {:?}",
                output_path.display(),
                e
            )
        });

        info!("Created BigWig coverage profile: {}", output_path.display());

        // Clean up temporary BedGraph
        let _ = std::fs::remove_file(&temp_path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    /// Helper: build ups_and_downs from (start, end, depth_change) triples
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

    /// Helper: read a BigWig file and return all intervals as (chrom, start, end, value)
    fn read_bigwig_intervals(path: &Path) -> Vec<(String, u32, u32, f32)> {
        let mut reader = bigtools::BigWigRead::open_file(path).unwrap();
        let chroms: Vec<_> = reader.chroms().to_vec();
        let mut intervals = vec![];
        for chrom in &chroms {
            let mut iter = reader.get_interval(&chrom.name, 0, chrom.length).unwrap();
            while let Some(Ok(val)) = iter.next() {
                intervals.push((chrom.name.clone(), val.start, val.end, val.value));
            }
        }
        intervals
    }

    #[test]
    fn test_uniform_coverage_bigwig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let ud = make_ups_and_downs(100, &[(0, 100, 10)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("contig_1", &ud);
        writer.finish();

        assert!(path.exists(), "BigWig file should exist");
        let intervals = read_bigwig_intervals(&path);
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0], ("contig_1".to_string(), 0, 100, 10.0));
    }

    #[test]
    fn test_rle_compression_bigwig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        // 100bp contig: [0..20) at 0x, [20..60) at 5x, [60..80) at 12x, [80..100) at 0x
        let ud = make_ups_and_downs(100, &[(20, 60, 5), (60, 80, 12)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        let intervals = read_bigwig_intervals(&path);
        // BigWig may omit 0-depth regions or include them, check non-zero at least
        let nonzero: Vec<_> = intervals.iter().filter(|i| i.3 > 0.0).collect();
        assert_eq!(nonzero.len(), 2);
        assert_eq!(nonzero[0], &("ctg".to_string(), 20, 60, 5.0));
        assert_eq!(nonzero[1], &("ctg".to_string(), 60, 80, 12.0));
    }

    #[test]
    fn test_multiple_contigs_bigwig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let ud1 = make_ups_and_downs(50, &[(10, 40, 8)]);
        let ud2 = make_ups_and_downs(30, &[(0, 30, 3)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("contig_A", &ud1);
        writer.write_contig("contig_B", &ud2);
        writer.finish();

        let reader = bigtools::BigWigRead::open_file(&path).unwrap();
        let chroms = reader.chroms().to_vec();
        assert_eq!(chroms.len(), 2);

        let chrom_names: Vec<_> = chroms.iter().map(|c| c.name.clone()).collect();
        assert!(chrom_names.contains(&"contig_A".to_string()));
        assert!(chrom_names.contains(&"contig_B".to_string()));
    }

    #[test]
    fn test_empty_contig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("empty", &[]);
        writer.finish();

        // No data → no BigWig file created
        assert!(!path.exists(), "No BigWig should be created for empty data");
    }

    #[test]
    fn test_zero_coverage_contig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let ud = vec![0i32; 100];

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("uncovered", &ud);
        writer.finish();

        assert!(path.exists());
        let _intervals = read_bigwig_intervals(&path);
        // All zero-depth — BigWig stores this as a single 0-value interval or no intervals
        // depending on implementation. Just check the file is valid.
        let reader = bigtools::BigWigRead::open_file(&path).unwrap();
        let chroms = reader.chroms().to_vec();
        assert_eq!(chroms.len(), 1);
        assert_eq!(chroms[0].name, "uncovered");
    }

    #[test]
    #[should_panic(expected = "already exists")]
    fn test_collision_detection() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let writer1 = CoverageProfileWriter::new(&path);
        writer1.finish();

        // This would try to create a temp file, but the .bw exists from
        // the first empty write... actually it won't exist (empty data).
        // Let's create a real file first.
        std::fs::write(&path, b"dummy").unwrap();
        let _writer2 = CoverageProfileWriter::new(&path);
    }

    #[test]
    fn test_overlapping_reads_depth_profile() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        // 3 overlapping reads: staircase pattern
        let ud = make_ups_and_downs(50, &[(10, 30, 1), (20, 40, 1), (25, 35, 1)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        let intervals = read_bigwig_intervals(&path);
        let nonzero: Vec<_> = intervals.iter().filter(|i| i.3 > 0.0).collect();
        // Should have depth 1, 2, 3, 2, 1 regions
        assert!(
            nonzero.len() >= 5,
            "Should have at least 5 nonzero regions, got {:?}",
            nonzero
        );
    }

    #[test]
    fn test_bigwig_file_is_valid() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let ud = make_ups_and_downs(1000, &[(10, 500, 5), (600, 900, 3)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("chr1", &ud);
        writer.finish();

        // Verify file is a valid BigWig
        let metadata = std::fs::metadata(&path).unwrap();
        assert!(metadata.len() > 0, "BigWig file should not be empty");

        let reader = bigtools::BigWigRead::open_file(&path).unwrap();
        let chroms = reader.chroms().to_vec();
        assert_eq!(chroms.len(), 1);
        assert_eq!(chroms[0].name, "chr1");
        assert_eq!(chroms[0].length, 1000);
    }

    #[test]
    fn test_no_temp_file_left() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");
        let temp_path = dir.path().join("test.bedgraph.tmp");

        let ud = make_ups_and_downs(100, &[(10, 90, 5)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        assert!(path.exists(), "BigWig should exist");
        assert!(!temp_path.exists(), "Temp BedGraph should be cleaned up");
    }

    #[test]
    fn test_contigs_in_non_lexicographic_order() {
        // Simulate BAM tid order: contigs are NOT in lexicographic order.
        // This is the real-world case for multi-contig MAGs.
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let mut writer = CoverageProfileWriter::new(&path);
        // Write contigs in non-lexicographic order (as a BAM would)
        writer.write_contig("contig_10", &make_ups_and_downs(100, &[(10, 90, 5)]));
        writer.write_contig("contig_2", &make_ups_and_downs(200, &[(0, 200, 3)]));
        writer.write_contig("contig_1", &make_ups_and_downs(150, &[(50, 100, 8)]));
        writer.write_contig("contig_20", &make_ups_and_downs(80, &[(10, 70, 2)]));
        writer.finish();

        assert!(
            path.exists(),
            "BigWig should be created despite non-lexicographic order"
        );

        // Verify all 4 contigs are present and readable
        let reader = bigtools::BigWigRead::open_file(&path).unwrap();
        let chroms = reader.chroms().to_vec();
        assert_eq!(chroms.len(), 4, "All 4 contigs should be in BigWig");

        let chrom_names: Vec<_> = chroms.iter().map(|c| c.name.clone()).collect();
        assert!(chrom_names.contains(&"contig_1".to_string()));
        assert!(chrom_names.contains(&"contig_2".to_string()));
        assert!(chrom_names.contains(&"contig_10".to_string()));
        assert!(chrom_names.contains(&"contig_20".to_string()));
    }

    #[test]
    fn test_many_contigs_non_lexicographic() {
        // Simulate a real MAG with 50 contigs in BAM tid order (numeric, not lexicographic)
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bw");

        let mut writer = CoverageProfileWriter::new(&path);
        for i in 0..50 {
            // BAM tid order: 0, 1, 2, ..., 49
            // Lexicographic order would be: 0, 1, 10, 11, ..., 19, 2, 20, ...
            let name = format!("scaffold_{}", i);
            let ud = make_ups_and_downs(500, &[(i * 5, i * 5 + 100, 3)]);
            writer.write_contig(&name, &ud);
        }
        writer.finish();

        assert!(path.exists());

        let reader = bigtools::BigWigRead::open_file(&path).unwrap();
        let chroms = reader.chroms().to_vec();
        assert_eq!(chroms.len(), 50, "All 50 contigs should be in BigWig");
    }
}
