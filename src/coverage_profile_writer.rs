use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use rust_htslib::bgzf;

/// Writes per-base coverage profiles in BedGraph format (bgzf compressed).
///
/// The BedGraph writer produces a raw, unfiltered coverage profile — no
/// contig_end_exclusion, no min_island_length. It's intentionally decoupled
/// from the spatial metrics estimators: the BedGraph is primary data, the
/// metrics are derived statistics.
///
/// Note: this writer performs its own cumulative sum over ups_and_downs,
/// independently from the estimators. This is intentional — the bottleneck
/// is BAM I/O, not an in-memory vector scan.
pub struct CoverageProfileWriter {
    writer: bgzf::Writer,
    path: PathBuf,
}

impl CoverageProfileWriter {
    /// Create a new BedGraph writer at the given path.
    ///
    /// The file is written in bgzf (block-gzip) format for compatibility
    /// with tabix, IGV, JBrowse, bedtools, etc.
    ///
    /// Panics if the file already exists (no silent overwrite).
    pub fn new(path: &Path) -> Self {
        if path.exists() {
            panic!(
                "Coverage profile output file already exists: {}. \
                 Refusing to overwrite.",
                path.display()
            );
        }
        let writer = bgzf::Writer::from_path(path).unwrap_or_else(|e| {
            panic!("Failed to create bgzf writer at {}: {}", path.display(), e)
        });
        CoverageProfileWriter {
            writer,
            path: path.to_path_buf(),
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
                // Depth changed — write the previous segment
                writeln!(
                    self.writer,
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
            self.writer,
            "{}\t{}\t{}\t{}",
            contig_name, segment_start, len, segment_depth
        )
        .unwrap_or_else(|e| panic!("Failed to write BedGraph: {}", e));
    }

    /// Close the bgzf writer and create a tabix index (.tbi).
    ///
    /// The tabix index enables random access by region (e.g. `tabix file.bedgraph.gz contig:start-end`).
    pub fn finish(self) {
        // Drop the writer to flush and close the bgzf file
        drop(self.writer);

        // Create tabix index
        let status = Command::new("tabix")
            .args(["-p", "bed", self.path.to_str().unwrap()])
            .status();

        match status {
            Ok(s) if s.success() => {
                info!("Created tabix index for {}", self.path.display());
            }
            Ok(s) => {
                warn!(
                    "tabix indexing failed with status {} for {}. \
                     The BedGraph file is still valid but not indexed. \
                     Run 'tabix -p bed {}' manually.",
                    s,
                    self.path.display(),
                    self.path.display()
                );
            }
            Err(e) => {
                warn!(
                    "Could not run tabix: {}. \
                     The BedGraph file is still valid but not indexed. \
                     Install tabix (part of htslib/samtools) and run \
                     'tabix -p bed {}' manually.",
                    e,
                    self.path.display()
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::TempDir;

    /// Helper: read a bgzf file back to a String
    fn read_bgzf(path: &Path) -> String {
        let mut reader = bgzf::Reader::from_path(path).unwrap();
        let mut contents = String::new();
        reader.read_to_string(&mut contents).unwrap();
        contents
    }

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

    #[test]
    fn test_uniform_coverage_bedgraph() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud = make_ups_and_downs(100, &[(0, 100, 10)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("contig_1", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        // Uniform 10x → single segment
        assert_eq!(contents, "contig_1\t0\t100\t10\n");
    }

    #[test]
    fn test_rle_compression() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        // 100bp contig: [0..20) at 0x, [20..60) at 5x, [60..80) at 12x, [80..100) at 0x
        let ud = make_ups_and_downs(100, &[(20, 60, 5), (60, 80, 12)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        let expected = "ctg\t0\t20\t0\nctg\t20\t60\t5\nctg\t60\t80\t12\nctg\t80\t100\t0\n";
        assert_eq!(contents, expected);
    }

    #[test]
    fn test_multiple_contigs() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud1 = make_ups_and_downs(50, &[(10, 40, 8)]);
        let ud2 = make_ups_and_downs(30, &[(0, 30, 3)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("contig_A", &ud1);
        writer.write_contig("contig_B", &ud2);
        writer.finish();

        let contents = read_bgzf(&path);
        assert!(contents.contains("contig_A\t0\t10\t0\n"));
        assert!(contents.contains("contig_A\t10\t40\t8\n"));
        assert!(contents.contains("contig_A\t40\t50\t0\n"));
        assert!(contents.contains("contig_B\t0\t30\t3\n"));
    }

    #[test]
    fn test_empty_contig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("empty", &[]);
        writer.finish();

        let contents = read_bgzf(&path);
        assert_eq!(contents, "");
    }

    #[test]
    fn test_zero_coverage_contig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud = vec![0i32; 100];

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("uncovered", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        // Entire contig at 0x → single segment
        assert_eq!(contents, "uncovered\t0\t100\t0\n");
    }

    #[test]
    #[should_panic(expected = "already exists")]
    fn test_collision_detection() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        // Create first writer and finish it
        let writer1 = CoverageProfileWriter::new(&path);
        writer1.finish();

        // Second writer should panic
        let _writer2 = CoverageProfileWriter::new(&path);
    }

    #[test]
    fn test_tabix_index_created() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");
        let tbi_path = dir.path().join("test.bedgraph.gz.tbi");

        let ud = make_ups_and_downs(100, &[(10, 50, 5)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("chr1", &ud);
        writer.finish();

        // Check that tabix index was created (only if tabix is available on PATH)
        if std::process::Command::new("tabix")
            .arg("--version")
            .output()
            .is_ok()
        {
            assert!(tbi_path.exists(), "tabix index file should exist");
        }
    }

    // =========================================================================
    // Edge case tests for BedGraph writer
    // =========================================================================

    #[test]
    fn test_overlapping_reads_depth_profile() {
        // 3 overlapping reads creating a staircase pattern
        // Read 1: [10..30), Read 2: [20..40), Read 3: [25..35)
        // Depth: [0..10)=0, [10..20)=1, [20..25)=2, [25..30)=3, [30..35)=2, [35..40)=1, [40..50)=0
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud = make_ups_and_downs(50, &[(10, 30, 1), (20, 40, 1), (25, 35, 1)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        let lines: Vec<&str> = contents.trim().split('\n').collect();
        assert_eq!(lines.len(), 7, "Should have 7 RLE segments");
        assert_eq!(lines[0], "ctg\t0\t10\t0");
        assert_eq!(lines[1], "ctg\t10\t20\t1");
        assert_eq!(lines[2], "ctg\t20\t25\t2");
        assert_eq!(lines[3], "ctg\t25\t30\t3");
        assert_eq!(lines[4], "ctg\t30\t35\t2");
        assert_eq!(lines[5], "ctg\t35\t40\t1");
        assert_eq!(lines[6], "ctg\t40\t50\t0");
    }

    #[test]
    fn test_bedgraph_segments_cover_full_contig_length() {
        // Verify that the sum of all segment lengths equals the contig length
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let contig_len = 500;
        let ud = make_ups_and_downs(contig_len, &[(10, 50, 5), (100, 200, 3), (300, 450, 8)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        let total_len: usize = contents
            .trim()
            .lines()
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                let start: usize = fields[1].parse().unwrap();
                let end: usize = fields[2].parse().unwrap();
                end - start
            })
            .sum();
        assert_eq!(
            total_len, contig_len,
            "Sum of segment lengths should equal contig length"
        );
    }

    #[test]
    fn test_bedgraph_segments_are_contiguous() {
        // Verify that segments are contiguous (no gaps, no overlaps)
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud = make_ups_and_downs(200, &[(20, 80, 5), (60, 150, 3)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        let mut prev_end: Option<usize> = None;
        for line in contents.trim().lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            let start: usize = fields[1].parse().unwrap();
            let end: usize = fields[2].parse().unwrap();

            assert!(end > start, "Segment end must be > start");
            if let Some(pe) = prev_end {
                assert_eq!(start, pe, "Segments must be contiguous");
            }
            prev_end = Some(end);
        }
    }

    #[test]
    fn test_bedgraph_no_adjacent_same_depth() {
        // Verify that RLE is correct: no two adjacent segments with same depth
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud = make_ups_and_downs(300, &[(10, 100, 5), (100, 200, 5), (200, 290, 5)]);
        // All at same depth 5 from 10..290, should be merged

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("ctg", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        let lines: Vec<&str> = contents.trim().split('\n').collect();
        // Should be: 0..10 at 0, 10..290 at 5, 290..300 at 0 → 3 segments
        assert_eq!(
            lines.len(),
            3,
            "Adjacent same-depth should be merged: {:?}",
            lines
        );
    }

    #[test]
    fn test_single_base_contig() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud = vec![5i32; 1]; // 1bp contig at depth 5

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("tiny", &ud);
        writer.finish();

        let contents = read_bgzf(&path);
        assert_eq!(contents, "tiny\t0\t1\t5\n");
    }

    #[test]
    fn test_many_contigs_ordering() {
        // Write many contigs and verify ordering is preserved
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let mut writer = CoverageProfileWriter::new(&path);
        for i in 0..20 {
            let name = format!("contig_{:03}", i);
            let ud = make_ups_and_downs(100, &[(10, 90, (i + 1) as i32)]);
            writer.write_contig(&name, &ud);
        }
        writer.finish();

        let contents = read_bgzf(&path);
        let contig_names: Vec<&str> = contents
            .lines()
            .map(|l| l.split('\t').next().unwrap())
            .collect();
        // First contig_000, then contig_001, etc.
        assert!(contig_names.starts_with(&["contig_000"]));
        assert!(contig_names.contains(&"contig_019"));
    }

    #[test]
    fn test_mixed_covered_uncovered_contigs() {
        // Simulate a real scenario: some contigs covered, some not
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let mut writer = CoverageProfileWriter::new(&path);

        // Covered contig
        let ud1 = make_ups_and_downs(200, &[(10, 190, 8)]);
        writer.write_contig("covered_ctg", &ud1);

        // Uncovered contig (all zeros)
        let ud2 = vec![0i32; 500];
        writer.write_contig("uncovered_ctg", &ud2);

        // Another covered contig
        let ud3 = make_ups_and_downs(100, &[(0, 100, 3)]);
        writer.write_contig("also_covered", &ud3);

        writer.finish();

        let contents = read_bgzf(&path);
        // Uncovered contig should still appear as single 0-depth segment
        assert!(contents.contains("uncovered_ctg\t0\t500\t0\n"));
        assert!(contents.contains("covered_ctg\t"));
        assert!(contents.contains("also_covered\t"));
    }

    #[test]
    fn test_bgzf_file_readable_by_zcat() {
        // Verify the output is valid bgzf (readable by bgzf::Reader and standard tools)
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.bedgraph.gz");

        let ud = make_ups_and_downs(100, &[(10, 90, 5)]);

        let mut writer = CoverageProfileWriter::new(&path);
        writer.write_contig("test", &ud);
        writer.finish();

        // Verify file exists and is not empty
        let metadata = std::fs::metadata(&path).unwrap();
        assert!(metadata.len() > 0, "bgzf file should not be empty");

        // Verify it's readable as bgzf
        let contents = read_bgzf(&path);
        assert!(contents.contains("test\t10\t90\t5"));
    }
}
