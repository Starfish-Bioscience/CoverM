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
