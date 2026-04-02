# Changelog

All notable changes to CoverM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.11.0] - 2026-04-02

### Added

#### BED region coverage (`genome` subcommand)

- `--regions-bed <FILE>`: load a 4-column BED file (`chrom start end label`) to
  compute per-region coverage statistics grouped by label.
  Appends `NB_LABELS × NB_METHODS` extra columns to the main output table,
  one per (label, method) combination, using the same `--methods` as the
  genome-level output.
- `--output-bedcov <DIR>`: write one multi-track bedGraph file per sample into
  DIR (`{DIR}/{sample_id}.bedgraph`). Each file contains one
  `track type=bedGraph` per label (alphabetical order), with mean coverage as
  the data value. Compatible with IGV, trackViewer (R/Bioconductor), and
  pyGenomeTracks. Note: `rtracklayer::import.bedGraph()` does not support
  multi-track files — use trackViewer instead.
- `--output-bedcov-compress`: gzip-compress bedGraph output files
  (`.bedgraph.gz`). Requires `--output-bedcov`.

### Changed

- `--contig-end-exclusion` is now incompatible with `--regions-bed` in the
  `genome` subcommand (runtime error if both are provided explicitly).

### Fixed

- `ReadsPerBaseCalculator`: returns `0.0` instead of `NaN` when called on a
  zero-length genome (previously caused a 0/0 division).

---

> **Pending integration from parallel branches (resolve at merge time):**
>
> #### Spatial coverage metrics — `feature/spatial-metrics` (v0.10.x)
> - `islands_per_mbp`, `max_gap`, `gap_fraction` spatial coverage metrics
> - `--min-island-length` option
> - `--coverage-profile <dir>`: per-sample BigWig output
>
> #### CRAM support — `develop-cram` (v0.9.0)
> - `--use-cram`: write cached alignment files in CRAM format instead of BAM
>   (reduces disk usage ~50%). Available in `contig`, `genome`, and `make`.

---

## [0.8.0] - 2026-01-19

### Added
- ANIr coverage method for average nucleotide identity estimation (#258)
- `minimap2-lr-hq` mapping preset option for high-quality long reads
- `strobealign-aemb` method for contig-level mapping
- Option to name cached BAM outputs with `--cache-unfiltered-bam-files`
- Support for macOS ARM64 (osx-arm64) platform
- Pixi package manager support with proper project configuration
- Demo content and walkthrough documentation

### Changed
- Renamed `--bam-file-cache-directory` to `--cache-unfiltered-bam-directory` (old flag kept as alias for backward compatibility) (#266)
- Migrated build system to use Pixi for dependency management
- Updated Galah dependency to include contig clustering support
- Restricted minimap2 to versions below v2.29 for compatibility
- Improved genome mode to not require checkm-tab-table when using min-completeness

### Fixed
- Large header BAM inputs handling in filter mode
- Tolerance adjustments for zero-coverage tests
- Various clippy warnings and code formatting improvements

### Documentation
- Updated README with run dev command using release mode
- Added AGENTS.md for development workflow
- Clarified that RPKM/TPM methods are for contigs/genomes, not genes (#248)
- Added badges, DOI, and arXiv citation information
- Improved help messages for covered_fraction and strobealign options

## [0.7.0] - 2023-XX-XX

Previous stable release. See git history for details.
