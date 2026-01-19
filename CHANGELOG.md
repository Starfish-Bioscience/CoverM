# Changelog

All notable changes to CoverM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
