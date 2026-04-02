# CoverM — Benchmark Results

Dataset: CL4-AERO (ONT long-read, WGS21_Bact_Bioaster)  
Reference: sfbio-mag-catalogue v0.1.0 (`.mmi` pre-built index)  
BED file: `annot_regions.bed` — 19 M regions / 1.7 GB (labels: cog, housekeeping, rrna + unlabeled)  
Host: 56-core server  
Parameters: `--threads 52 --minimap2-reference-is-index --mapper minimap2-lr-hq`  
Date: 2026-04-02

## Results

| Version | Conditions | CPU% | Wall time | Peak RAM |
|---------|-----------|------|-----------|----------|
| v0.9.0-cram  | no BED (baseline) | 1946% | 1:22 | 4.1 GB |
| v0.10.1-spatial | no BED | 1880% | 1:25 | 4.1 GB |
| v0.11.0-bedcov-labels | no BED | 1927% | 1:23 | 4.1 GB |
| v0.11.0-bedcov-labels | + BED + unlabeled | 1427% | **1:57** | 6.2 GB |
| v0.11.0-bedcov-labels | + BED + unlabeled + bedgraph | 1397% | **1:59** | 6.8 GB |

## Notes

- No regression between v0.9, v0.10 and v0.11 on the base mapping path (~1:23).
- BED processing overhead on a 19 M-region file: +34 s (+41%), +2.1 GB RAM.
- bedGraph output adds ~2 s on top of BED processing — negligible.
- Peak RAM with BED is bounded by rayon thread buffers + interned chrom table;
  the BED file itself is memory-mapped (memmap2, added in v0.11.0) so it does
  not contribute a full heap copy.

## How to reproduce

```bash
# Single scenario
bash benchmark/run_bench.sh \
    --sif /HDD/singularity_apptainer/coverm_0.11.0-bedcov-labels.sif \
    --label v0.11.0-bedcov-labels \
    --scenarios bed

# All scenarios
bash benchmark/run_bench.sh \
    --sif /HDD/singularity_apptainer/coverm_0.11.0-bedcov-labels.sif
```
