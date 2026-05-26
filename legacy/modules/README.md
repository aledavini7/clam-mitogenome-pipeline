# Legacy Modules

These modules are not imported by the active `clam.nf` workflow.

They are preserved as historical CLAM reference material while the pipeline is
modernized. Useful ideas can be reintroduced as clean, tested modules later, but
new development should happen in `modules/`.

Archived modules:

- `9_vcf2maf.nf`: older local vcf2maf/VEP-based MAF conversion.
- `10_numt_exploration.nf`: experimental NUMT-read exploration workflow.
- `11_numt_read_count.nf`: chromosome-level read-count exploration for NUMT reads.
- `12_wxs_cell_lines.nf`: older WXS/cell-line before/after calling helpers.
- `7_mutserve_annotation_legacy.nf`: mutserve annotation helpers removed from the
  active variant-calling module.
