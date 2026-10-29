# CLAM Manta container

This image uses the Bioconda Manta 1.6.0/Python 2.7 runtime and adds `bgzip`
and `tabix` from an isolated HTSlib environment for preparing Manta's indexed
call-region BED. Manta is kept in a separate image because its runtime is
separate from the CLAM nuclear tools.

The workflow configures Manta with `--exome` and call regions covering the
primary GRCh38 chromosomes. The broad regions avoid limiting breakpoints to
captured exons, while excluding decoy and unplaced contigs that can add runtime.

The image is built and published as:

```text
ghcr.io/aledavini7/clam-manta:0.1.0
```
