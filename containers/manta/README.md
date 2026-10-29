# CLAM Manta container

This image contains Manta 1.6.0 and HTSlib for optional structural-variant
calling on nuclear WES BAMs. Manta is kept in a separate image because its
Bioconda build uses a Python 2.7 runtime, separate from the CLAM nuclear tools.

The workflow configures Manta with `--exome` and call regions covering the
primary GRCh38 chromosomes. The broad regions avoid limiting breakpoints to
captured exons, while excluding decoy and unplaced contigs that can add runtime.

The image is built and published as:

```text
ghcr.io/aledavini7/clam-manta:0.1.0
```
