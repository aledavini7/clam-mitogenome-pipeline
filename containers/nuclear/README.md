# Nuclear WES runtime

This image contains pinned versions of BWA-MEM2, samtools, CNVkit, Ensembl VEP,
and Python for the independent nuclear WES branch. Large reference files,
capture BEDs, VEP caches, population VCFs, and COSMIC resources remain external
workflow inputs and are not included in the image.

Build locally with:

```bash
docker buildx build --platform linux/amd64 --load -t clam-nuclear:0.1.0 containers/nuclear
```

The offline VEP cache must match GRCh38 and VEP cache release 112 used by this
image. Optional gnomAD/COSMIC VCFs must use compatible GRCh38 coordinates.
