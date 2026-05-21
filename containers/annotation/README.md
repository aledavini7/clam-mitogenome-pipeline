# CLAM Annotation Container

This image contains the lightweight runtime used for CLAM mitochondrial variant
summary and annotation helper scripts.

Current contents:

- Python 3.11

Build locally:

```bash
docker buildx build --platform linux/amd64 --load -t clam-annotation:0.1.0 containers/annotation
```

Suggested publish target for GitHub Container Registry:

```bash
docker buildx build \
  --platform linux/amd64 \
  --push \
  -t ghcr.io/aledavini7/clam-annotation:0.1.0 \
  -t ghcr.io/aledavini7/clam-annotation:latest \
  containers/annotation
```

Future revisions will add the vcf2maf/VEP runtime once the VEP cache strategy is
finalized for the cluster.
