# CLAM Mutect2 Container

This image wraps the pinned GATK4 Biocontainer used for the CLAM Mutect2 step:

- Base image: `quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0`
- Target platform: `linux/amd64`

Build locally:

```bash
docker buildx build --platform linux/amd64 --load -t clam-mutect2:0.1.0 containers/mutect2
```

Publish to GitHub Container Registry:

```bash
docker buildx build \
  --platform linux/amd64 \
  --push \
  -t ghcr.io/aledavini7/clam-mutect2:0.1.0 \
  -t ghcr.io/aledavini7/clam-mutect2:latest \
  containers/mutect2
```
