# CLAM Core Container

This image contains the runtime used by the CLAM mitochondrial pipeline core:

- GSNAP/GMAP for rCRS and nuclear-only remapping
- samtools and htslib/bgzip
- Java 17
- MitoScape 1.0 jar and bundled model/resources pinned from upstream commit `118d4d5dfcf66774296017130b8886eee0290f49`
- mutserve 2.0.3 with bundled rCRS resources
- HaploGrep3 3.2.2 with `phylotree-rcrs@17.2` installed at build time

Build locally:

```bash
docker buildx build --platform linux/amd64 --load -t clam-core:0.1.0 containers/clam-core
```

The Dockerfile defaults to `linux/amd64` because the cluster/Seqera runtime is
expected to run on Linux x86_64 nodes. On Apple Silicon, keep the explicit
`--platform linux/amd64` flag to avoid accidentally creating a local ARM image.

```bash
docker image inspect clam-core:0.1.0 --format '{{.Os}}/{{.Architecture}}'
```

Suggested publish target for GitHub Container Registry:

```bash
docker buildx build \
  --platform linux/amd64 \
  --push \
  -t ghcr.io/aledavini7/clam-core:0.1.0 \
  -t ghcr.io/aledavini7/clam-core:latest \
  containers/clam-core
```

The same `linux/amd64` build is automated by
`.github/workflows/build-clam-core.yml`.

The pipeline currently expects these paths inside the container:

- `/opt/clam/share/mitoscape/MitoScapeClassify.jar`
- `/opt/clam/share/mitoscape/mitomap.ld`
- `/opt/clam/share/mitoscape/MTClassifierModel.RF`
- `/opt/clam/share/mitoscape/NUMTs_hg38.txt`
- `/opt/clam/mutserve/rCRS.fasta`
- `/opt/clam/mutserve/rCRS_annotation_2020-08-20.txt`
- `/usr/local/bin/haplogrep3`
- `/usr/local/bin/mutserve`
