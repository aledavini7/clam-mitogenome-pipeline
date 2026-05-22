# CLAM MitoScape Container

This image contains the Java 8 runtime used by the MitoScape NUMT
classification step.

MitoScape 1.0 bundles an older Spark/Scala runtime that is not fully compatible
with Java 17. The rest of CLAM can continue to use `clam-core`, while this
process runs in a focused Java 8 container.

Build locally:

```bash
docker buildx build --platform linux/amd64 --load -t clam-mitoscape:0.1.0 containers/mitoscape
```

Suggested publish target for GitHub Container Registry:

```bash
docker buildx build \
  --platform linux/amd64 \
  --push \
  -t ghcr.io/aledavini7/clam-mitoscape:0.1.0 \
  -t ghcr.io/aledavini7/clam-mitoscape:latest \
  containers/mitoscape
```

The pipeline expects these paths inside the container:

- `/opt/clam/share/mitoscape/MitoScapeClassify.jar`
- `/opt/clam/share/mitoscape/mitomap.ld`
- `/opt/clam/share/mitoscape/NUMTs_hg38.txt`
- `/opt/clam/share/mitoscape/MTClassifierModel.RF`
