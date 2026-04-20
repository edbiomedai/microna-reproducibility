
Compendium notebooks for preprint found at https://doi.org/10.21105/joss.08324

Setup
---
Download the repository's content
```bash
git clone https://github.com/edbiomedai/microrna-reproducibility && cd micorna-reproducibility
```

Create a conda environment with required dependencies
```bash
conda create -y -n micrornas -f env.yaml
```

Preprocessed data is available for five cell lines: C2C12, N2a, HEK293T, HepG2,
and SH-SY5Y. You can download data from [CellPainting
Gallery](https://broadinstitute.github.io/cellpainting-gallery/overview.html).
See [this
page](https://broadinstitute.github.io/cellpainting-gallery/download_instructions.html)
for more information on how to download from the CellPainting Gallery:
```bash
CELLLINE="C2C12"
aws s3 cp "s3://cellpainting-gallery/cpg0046-microrna/edinburgh_igc/workspace/profiles_assembled/${CELLLINE}_allPlates_postQC_featureSelected.h5ad" profiles/ --no-sign-request
```

This will download the preprocessed, single-cell morphological profiles of one
cell line. To get data for the other cell lines, just change the line
`CELLLINE="C2C12"` to e.g. `CELLLINE="N2a"`.

