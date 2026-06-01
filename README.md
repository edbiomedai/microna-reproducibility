
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

Data is available for five cell lines: C2C12, N2a, HEK293T, HepG2,
and SH-SY5Y. You can download data from [CellPainting
Gallery](https://broadinstitute.github.io/cellpainting-gallery/overview.html).
See [this
page](https://broadinstitute.github.io/cellpainting-gallery/download_instructions.html)
for more information on how to download from the CellPainting Gallery. If you
wish to follow all notebooks in this repository, you will need the whole
`workspace` folder. Therefore, download it using

```bash
aws s3 cp "s3://cellpainting-gallery/cpg0046-microrna/edinburgh_igc/workspace workspace --no-sign-request
```

This will download the various image-level and single-cell datasets that are
used in the analysis.

