# miRNA analysis reproducibility

This repository contains code to reproduce findings from "Systematic
morphological profiling of human and mouse miRNAs in 24M single cells",
a preprint available at: https://www.biorxiv.org/content/10.1101/2025.11.14.687149v1

It utilises `scmorph`, a Python package for analysis of morphological
profiles. The corresponding paper can be found here: https://doi.org/10.21105/joss.08324.

The notebooks herein aim at reproducibility, but are not an introduction to the
field or scmorph. For those, please see
https://doi.org/10.1038/s44320-026-00197-7 and https://scmorph.readthedocs.io/.

To get started with this repository, we first need to make sure that git, aws and
conda are available. Then, we will download the repository's contents, create a
conda environment, and download the required data.

What is available
---
Data is available in various stages of processing through the [Cell Painting
Gallery](https://broadinstitute.github.io/cellpainting-gallery/overview.html),
under the `cpg0046-microrna` project. Available data includes the folders:
- `images`: raw images as they come off the microscope.
- `workspace`: processed data, with the subfolders:
    - `metadata`: CSV files with platemaps per cell line and biological replicate.
    - `load_data_csv`: CSV files that can be used to load the data into CellProfiler.
    - `qc`: image-level QC features, used for QC. In
[AnnData](https://anndata.scverse.org/en/stable/)
format.
    - `profiles_assembled`, with the subfolders:
        - `profiles_raw`: Single-cell features, one AnnData file per cell line
        - `profiles_qc`: Single-cell features, one AnnData file per cell line
          and biological replicate. These features are after image-level quality
          control, and thus contain fewer cells than the data in `profiles_raw`.
          - `profiles_processed`: Single-cell features after QC and feature
            selection, ready to be used for hit calling.

Prerequisites
---
- Install miniconda, if you have not already:
  https://www.anaconda.com/download/success?reg=skipped
- If you are on Windows, set up git: https://git-scm.com/install/windows
- Install aws cli: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html

Setup
---
Open a terminal and download the repository's content:
```bash
git clone https://github.com/edbiomedai/microrna-reproducibility && cd micorna-reproducibility
```

Create a conda environment with required dependencies:
```bash
conda create -y -n micrornas -f env.yaml
```

Download data
---
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

You can also download the raw images, but they will not be used by the code in
this repository. If you want them anyway, you
can download them (note: large download, ~15 Tb) using:

```bash
aws s3 cp "s3://cellpainting-gallery/cpg0046-microrna/edinburgh_igc/images images --no-sign-request
```

Getting Started
---
The notebooks cover various processing steps. They can be run independently, but
the outputs will overwrite the data you downloaded from the Cell Painting Gallery,
unless you change the output path in the notebook.
- `00_imageQC.ipynb` covers the step from single-cell features + image-level QC
  features -> QC'ed single-cell features
- `01_feature_selection.ipynb` goes from QC'ed single-cell features -> feature
  selected single-cell features
- `02_rank_mirnas.ipynb` goes from feature
  selected single-cell features -> hits
- `03_data_access.ipynb`: this notebook is useful **if you are interested in using the
data for machine learning/AI**. It demonstrates how
`workspace/profiles_assembled/profiles_processed/all_adata_featFilt_batchCorr.h5ad`,
which contains single-cell features after batch correction and quality control,
can be read in and used.
