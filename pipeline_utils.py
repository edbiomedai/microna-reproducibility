import anndata as ad
import numpy as np
import pandas as pd
import scmorph as sm
from natsort import natsort_keygen

LOPACTARGETS = "/exports/igmm/datastore/khamseh-lab/jwagner/projects/micromorph/data/LOPAC_Targets_clean.csv"


def read_platemap(file, organism="mouse", lopac=LOPACTARGETS):
    """Read plate layout maps and name columns to match single-cell metadata"""
    df = pd.read_csv(file)

    # Filter organism
    df.query(f"organism == '{organism}'", inplace=True)

    # Adhere to naming convention
    df["PlateLayout"] = df["plate"].str.replace(f"{organism}_", "P")

    # Rename columns to match metadata
    df.rename(columns={"well": "Well"}, inplace=True)

    df = df[["PlateLayout", "Well", "treatment", "ttype"]].copy()

    if "ttype" in df.columns and "tr_dmso" in df["ttype"].unique():
        df["ttype"] = df["ttype"].astype(str)
        df.loc[df["ttype"] == "tr_dmso", "ttype"] = "DMSO"
        df["ttype"] = df["ttype"].astype("category")

    if "treatment" in df.columns and "control_tr_dmso" in df["treatment"].unique():
        df["treatment"] = df["treatment"].astype(str)
        df.loc[df["treatment"] == "control_tr_dmso", "treatment"] = "DMSO"
        df["treatment"] = df["treatment"].astype("category")

    return df

def read_adata(file, backed=False):
    """Read in single-cell data and rename columns to match metadata"""
    adata = ad.read_h5ad(file, backed="r") if backed else ad.read_h5ad(file)
    if "Experiment" in adata.obs.columns:
        adata.obs.rename(columns={"Experiment": "PlateLayout"})
    return adata


def process_metadata(adata, cellline=""):
    """Massage metadata to adhere to a common standard"""
    obs = adata.obs
    METAREGEX = r"CS-[A-Z0-9]{3}-R(?P<Replicate>[12])-(?P<PlateLayout>.*)_(?P<CellLine>.*)_(?P<Date>.*)_(?P<PlateID>[0-9]*)_(?P<Well>[A-Z][0-9]{2})_(?P<Site>[0-9])"

    st = obs["SampleName"].str

    meta = st.extract(METAREGEX, expand=True)
    meta.insert(0, "SampleName", obs["SampleName"])
    meta["CellLine"] = meta["CellLine"].astype(str).astype("category")
    meta["PlateLayout"] = meta["PlateLayout"].astype(str).astype("category")
    meta["CellLine"] = meta["CellLine"].astype(str).astype("category")
    meta["PlateID"] = meta["PlateID"].astype(int).astype("category")
    meta["Well"] = meta["Well"].astype(str).astype("category")
    meta["Site"] = meta["Site"].astype(int)

    if "Replicate" in meta.columns:
        meta["Replicate"] = meta["Replicate"].astype(int)

    df = obs.drop("SampleName", axis=1)

    new = pd.concat([meta, df], axis=1)

    adata.obs = new
    if not adata[0].obs["SampleName"].str.startswith("E")[0]:
        # Filter out control plates
        adata = adata[~adata.obs["PlateLayout"].str.startswith("C")].copy()
    return adata


def add_platemap(adata, platemap):
    """Join platemap to single-cell data"""
    new = adata.obs.merge(platemap, on=["PlateLayout", "Well"], how="left")
    new.index = adata.obs.index
    adata.obs = new
    return adata


def sort_adata(adata):
    """Sort AnnData by metadata columns"""
    idx = adata.obs.sort_values(
        ["CellLine", "PlateLayout", "PlateID", "Well", "Site"],
        key=natsort_keygen(),
    ).index
    return adata[idx].copy()  # avoid view due to slow IO

def filter_adata(scadata, qcadata):
    return (
        scadata[scadata.obs["PassQC"] == "True"].copy(),
        qcadata[qcadata.obs["PassQC"] == "True"].copy(),
    )

def read_miRNA_qc(qc_file, cellline=""):
    """miRNA-specific QC parser"""
    qc = ad.read_h5ad(qc_file)
    qc = process_metadata(qc, cellline=cellline)
    qc = sort_adata(qc)
    return qc


def count_cells_per_img(adata):
    """Helper to count (rounded) number of cells per image"""
    obs = adata.obs
    counts = (
        obs.groupby(["PlateID", "Well", "treatment"])
        .apply(lambda x: np.rint(len(x) / len(x["Site"].unique())))
        .astype(int)
        .rename("CellsPerImage")
        .to_frame()
        .reset_index()
    )
    return counts


def batch_correct(adata, copy=False, neg_control="DMSO"):
    assert neg_control in adata.obs["treatment"].unique()
    if copy:
        adata = sm.pp.remove_batch_effects(
            adata,
            batch_key="PlateID",
            treatment_key="treatment",
            control=neg_control,
            copy=copy,
        )
    else:
        sm.pp.remove_batch_effects(
            adata,
            batch_key="PlateID",
            treatment_key="treatment",
            control=neg_control,
            copy=copy,
        )

    adata.uns["batch_effects"].columns = adata.uns["batch_effects"].columns.astype(str)
    return adata
