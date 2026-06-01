import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotnine import *
from matplotlib import colors as mcolors
from anndata import AnnData

def get_data(
    adata,
    layer="X_pca",
    obs_cols=["Replicate", "PlateLayout", "Treatment", "Row", "Column"],
):
    def get_layer(adata, layer):
        if layer == "X":
            X = adata.to_df()
        else:
            X = adata.obsm[layer]
            column_names = [f"{layer}_{i+1}" for i in range(X.shape[1])]
            if layer == "X_umap":
                column_names = [f"UMAP_{i+1}" for i in range(X.shape[1])]
            else:
                column_names = [
                    f"PC{i + 1} ({round(v * 100, 2)}%)"
                    for i, v in enumerate(adata.uns["pca"]["variance_ratio"])
                ]
            X = pd.DataFrame(X, index=adata.obs_names, columns=column_names)
        return X

    X = get_layer(adata, layer)
    obs_cols = [col for col in obs_cols if col in adata.obs.columns]
    df = pd.concat([adata.obs[obs_cols], pd.DataFrame(X)], axis=1)
    return df


def plot_latent(adata, n_dim=(1, 2), color="Replicate", layer="X_pca", draw=False):
    df = get_data(adata, layer=layer, obs_cols=[color])
    x_col = df.columns[1]
    y_col = df.columns[2]
    fig = (
        ggplot(df, aes(x=x_col, y=y_col, color=color))
        + geom_point()
        + theme_bw()
        + theme(  # legend_position=(0.8,0.8),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            axis_line_x=element_line(size=1),
            axis_line_y=element_line(size=1),
            axis_text=element_blank(),
            #   aspect_ratio=1,
            dpi=150,
            #   figure_size=(5,5),
            legend_box="vertical",
            axis_ticks=element_blank(),
            legend_key=element_blank(),
        )
    )
    if f"{color}_colors" in adata.uns:
        fig = fig + scale_color_manual(values=adata.uns[f"{color}_colors"])
    if draw:
        fig.draw()
    return fig


def plot_latent_dims(
    adata, obs_cols=["CellLine", "Replicate", "PlateLayout", "Treatment"], layer="X_pca"
):
    map_cell_lines(adata)
    color_factor(adata)
    figs = []
    for c in obs_cols:
        figs.append(plot_latent(adata, color=c, layer=layer, draw=True))
    return figs


def map_cell_lines(x):
    mapping = {
        "C2C12": "C2C12",
        "N2As": "N2a",
        "Control Testing": "HEK293T",
        "Hepg2": "HepG2",
        "SHSY5Y": "SH-SY5Y",
    }
    if isinstance(x, AnnData):
        return map_cell_lines(x.obs)

    if not isinstance(x, pd.DataFrame):
        raise ValueError(
            f"Expected AnnData or DataFrame, got {type(x)}"
        )

    if "CellLine" not in x:
        return
    if "SH-SY5Y" not in x["CellLine"].values:
        cline = x["CellLine"]
        x["CellLine"] = cline.map(mapping, na_action=None)
    if not isinstance(x["CellLine"], pd.CategoricalDtype):
        x["CellLine"] = x["CellLine"].astype("category")

    x["CellLine"] = (
        x["CellLine"]
        .cat
        .set_categories(list(mapping.values()), ordered=True)
    )


def color_factor(adata):
    cell_line_colors(adata)
    plate_colors(adata)
    replicate_colors(adata)


def replicate_colors(adata):
    adata.uns["Replicate_colors"] = ["#000000", "#BB09BE"]
    return


def plate_colors(adata):
    from natsort import natsort_key

    plates = adata.obs["PlateLayout"].unique()
    plates = sorted(plates, key=natsort_key)
    colors = plt.cm.viridis(np.linspace(0, 1, len(plates)))
    hex_colors = [mcolors.rgb2hex(color) for color in colors]
    plate_colors = dict(zip(plates, hex_colors))
    adata.obs["PlateLayout"] = (
        adata.obs["PlateLayout"]
        .astype("category")
        .cat.set_categories(plates, ordered=True)
    )
    adata.uns["PlateLayout_colors"] = [plate_colors[plate] for plate in plates]
    return


def cell_line_colors(adata):
    mapping = {
        "C2C12": "blue",
        "N2a": "#61E8FF",
        "HEK293T": "#FD8008",
        "SH-SY5Y": "#804003",
        "HepG2": "#FECC66",
    }
    adata.uns["CellLine_colors"] = [
        mapping[cell] for cell in adata.obs["CellLine"].cat.categories
    ]
    return


def rorder(adata, fraction=1, copy=False):
    from scanpy.pp import subsample

    np.random.seed(2025)
    if copy:
        return subsample(adata, fraction=fraction, copy=copy)
    else:
        subsample(adata, fraction=fraction, copy=copy)
        return adata


def get_plates_with_treatments_of_interest(adata, search_treatments, keep_treatments=["control_tr_dmso"]):
    keep_treatments += search_treatments
    plates = []
    for cell_line in adata.obs["CellLine"].unique():
        p = (
            adata[
                (adata.obs["CellLine"] == cell_line)
                & (adata.obs["Treatment"].isin(search_treatments))
            ]
            .obs["PlateID"]
            .unique()
        )
        plates.extend(p)
    plates = list(set(plates))
    return adata[(adata.obs["PlateID"].isin(plates)) & (adata.obs["Treatment"].isin(keep_treatments))]
