# -*- coding: utf-8 -*-
import dash
from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Input, Output
from dash import no_update
import plotly.graph_objs as go

# from plotly import tools
from plotly import subplots

import numpy as np
import pandas as pd
import uuid
import colorlover as cl
import feather
import re
from sklearn.preprocessing import scale, minmax_scale

import annotations as an

from flask_caching import Cache

import matplotlib

matplotlib.use("Agg")
import seaborn as sns

import logging

# Configure logging
logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# fix back once done testing
# data_dir = '/efs/'
data_dir = "assets"


# Create append
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "CD4 Differentiation Trajectory Browser"
application = app.server

tab_selected_style = {"padding": "4px"}


# Cache
cache = Cache(
    app.server,
    config={
        "CACHE_TYPE": "redis",
        # Note that filesystem cache doesn't work on systems with ephemeral
        # filesystems like Heroku.
        "CACHE_TYPE": "filesystem",
        "CACHE_DIR": "cache-directory",
        # should be equal to maximum number of users on the app at a single time
        # higher numbers will store more data in the filesystem / redis de
        "CACHE_THRESHOLD": 200,
    },
)

# CSS
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True


# High resolution image
def HighResolutionGraph(**kwargs):
    return dcc.Graph(
        config={"toImageButtonOptions": {"filename": "treg_dep_image", "scale": 5}},
        **kwargs,
    )


app.layout = html.Div(
    children=[
        html.Link(rel="stylesheet", href="/assets/style.css"),
        # UUid
        html.Div(str(uuid.uuid4()), id="session-id", style={"display": "none"}),
        # Title
        html.H2(
            children="Differentiation of precursor central memory cells from heterogeneous naïve CD4+ T cells",
            style={"text-align": "center"},
        ),
        # Sidebar
        html.Div(
            children=[
                html.H6("Dataset"),
                dcc.Dropdown(id="dataset-dropdown", value="d1d2_normlog"),
                html.Br(),
                html.H6("Gene"),
                dcc.Dropdown(id="gene-name", value="MX1"),
                html.Br(),
                html.H6("Multiple genes"),
                dcc.Dropdown(
                    id="multi-gene-name", value="ISG15", disabled=True, multi=True
                ),
                html.Br(),
                dcc.Markdown(
                    """
`5/1/2023`: Data browser is now live!
                    """
                ),
            ],
            className="three columns",
        ),
        # Body
        html.Div(
            children=[
                dcc.Tabs(
                    id="tabs",
                    children=[
                        # Introduction tab
                        dcc.Tab(
                            label="Home",
                            children=[
                                html.Br(),
                                html.Img(
                                    src="assets/banner.jpg",
                                    style={"width": "80%", "align": "center"},
                                ),
                                html.Br(),
                                html.Br(),
                                dcc.Markdown(
                                    """
Browser description...
This site was developed using the [```dash```](https://dash.plot.ly) framework.

     """
                                ),
                                dcc.Markdown(
                                    """
##### Study summary

We ... # fill in
    """
                                ),
                                dcc.Markdown(
                                    """


        """
                                ),
                                dcc.Markdown(
                                    """
##### Datasets

The following datasets are available through this website:

**Dataset description...**


    """
                                ),
                                dcc.Markdown(
                                    """


        """
                                ),
                                dcc.Markdown(
                                    """

##### Tutorials


We recommend using the Google Chrome browser for optimal experience. See the ```Tutorial``` tab for demonstration of the functionality.


Navigate to the ```Gene expression``` tab and select the dataset of interest from the side bar.
Select a gene from corresponding drop down bar to visualize normalized expression on t-SNE and as violin plots. The expression levels displayed are library size normalized and log transformed (log (x + 1)).


```plot.ly``` is used as the plotting engine. Hover for the plot for various options including zooming and downloading the plot.
See the [plot.ly tutorial](https://help.plot.ly/zoom-pan-hover-controls/) for a detailed description of these options.

        """
                                ),
                                dcc.Markdown(
                                    """


        """
                                ),
                                dcc.Markdown(
                                    """


        """
                                ),
                                dcc.Markdown(
                                    """
##### Data access

Data is available in the `Downloads` tab. The following datasets are available
* Count matrix of all cells with metadata.
    """
                                ),
                                dcc.Markdown(
                                    """


        """
                                ),
                                dcc.Markdown(
                                    """
##### Citation

If you use the data, please cite our [manuscript](link to paper).
```
Differentiation of precursor central memory cells from heterogeneous naïve CD4+ T cells. 

Author list...




```

    """
                                ),
                                dcc.Markdown(
                                    """


        """
                                ),
                                html.Br(),
                                html.Br(),
                                # End tab
                            ],
                        ),
                        # Gene expression tab
                        dcc.Tab(
                            label="Gene expression",
                            id="gene-expression-tab",
                            children=[
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [
                                                dcc.RadioItems(
                                                    id="meta-field-1",
                                                    options=[
                                                        "cluster",
                                                        "condition",
                                                        "sample",
                                                    ],
                                                    value="cluster",
                                                    labelStyle={
                                                        "display": "inline-block",
                                                        "margin-right": "20px",
                                                        "margin-left": "5px",
                                                    },
                                                ),
                                                HighResolutionGraph(id="meta-graph-1"),
                                            ]
                                        ),
                                        dbc.Col(
                                            [
                                                dcc.RadioItems(
                                                    id="meta-field-2",
                                                    options=[
                                                        "cluster",
                                                        "condition",
                                                        "sample",
                                                    ],
                                                    value="condition",
                                                    labelStyle={
                                                        "display": "inline-block",
                                                        "margin-right": "20px",
                                                        "margin-left": "5px",
                                                    },
                                                ),
                                                HighResolutionGraph(id="meta-graph-2"),
                                            ]
                                        ),
                                        dbc.Col(
                                            [HighResolutionGraph(id="expression-graph")]
                                        ),
                                    ]
                                ),
                                html.H6("Cluster expression"),
                                HighResolutionGraph(id="cluster-expression"),
                                # End tab
                            ],
                        ),
                        # Pseudo-time
                        # TODO
                        dcc.Tab(
                            label="Pseudo-time",
                            id="pseudo-time-tab",
                            children=[
                                html.H6("Palantir Pseudo-time and diff. potential"),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [
                                                HighResolutionGraph(
                                                    id="pseudo-time-graph"
                                                )
                                            ]
                                        ),
                                        dbc.Col(
                                            [HighResolutionGraph(id="entropy-graph")]
                                        ),
                                    ]
                                ),
                                html.H6("Branch probabilities"),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [
                                                HighResolutionGraph(
                                                    id="branch-probs-graph"
                                                )
                                            ]
                                        ),
                                    ]
                                ),
                                html.H6("Gene exression trends"),
                                dcc.RadioItems(
                                    id="trend-group-by",
                                    value="Branch",
                                    options=[
                                        {"label": "All branches", "value": "Branch"},
                                        {"label": "Th1", "value": "Th1"},
                                        {"label": "Tfh", "value": "Tfh"},
                                    ],
                                    labelStyle={
                                        "display": "inline-block",
                                        "margin-right": "20px",
                                        "margin-left": "5px",
                                    },
                                ),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [HighResolutionGraph(id="trends-graph")]
                                        ),
                                    ]
                                ),
                                # End tab
                            ],
                        ),
                        # Downloads tab
                        dcc.Tab(
                            label="Downloads",
                            id="download",
                            children=[
                                html.Br(),
                                html.Br(),
                                html.H5("GEO"),
                                dcc.Markdown(
                                    """
Raw count data and `scanpy Anndata` objects can be downloaded from GEO: add data download link
                    """
                                ),
                                html.Br(),
                                html.Br(),
                                html.Br(),
                                html.Br(),
                                # End tab
                            ],
                        ),
                        # Updates tab
                        dcc.Tab(
                            label="Updates",
                            id="update",
                            children=[
                                # End tab
                            ],
                        ),
                    ],
                ),
            ],
            className="eight-five columns",
        ),
    ]
)


# Functions
def matplotlib_to_plotly(colormap, pl_entries=255):
    cmap = matplotlib.colormaps.get_cmap(colormap)

    h = 1.0 / (pl_entries - 1)
    pl_colorscale = []

    for k in range(pl_entries):
        C = (np.array(cmap(k * h)[:3]) * 255).astype(np.uint8)
        pl_colorscale.append([k * h, "rgb" + str((C[0], C[1], C[2]))])

    return pl_colorscale


def load_data(dataset):
    store = pd.HDFStore(f"{data_dir}/{dataset}_metadata.h5", "r")
    meta_df = store["meta"]
    store.close()

    return meta_df


def construct_layout(data_df):
    layout = data_df.loc[:, ["fa_1", "fa_2"]]
    layout.columns = [0, 1]
    return layout


def plot_discrete(layout, values, ct_color=True):
    # Traces
    traces = []
    for i in np.sort(values.astype("str").unique()):
        marker = {"size": 5, "opacity": 0.75, "line": dict(width=0.3)}
        # ct = re.sub('[1-9].', '', i)
        # if ct_color and ct in an.ct_colors.index:
        #     marker['color'] = an.ct_colors[ct]

        traces.append(
            go.Scattergl(
                x=layout.loc[values == i, 0],
                y=layout.loc[values == i, 1],
                mode="markers",
                marker=marker,
                name=i,
                hoverinfo=None,
            )
        )

    return {
        "data": traces,
        "layout": go.Layout(
            autosize=False,
            width=500,
            height=500,
            xaxis={
                "title": "",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            yaxis={
                "title": "",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            legend={"x": 1, "y": 1},
        ),
    }


def plot_continuous(layout, values, gene, colormap="magma"):
    colorscale = matplotlib_to_plotly(colormap)
    # else:
    #     colors = cl.scales['11']['div'][colormap][::-1]
    #     colorscale = [[i, j] for i, j in zip(np.linspace(0, 1.0, len(colors)), colors)]
    values = values[layout.index]

    # Traces
    traces = [
        go.Scattergl(
            x=layout[0],
            y=layout[1],
            mode="markers",
            marker={
                "size": 5,
                "opacity": 0.75,
                "cmin": np.percentile(values, 0),
                "cmax": np.percentile(values, 100),
                "color": values.values,
                "colorscale": colorscale,
                "showscale": True,
                "line": dict(width=0.5),
            },
            hoverinfo=None,
        )
    ]

    return {
        "data": traces,
        "layout": go.Layout(
            autosize=False,
            width=500,
            height=500,
            xaxis={
                "title": "",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            yaxis={
                "title": "",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            title=gene,
        ),
    }


def plot_violin(
    values, metadata, condition_data, gene, secondary=None, xlabel="", ct_color=True
):
    traces = []

    # Reorder based on metadata
    condition_data = condition_data.sort_values()
    values = values[condition_data.index]

    # Row selection
    row_selection = condition_data
    if secondary is not None:
        row_selection = secondary[condition_data.index]

    for i in np.sort(row_selection.unique()):
        # plot differences for each condition in each cell type

        trace_dict = {
            "type": "violin",
            "x": metadata[row_selection == i],
            "y": values[row_selection == i],
            "color": condition_data,
            "name": i,
            "legendgroup": i,
            "scalegroup": i,
            "box": {"visible": True},
            "meanline": {"visible": True},
        }

        # trace_dict = {
        #     'type': 'violin',
        #     'x': metadata[row_selection == i],
        #     'y': values[row_selection == i],
        #     'name': i,
        #     "legendgroup": i,
        #     "scalegroup": i,
        #     'box': {'visible': True},
        #     'meanline': {'visible': True}
        # }

        # ct = re.sub('[1-9].', '', i)
        # cond = i
        # if ct_color and cond in an.ct_colors.index:
        #    trace_dict['line_color'] = {'color': an.ct_colors[ct]}

        traces.append(trace_dict)

    # Layout
    layout = go.Layout(
        autosize=True,
        yaxis={"title": "Normalized log gene expression", "zeroline": False},
        xaxis={"title": xlabel},
        hovermode="closest",
        title=gene,
    )
    # if secondary is not None:
    layout = layout.update({"violinmode": "group"})

    return {"data": traces, "layout": layout}


def load_trends(dataset, probs, genes):
    logging.info("loading trends")
    # Load trends dataframe
    trends = pd.Series()
    stds = pd.Series()
    branches = probs.columns.str.replace("^Prob_", "", regex=True)

    for br in branches:
        trends[br] = (
            feather.read_dataframe(
                f"{data_dir}/{dataset}_Prob_{br}_trends.feather", columns=genes
            )[genes].T
            + 3.2
        )
        stds[br] = feather.read_dataframe(
            f"{data_dir}/{dataset}_Prob_{br}_std.feather", columns=genes
        )[genes].T

        # column names
        bins = pd.read_pickle(f"{data_dir}/{dataset}_Prob_{br}_bins.p")
        trends[br].columns = bins
        stds[br].columns = bins

    return trends, stds


def get_dataframe(dataset, session_id):
    @cache.memoize()
    def query_and_serialize_data(dataset, session_id):
        # Load data
        df = load_data(dataset)

        return df.to_json()

    # Dataframe
    df = pd.read_json(query_and_serialize_data(dataset, session_id))
    layout = construct_layout(df)

    return layout, df


# plot gene trend across branches
def plot_trends_groupby_branch(dataset, probs, genes):
    logging.info("plotting trends group by branch")

    # Colors
    colors = (np.array(sns.color_palette(n_colors=probs.shape[1])) * 255).astype(
        np.uint8
    )

    trends, stds = load_trends(dataset, probs, genes)

    # Setup figure
    n_genes = len(genes)

    fig = subplots.make_subplots(rows=n_genes, cols=1, subplot_titles=genes)
    # Gene trace list
    showlegend = True

    # modify branch names
    branches = probs.columns.str.replace("^Prob_", "", regex=True)
    for i, gene in enumerate(genes):
        if i > 0:
            showlegend = False

        for j, br in enumerate(branches):
            # Mean trace
            trace = go.Scatter(
                x=trends[br].columns,
                y=trends[br].loc[gene, :],
                mode="lines",
                name=br,
                showlegend=showlegend,
                hoverinfo=None,
                line=dict(color=f"rgb({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]})"),
            )
            fig.append_trace(trace, i + 1, 1)

            # Std trace
            gene_std = stds[br].loc[gene, :]
            trend_stds = pd.Series()
            trend_stds["upper"] = list(trends[br].loc[gene, :] + gene_std)
            trend_stds["lower"] = list(trends[br].loc[gene, :] - gene_std)
            for k in trend_stds.index:
                trace = go.Scatter(
                    x=trends[br].columns,
                    y=trend_stds[k],
                    mode="lines",
                    showlegend=False,
                    hoverinfo=None,
                    fill="tonextx",
                    fillcolor=f"rgba({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]}, 0.2)",
                    line=dict(
                        color=f"rgba({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]}, 0)"
                    ),
                )
                fig.append_trace(trace, i + 1, 1)

    # Layout properties
    fig["layout"].update(height=400 * n_genes, width=800, autosize=False)
    for i in range(1, n_genes + 1):
        fig["layout"][f"xaxis{i}"].update(
            {
                "title": "Pseudo-time",
                # "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            }
        )
        fig["layout"][f"yaxis{i}"].update(
            {
                "title": "Gene expression",
                "zeroline": False,
                # "showgrid": False,
                "showticklabels": True,
            }
        )

    return fig


# plot gene trends colored by branches
def plot_trends_groupby_gene(dataset, probs, genes, branches=None):
    logging.info("plotting trends group by gene")

    # Colors
    colors = (np.array(sns.color_palette(n_colors=len(genes))) * 255).astype(np.uint8)

    trends, stds = load_trends(dataset, probs, genes)

    # Setup figure
    if branches is None:
        branches = probs.columns.str.replace("^Prob_", "", regex=True)
    n_branches = len(branches)
    fig = subplots.make_subplots(rows=n_branches, cols=1, subplot_titles=branches)
    # Gene trace list
    showlegend = True
    for i, br in enumerate(branches):
        if i > 0:
            showlegend = False

        for j, gene in enumerate(genes):
            # Mean trace
            trace = go.Scatter(
                x=trends[br].columns,
                y=minmax_scale(trends[br].loc[gene, :]),
                mode="lines",
                name=gene,
                showlegend=showlegend,
                hoverinfo=None,
                line=dict(color=f"rgb({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]})"),
            )
            fig.append_trace(trace, i + 1, 1)

            # Std trace
            gene_std = stds[br].loc[gene, :] / (
                np.max(trends[br].loc[gene, :]) - np.min(trends[br].loc[gene, :])
            )
            trend_stds = pd.Series()
            trend_stds["upper"] = minmax_scale(trends[br].loc[gene, :]) + gene_std
            trend_stds["lower"] = minmax_scale(trends[br].loc[gene, :]) - gene_std
            for k in trend_stds.index:
                trace = go.Scatter(
                    x=trends[br].columns,
                    y=trend_stds[k],
                    mode="lines",
                    showlegend=False,
                    hoverinfo=None,
                    fill="tonextx",
                    fillcolor=f"rgba({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]}, 0.2)",
                    line=dict(
                        color=f"rgba({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]}, 0)"
                    ),
                )
                fig.append_trace(trace, i + 1, 1)

    # Layout properties
    fig["layout"].update(height=400 * n_branches, width=800, autosize=False)
    for i in range(1, n_branches + 1):
        fig["layout"][f"xaxis{i}"].update(
            {
                "title": "Pseudo-time",
                # "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            }
        )
        fig["layout"][f"yaxis{i}"].update(
            {
                "title": "Gene expression",
                "zeroline": False,
                # "showgrid": False,
                "showticklabels": True,
            }
        )

    return fig


#### Callbacks ####


@app.callback(
    [
        # Output("dataset-dropdown", "disabled"),
        Output("gene-name", "disabled"),
        Output("multi-gene-name", "disabled"),
        Output("multi-gene-name", "value"),
        Output("dataset-dropdown", "options"),
        # Output("dataset-dropdown", "value"),
    ],
    [Input("tabs", "value")],
)
def tab_selection(selected_tab):

    if selected_tab in ["tab-1", "tab-2", "tab-4", "tab-5"]:

        # DAtaset options
        ds_options = [
            {"label": "log normalized expression", "value": "d1d2_normlog"},
            {"label": "MAGIC imputed expression", "value": "d1d2_imputed"},
        ]
        logging.info("returning options")
        # return False, False, True, "ISG15", ds_options  # , "d1d2_normlog"
        return False, True, "ISG15", ds_options  # , "d1d2_normlog"

    elif selected_tab in ["tab-3"]:  # , 'tab-2', 'tab-3']:
        # DAtaset options
        ds_options = [
            {"label": "Effector", "value": "d1d2_effector"},
            {"label": "TCR dependent", "value": "d1d2_TCRhilo"},
        ]
        logging.info("returning options tab 3")
        # return False, True, False, ["ISG15", "MX1"], ds_options  # , "d1d2_effector"
        return True, False, ["ISG15", "MX1"], ds_options  # , "d1d2_effector"


# dropdown options
@app.callback(
    [
        Output("gene-name", "options"),
        Output("multi-gene-name", "options"),
        Output("trend-group-by", "options"),
    ],
    [
        Input("dataset-dropdown", "value"),
    ],
)
def update_options(dataset):
    logging.info("update options")

    # Genes in dropdown
    genes = pd.read_csv(
        f"{data_dir}/{dataset}_genes.csv", header=None, index_col=None
    ).iloc[:, 0]

    gene_options = [{"label": i, "value": i} for i in np.sort(genes)]

    # dropdown options
    if dataset == "d1d2_effector":
        branch_options = [
            {"label": "All branches", "value": "Branch"},
            {"label": "Th1", "value": "Th1"},
            {"label": "Tfh", "value": "Tfh"},
        ]
    elif dataset == "d1d2_TCRhilo":
        branch_options = [
            {"label": "All branches", "value": "Branch"},
            {"label": "TCR independent", "value": "TCR_independent"},
            {"label": "TCR dependent", "value": "TCR_dependent"},
        ]
    else:
        branch_options = [
            {"label": "All branches", "value": "Branch"},
            {"label": "Th1", "value": "Th1"},
            {"label": "Tfh", "value": "Tfh"},
        ]

    return gene_options, gene_options, branch_options


# plotting on FDL
@app.callback(
    [
        Output("expression-graph", "figure"),
        Output("cluster-expression", "figure"),
        Output("meta-graph-1", "figure"),
        Output("meta-graph-2", "figure"),
    ],
    [
        Input("dataset-dropdown", "value"),
        Input("gene-name", "value"),
        Input("meta-field-1", "value"),
        Input("meta-field-2", "value"),
        Input("session-id", "children"),
    ],
)
def plot_gene_on_fdl(dataset, gene, meta_val_1, meta_val_2, session_id):
    logging.info("plotting on FDL")

    if dataset in ["d1d2_effector", "d1d2_TCRhilo"]:

        return go.Figure(), go.Figure(), go.Figure(), go.Figure()

    else:

        # Load layout
        layout, data_df = get_dataframe(dataset, session_id)

        # Gene data
        exprs = feather.read_dataframe(
            f"{data_dir}/{dataset}_data.feather", columns=[gene]
        )[gene]

        exprs.index = layout.index

        # Continuous plot
        scatter = plot_continuous(layout, exprs, gene)

        # Violin plot
        ct_violin = plot_violin(
            exprs,
            data_df["cluster"],
            data_df["condition"],
            gene,
            xlabel="Cluster",
            ct_color=False,
        )

        # metadata plot 1
        meta1_scatter = plot_discrete(layout, data_df[meta_val_1])

        # metadata plot 2
        meta2_scatter = plot_discrete(layout, data_df[meta_val_2])

        return scatter, ct_violin, meta1_scatter, meta2_scatter


# pseudotime FDL
@app.callback(
    [
        Output("pseudo-time-graph", "figure"),
        Output("entropy-graph", "figure"),
        Output("branch-probs-graph", "figure"),
    ],
    [Input("dataset-dropdown", "value"), Input("session-id", "children")],
)
def plot_pseudotime(dataset, session_id):
    if dataset not in ["d1d2_effector", "d1d2_TCRhilo"]:
        return go.Figure(), go.Figure(), go.Figure()

    logging.info("getting to pseudotime")

    # Load data
    layout, data_df = get_dataframe(dataset, session_id)

    # Trajectory and entropy
    pt_graph = plot_continuous(layout, data_df["pseudotime"], "Pseudo-time", "plasma")
    ent_graph = plot_continuous(layout, data_df["entropy"], "Diff. potential", "plasma")

    # Branch probabilities
    trace_list = dict()
    for i in data_df.columns[data_df.columns.str.contains("^Prob")]:
        trace_list[i] = plot_continuous(layout, data_df[i], i, "plasma")["data"][0]
        trace_list[i].update(showlegend=False)

    # Set up figure
    fig = subplots.make_subplots(
        rows=1, cols=len(trace_list), subplot_titles=list(trace_list.keys())
    )
    for i, k in enumerate(trace_list.keys()):
        fig.append_trace(trace_list[k], 1, i + 1)
        fig["layout"][f"xaxis{i+1}"].update(
            {"title": "", "showgrid": False, "zeroline": False, "showticklabels": False}
        )
        fig["layout"][f"yaxis{i+1}"].update(
            {"title": "", "showgrid": False, "zeroline": False, "showticklabels": False}
        )
    fig["layout"].update(height=400, width=400 * len(trace_list))

    return pt_graph, ent_graph, fig


# gene expression trends in pseudotime
@app.callback(
    Output("trends-graph", "figure"),
    [
        Input("dataset-dropdown", "value"),
        Input("trend-group-by", "value"),
        Input("multi-gene-name", "value"),
        Input("session-id", "children"),
    ],
)
def update_trends(dataset, group_by, genes, session_id):
    logging.info("updating trends")

    if dataset not in ["d1d2_effector", "d1d2_TCRhilo"]:
        return go.Figure()

    # Load data
    layout, data_df = get_dataframe(dataset, session_id)

    # Load gene expression
    if type(genes) is str:
        genes = [genes]

    # Plot
    probs = data_df.loc[:, data_df.columns[data_df.columns.str.contains("^Prob")]]
    probs.columns = probs.columns.str.replace("^Prob_", "", regex=True)

    if len(genes) > 0:
        if group_by == "Branch":
            graph = plot_trends_groupby_branch(dataset, probs, genes)
        else:
            branch = group_by
            graph = plot_trends_groupby_gene(dataset, probs, genes, [branch])
    else:
        graph = go.Figure()

    return graph


if __name__ == "__main__":
    logger.info("Starting the dash app")
    app.run_server(debug=True, host="0.0.0.0", port=8050)
