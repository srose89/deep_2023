# -*- coding: utf-8 -*-
import dash
from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Input, Output
import plotly.graph_objs as go
from plotly import tools

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
                        id='multi-gene-name',
                        value='Isg15',
                        disabled=True,
                        multi=True
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
                        # dcc.Tab(
                        #     label="Pseudo-time",
                        #     id="pseudo-time-tab",
                        #     children=[
                        #         html.H6("Palantir Pseudo-time and diff. potential"),
                        #         dbc.Row(
                        #             [
                        #                 dbc.Col(
                        #                     [
                        #                         HighResolutionGraph(
                        #                             id="pseudo-time-graph"
                        #                         )
                        #                     ]
                        #                 ),
                        #                 dbc.Col(
                        #                     [HighResolutionGraph(id="entropy-graph")]
                        #                 ),
                        #             ]
                        #         ),
                        #         html.H6("Branch probabilities"),
                        #         dbc.Row(
                        #             [
                        #                 dbc.Col(
                        #                     [
                        #                         HighResolutionGraph(
                        #                             id="branch-probs-graph"
                        #                         )
                        #                     ]
                        #                 ),
                        #             ]
                        #         ),
                        #         html.H6("Gene exression trends"),
                        #         dcc.RadioItems(
                        #             id="trend-group-by",
                        #             value="Branch",
                        #             options=[
                        #                 {"label": "All branches", "value": "Branch"},
                        #                 {"label": "Epi", "value": "EPI"},
                        #                 {"label": "PrE/VE", "value": "PrE"},
                        #             ],
                        #             labelStyle={
                        #                 "display": "inline-block",
                        #                 "margin-right": "20px",
                        #                 "margin-left": "5px",
                        #             },
                        #         ),
                        #         dbc.Row(
                        #             [
                        #                 dbc.Col(
                        #                     [HighResolutionGraph(id="trends-graph")]
                        #                 ),
                        #             ]
                        #         ),
                        #         # End tab
                        #     ],
                        # ),
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
    cmap = matplotlib.cm.get_cmap(colormap)

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


# def plot_violin(values, metadata, gene, secondary=None, xlabel=''):
#     traces = []

#     # Reorder based on metadata
#     metadata = metadata.sort_values()
#     values = values[metadata.index]

#     # Row selection
#     row_selection = metadata
#     if secondary is not None:
#         row_selection = secondary[metadata.index]

#     for i in np.sort(row_selection.unique()):
#         traces.append({
#             'type': 'violin',
#             'x': metadata[row_selection == i],
#             'y': values[row_selection == i],
#             'name': i,
#             "legendgroup": i,
#             "scalegroup": i,
#             'box': {'visible': True},
#             'meanline': {'visible': True}
#         })

#     # Layout
#     layout = go.Layout(
#         autosize=True,
#         yaxis={'title': 'Imputed log gene expression', 'zeroline': False},
#         xaxis={'title': xlabel},
#         hovermode='closest',
#         title=gene
#     )
#     if secondary is not None:
#         layout = layout.update({'violinmode': 'group'})

#     return {
#         'data': traces,
#         'layout': layout
#     }



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


# def plot_heatmap(trends):
#     print(trends)
#     colors = cl.scales["11"]["div"]["Spectral"][::-1]
#     colorscale = [[i, j] for i, j in zip(np.linspace(0, 1.0, len(colors)), colors)]

#     #  Append each genes
#     z = []
#     for i in trends.index:
#         z.append(list(scale(trends.loc[i, :])))
#     print(z)

#     data = [
#         go.Heatmap(
#             z=z,
#             x=np.repeat("", trends.shape[1]),
#             y=list(trends.index),
#             colorscale=colorscale,
#             zmin=-0.25,
#             zmax=2.75,
#         )
#     ]

#     layout = go.Layout(
#         title="Gene expression along E8.75 guttube pseudo-space",
#         width=700,
#         height=900,
#         xaxis=dict(ticks=""),
#         yaxis=dict(ticks=""),
#     )

#     return {"data": data, "layout": layout}


def load_trends(dataset, probs, genes):
    # Load trends dataframe
    trends = pd.Series()
    stds = pd.Series()
    for br in probs.columns:
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


# Callbacks

# gene options
@app.callback(
    [Output("gene-name", "options"),
    Output('multi-gene-name', 'options')]
    #  Output('timepoint-element', 'children')
    ,
    [Input("dataset-dropdown", "value"), Input("tabs", "value")],
)
def update_options(dataset, selected_tab):
    # Genes in dropdown
    genes = pd.read_csv(
        f"{data_dir}/{dataset}_genes.csv", header=None, index_col=None
    ).iloc[:, 0]
    
    # Radio items
    # if dataset == 'E85_guttube':
    # meta_options = [
    #     {'label': 'Cell type   ', 'value': 'CellType'},
    #     {'label': 'Lineage  ', 'value': 'Timepoint'},
    #     {'label': 'Clusters   ', 'value': 'clusters'}
    # ]
    # timepoint_element = 'Lineage expression'
    # meta_options = [
    #         {'label': 'Cell type   ', 'value': 'Celltype'},
    #         {'label': 'Condition  ', 'value': 'condition'},
    # ]
    # else:
    #     meta_options = [
    #         {'label': 'Cell type   ', 'value': 'CellType'},
    #         {'label': 'Time point  ', 'value': 'Timepoint'},
    #         {'label': 'Clusters   ', 'value': 'clusters'}
    #     ]
    #     timepoint_element = 'Timepoint expression'

    gene_options = [{"label": i, "value": i} for i in np.sort(genes)]
    # return gene_options #, meta_options#, timepoint_element
    return gene_options, gene_options

# dropdown options
@app.callback(
    [
        Output("dataset-dropdown", "options"),
       # Output("dataset-dropdown", "value"),
        Output("gene-name", "disabled"),
        Output('multi-gene-name', 'disabled'),
        Output('multi-gene-name', 'value'),
    ],
    [Input("tabs", "value")],
)
def update_dropdown(selected_tab):
    if selected_tab in ["tab-1", "tab-2"]:  # , 'tab-2', 'tab-3']:
        # DAtaset options
        ds_options = [
            {"label": "log normalized expression", "value": "d1d2_normlog"},
            {"label": "MAGIC imputed expression", "value": "d1d2_imputed"},
        ]

        #return  ds_options, "d1d2_normlog", False
        return  ds_options, False, False, 'Isg15'

    elif selected_tab in ["tab-3"]:  # , 'tab-2', 'tab-3']:
        # DAtaset options
        ds_options = [
            {"label": "Effector", "value": "d1d2_effector"},
            {"label": "TCR dependent", "value": "d1d2_TCRhilo"},
        ]
        #return  ds_options, "d1d2_effector", False
        return  ds_options,  False, True, 'Isg15'

    # else:  # , 'tab-2', 'tab-3']:
    #     # DAtaset options
    #     ds_options = [
    #         {"label": "log normalized expression", "value": "d1d2_normlog"},
    #         {"label": "MAGIC imputed expression", "value": "d1d2_imputed"},
    #     ]
    #     return  ds_options, "d1d2_normlog", False


# @app.callback([Output('dataset-dropdown', 'disabled'),
#                Output('gene-name', 'disabled'),
#                Output('multi-gene-name', 'disabled'),
#                Output('multi-gene-name', 'value'),
#                Output('dataset-dropdown', 'options'),
#                Output('dataset-dropdown', 'value')],
#               [Input('tabs', 'value'),
#                Input('session-id', 'children')])
# def tab_selection(selected_tab, session_id):

#     print("getting to tab selection")
#     print(session_id)
#     if selected_tab in ['tab-1', 'tab-2']: #, 'tab-2', 'tab-3']:
#         # DAtaset options
#         ds_options = [
#             {'label': 'log normalized expression', 'value': 'd1d2_normlog'},
#             {'label': 'MAGIC imputed expression', 'value': 'd1d2_imputed'}

#         ]
#         return False, False, True, 'Isg15', ds_options, 'd1d2_normlog'

#     elif selected_tab in ['tab-3']: #, 'tab-2', 'tab-3']:
#         # DAtaset options
#         ds_options = [
#             {'label': 'Effector', 'value': 'd1d2_effector'},
#             {'label': 'TCR High/Low', 'value': 'd1d2_TCRhilo'}
#         ]
#         return False, False, True, 'Isg15', ds_options, 'd1d2_effector'

#     else: #, 'tab-2', 'tab-3']:
#         # DAtaset options
#         ds_options = [
#             {'label': 'log normalized expression', 'value': 'd1d2_normlog'},
#             {'label': 'MAGIC imputed expression', 'value': 'd1d2_imputed'}

#         ]
#         return False, False, True, 'Isg15', ds_options, 'd1d2_normlog'
# elif selected_tab in ['tab-4']:
#     ds_options = [
#         {'label': 'E3.5-E4.5 (EPI, PrE)', 'value': 'E35_E45_sub'},
#         {'label': 'E3.5-E5.5 (EPI, VE)', 'value': 'E35_E55_no_bridge'}
#     ]
#     return False, True, False, 'Fgf4', ds_options, 'E35_E45_sub'
# else:
#     # DAtaset options
#     ds_options = [
#         {'label': 'KP endothelial', 'value': 'KP_endo'},
#     ]
#     return False, True, False, 'Csf3', ds_options, 'KP_endo'





# @app.callback(
#     Output("meta-graph-1", "figure"),
#     [
#         Input("dataset-dropdown", "value"),
#         Input("meta-field", "value"),
#         Input("session-id", "children"),
#     ],
# )
# def plot_celltype(dataset, field, session_id):
#     # Load data
#     layout, data_df = get_dataframe(dataset, session_id)
#     ct_color = True
#     graph = plot_discrete(layout, data_df[field], ct_color)

#     return graph


# def plot_gene_on_fdl(dataset, gene, session_id):

#     # Load layout
#     layout, data_df = get_dataframe(dataset, session_id)

#     # Gene data
#     exprs = feather.read_dataframe(f'{data_dir}/{dataset}_data.feather', columns=[gene])[gene]
#     exprs = exprs + 3.3222656
#     exprs.index = layout.index

#     # Continuous plot
#     scatter = plot_continuous(layout, exprs, gene)

#     # Violin plot
#     tp_xlabel = 'Time point'
#     if dataset == 'E85_guttube':
#         tp_xlabel = 'Lineage'
#     ct_color = True
#     if dataset == 'E35_E85_VE':
#         ct_color = False
#     ct_violin = plot_violin(exprs, data_df['CellType'], gene, data_df['Timepoint'],
#                             xlabel='Cell type', ct_color=ct_color)
#     tp_violin = plot_violin(exprs, data_df['Timepoint'], gene, data_df['CellType'],
#                             xlabel=tp_xlabel, ct_color=ct_color)
#     cl_violin = plot_violin(exprs, data_df['clusters'], gene,
#                             xlabel='Cluster', ct_color=ct_color)

#     return scatter, ct_violin, tp_violin, cl_violin


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
    # Load layout
    layout, data_df = get_dataframe(dataset, session_id)

    #print(dataset)
    #print(data_df.head())
    # Gene data
    exprs = feather.read_dataframe(
        f"{data_dir}/{dataset}_data.feather", columns=[gene]
    )[gene]

    # exprs = exprs + 3.3222656
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
    # ct_color=ct_color)

    # metadata plot 1
    meta1_scatter = plot_discrete(layout, data_df[meta_val_1])

    # metadata plot 2
    meta2_scatter = plot_discrete(layout, data_df[meta_val_2])

    return scatter, ct_violin, meta1_scatter, meta2_scatter


# @app.callback(
#     [Output('ps-meta-graph', 'figure'),
#      Output('ps-graph', 'figure')],
#     [Input('ps-meta-field', 'value'),
#      Input('session-id', 'children')])
# def plot_pseudospace(field, session_id):
#     # Load data
#     layout, data_df = get_dataframe('E85_guttube', session_id)
#     meta_graph = plot_discrete(layout, data_df[field])

#     # Pseudo-space
#     ps_graph = plot_continuous(layout, data_df['pseudospace'], 'Pseudo-space')

#     return meta_graph, ps_graph


# @app.callback(
#     Output('ps-heatmap', 'figure'),
#     [Input('dataset-dropdown', 'value'),
#      Input('multi-gene-name', 'value'),
#      Input('session-id', 'children')])
# def plot_pseudospace_heatmap(dataset, genes, session_id):
#     if dataset != 'E85_guttube':
#         return go.Figure()

#     # Load data
#     layout, data_df = get_dataframe('E85_guttube', session_id)

#     # Load gene expression
#     if type(genes) is str:
#         genes = [genes]
#     exprs = feather.read_dataframe(f'{data_dir}/E85_guttube_data.feather', columns=genes)[genes[::-1]]
#     exprs.index = data_df.index

#     # Plot heatmap
#     # Load trends dataframe
#     trends = feather.read_dataframe(f'{data_dir}/{dataset}_pseudospace_trends.feather', columns=genes)[genes].T

#     # column names
#     bins = pd.read_pickle(f'{data_dir}/{dataset}_pseudospace_bins.p')
#     trends.columns = bins

#     graph = plot_heatmap(trends)

#     return graph


# Callbacks
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
    fig = tools.make_subplots(
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
    # def update_trends(dataset, group_by, genes, session_id):

    if dataset not in ["d1d2_effector", "d1d2_TCRhilo"]:
        return go.Figure(), go.Figure()

    # Load data
    layout, data_df = get_dataframe(dataset, session_id)

    # Load gene expression
    if type(genes) is str:
        genes = [genes]

    # Plot
    probs = data_df.loc[:, data_df.columns[data_df.columns.str.contains("^Prob")]]
    probs.columns = probs.columns.str.replace("^Prob_", "")
    print(probs.head())
    if len(genes) > 0:
        if group_by == "Branch":
            graph = plot_trends_groupby_branch(dataset, probs, genes)
        else:
            branch = group_by
            # if branch == 'PrE' and dataset == 'E35_E55_no_bridge':
            #     branch = 'VE'
            graph = plot_trends_groupby_gene(dataset, probs, genes, [branch])
    else:
        graph = go.Figure()
    
    return graph


def plot_trends_groupby_branch(dataset, probs, genes):
    # Trends grouped by gene

    # Colors
    colors = (np.array(sns.color_palette(n_colors=probs.shape[1])) * 255).astype(
        np.uint8
    )

    trends, stds = load_trends(dataset, probs, genes)

    # Setup figure
    n_genes = len(genes)
    fig = tools.make_subplots(rows=n_genes, cols=1, subplot_titles=genes)
    # Gene trace list
    showlegend = True
    for i, gene in enumerate(genes):
        if i > 0:
            showlegend = False
        for j, br in enumerate(probs.columns):
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
    # fig["layout"].update(height=400 * n_genes, width=800)
    # for i in range(1, n_genes + 1):
    #     fig["layout"][f"xaxis{i}"].update({"title": "Pseudo-time"})
    #     fig["layout"][f"yaxis{i}"].update(
    #         {"title": "Gene expression", "zeroline": False}
    #     )
    return {
        "data": fig,
        "layout": go.Layout(
            autosize=False,
            width=800,
            height=400,
            xaxis={
                "title": "Pseudo-time",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            yaxis={
                "title": "Gene expression",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            legend={"x": 1, "y": 1},
        )
        
    }
    #return fig


def plot_trends_groupby_gene(dataset, probs, genes, branches=None):
    # Trends grouped by gene

    # Colors
    colors = (np.array(sns.color_palette(n_colors=len(genes))) * 255).astype(np.uint8)

    trends, stds = load_trends(dataset, probs, genes)

    # Setup figure
    if branches is None:
        branches = probs.columns
    n_branches = len(branches)
    fig = tools.make_subplots(rows=n_branches, cols=1, subplot_titles=branches)
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
                mode='lines',
                name=gene,
                showlegend=showlegend,
                hoverinfo=None,
                line=dict(color=f'rgb({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]})'),
            )
            fig.append_trace(trace, i + 1, 1)

            # Std trace
            gene_std = stds[br].loc[gene, :] / (np.max(trends[br].loc[gene, :]) - np.min(trends[br].loc[gene, :]))
            trend_stds = pd.Series()
            trend_stds['upper'] = minmax_scale(trends[br].loc[gene, :]) + gene_std
            trend_stds['lower'] = minmax_scale(trends[br].loc[gene, :]) - gene_std
            for k in trend_stds.index:
                trace = go.Scatter(
                    x=trends[br].columns,
                    y=trend_stds[k],
                    mode='lines',
                    showlegend=False,
                    hoverinfo=None,
                    fill='tonextx',
                    fillcolor=f'rgba({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]}, 0.2)',
                    line=dict(color=f'rgba({colors[j, 0]}, {colors[j, 1]}, {colors[j, 2]}, 0)'),
                )
                fig.append_trace(trace, i + 1, 1)

    # Layout properties
    # fig['layout'].update(height=400 * n_branches, width=800)
    # for i in range(1, n_branches + 1):
    #     fig['layout'][f'xaxis{i}'].update({'title': 'Pseudo-time'})
    #     fig['layout'][f'yaxis{i}'].update({'title': 'Normalized expression', 'zeroline': False})

    return {
        "data": fig,
        "layout": go.Layout(
            autosize=False,
            width=800,
            height=400,
            xaxis={
                "title": "Pseudo-time",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            yaxis={
                "title": "Gene expression",
                "showgrid": False,
                "zeroline": False,
                "showticklabels": False,
            },
            legend={"x": 1, "y": 1},
        )
        
    }

    #return fig


# load gene weights for the factor
# write a function that takes in the dataset, factor name, and session id
# the function will then load the gene weights for the factor
# then the function will return a dataframe with the gene weights ranked by their weights for that factor
# @app.callback(Output('gene-weight-df', 'children'),
#                 [Input('dataset-dropdown', 'value'),
#                  Input('factor-name', 'value'),
#                 Input('session-id', 'children')])
# def load_gene_weights(dataset, factor_name, session_id):
#     # Load gene weights
#     data_df = get_gw_dataframe(dataset, session_id)

#     # extract column of data_df that corresponds to the factor
#     factor_df = data_df[factor_name]

#     # return the factor_df
#     return factor_df

# function that will load a figure from a png file
# and return it as a html.Img object
# based on the dataset selected
# @app.callback(
#     Output('DA-plot', 'src'),
#     [Input('dataset-dropdown', 'value')])
# def DAimage(dataset):
#     if dataset == 'KP_endo':
#         return '/assets/KP_endo_DA_plot.png'
#     elif dataset == 'KP_fib':
#         return '/assets/KP_fib_DA_plot.png'
#     elif dataset == 'KP_myl':
#         return '/assets/KP_myl_DA_plot.png'

# function that will load a figure from a png file
# and return it as a html.Img object
# based on the dataset selected
# this is for the visium tissue sections and tumor states
# @app.callback(
#     Output('tumor-state-img', 'src'),
#     [Input('dataset-dropdown', 'value')])
# def tumor_state_image(dataset):
#     if dataset == 'KP_ctrl1_A1':
#         return '/assets/KP_ctrl1_A1_plot.png'
#     elif dataset == 'KP_ctrl2_C1':
#         return '/assets/KP_ctrl2_C1_plot.png'
#     elif dataset == 'KP_DT_A1':
#         return '/assets/KP_DT_A1_plot.png'
#     elif dataset == 'KP_DT_C1':
#         return '/assets/KP_DT_C1_plot.png'


if __name__ == "__main__":
    #application.run(host = '0.0.0.0', debug=True, port=7777)
    app.run_server(debug=False, host='0.0.0.0', port=8050)
    #app.run_server(debug=False, port=7777)
