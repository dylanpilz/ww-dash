import dash
from dash import dcc, html, Input, Output, callback, dash_table
import pandas as pd
import plotly.express as px
import pickle
import json
import plotly.graph_objects as go
from scipy import signal
from datetime import date

dash.register_page(__name__, path='/')

dat = pickle.load(open("samples_deconv_dict.pkl",'rb'))
df_meta = pd.read_csv('wastewater_ncbi_ALL.csv',index_col=0)

layout = html.Div([ 
            html.Div([
                html.P(['This is a simple resource to query for specific mutations and lineages across SRA wastewater samples.',html.Br(),
                'All data is publicly available and analyses are performed using the Freyja package.']),
                ], style = {'textAlign': 'center'}),
            html.Br(),
            html.Div([
                html.Div(dcc.Graph(id='lineplot',config={'displayModeBar': False}),
                        style={
                                # "width": "68%",
                                # "height": "800px",
                                "textAlign": "center",
                                "display": "inline-block",
                                "padding-top": "5px",
                                "padding-left": "1px",
                                "overflow": "hidden",
                                'marginLeft': 'auto',
                                'marginRight': 'auto',
                        }
                )
            ], style = {'textAlign': 'center'}),
            html.Br(),
            html.Div([
                html.Div(dcc.Graph(id='stackplot',config={'displayModeBar': False}),
                        style={
                                # "width": "68%",
                                # "height": "800px",
                                "textAlign": "center",
                                "display": "inline-block",
                                "padding-top": "5px",
                                "padding-left": "1px",
                                "overflow": "hidden",
                                'marginLeft': 'auto',
                                'marginRight': 'auto',
                        }
                )
            ], style = {'textAlign': 'center'}),
            html.Br(),
            html.Div([html.Img(src=dash.get_asset_url('../assets/ncbi-logo.jpeg')),html.Img(src=dash.get_asset_url('../assets/search-logo.png'),style={'height':'25%', 'width':'25%'}),
                ], style={'textAlign': 'center'}),
        ])

def aggregate_viral_load(metadata):
    
    metadata['collection_date'] = pd.to_datetime(metadata['collection_date'], format='mixed')

    # Drop rows where ww_surv_target_1_conc is NaN or empty
    metadata = metadata.dropna(subset=['ww_surv_target_1_conc_unit'])
    metadata = metadata[metadata['ww_surv_target_1_conc_unit'].str.lower().str.contains('copies/l')]
    metadata = metadata[metadata['ww_surv_target_1_conc'] != '']
    metadata = metadata[metadata['ww_surv_target_1_conc'] != 'not provided']
    metadata = metadata[metadata['ww_surv_target_1_conc'] != 'missing']
    metadata = metadata[metadata['ww_surv_target_1_protocol'] != 'Multiplex SARS-Flu-RSV-Mpox dPCR'] 

    metadata['ww_population'] = metadata['ww_population'].astype(str)
    metadata = metadata[~metadata['ww_population'].str.contains('<')]
    metadata = metadata[~metadata['ww_population'].str.contains('>')]

    metadata['ww_surv_target_1_conc_unit'] = metadata['ww_surv_target_1_conc_unit'].str.lower()
    metadata['ww_surv_target_1_conc'] = metadata['ww_surv_target_1_conc'].astype(float)

    #metadata.loc[metadata['ww_surv_target_1_conc_unit'].str.contains('copies/ul'), 'ww_surv_target_1_conc'] *= 1000000

    # print(metadata['ww_surv_target_1_conc_unit'].value_counts())
    # print(metadata['collection_date'].value_counts())
    metadata['ww_population'] = metadata['ww_population'].astype(float)
    metadata['viral_load (copies/L)'] = metadata['ww_surv_target_1_conc']

    return metadata

@callback(
    Output('lineplot', 'figure'),
    Input('lineplot', 'figure'))
def update_viral_load(input):
    import plotly.graph_objs as go

    agg_df = aggregate_viral_load(df_meta)

    print(agg_df['viral_load (copies/L)'])
    agg_df = agg_df.groupby('collection_date').agg({'viral_load (copies/L)': 'mean'}).reset_index()

    fig = go.Figure(layout=dict(template='plotly'))
    # signal.savgol_filter(agg_df['viral_load (copies/L)'], 30, 2)
    #fig = go.Figure(data=go.Scatter(x=agg_df['collection_date'], y=agg_df['Copies/L-person'], mode='lines+markers'))
    fig = px.line(agg_df, x='collection_date', y=signal.savgol_filter(agg_df['viral_load (copies/L)'], 30, 2), title='Viral Load in USA Wastewater Samples', markers=True, line_shape='spline')
    fig.update_layout(
        transition_duration=100,
        title_text='Viral Load in USA Wastewater Samples',
        title_x=0.5,
        showlegend=False,
        xaxis_title='Collection Date',
        yaxis_title='Viral Load (copies/L)'
    )

    return fig

@callback(
    Output('stackplot', 'figure'),
    Input('stackplot', 'figure'))
def update_lineage_prevalences(input):
    df_demix = pd.read_csv('aggregate_demix.tsv',sep='\t',index_col=0)

    df_demix.index = df_demix.index.str.split('.').str[0]

    df_demix = df_demix[df_demix['coverage'] > 60]

    df_metadata = aggregate_viral_load(df_meta)

    df_demix = df_demix[df_demix.index.isin(df_metadata.index)]


    # Add metadata columns to df_demix
    df_demix['collection_date'] = df_metadata.loc[df_demix.index, 'collection_date'].values

    df_demix['ww_population'] = df_metadata.loc[df_demix.index, 'ww_population'].values
    df_demix['ww_surv_target_1_conc'] = df_metadata.loc[df_demix.index, 'ww_surv_target_1_conc'].values
    df_demix['viral_load (copies/L)'] = df_metadata.loc[df_demix.index, 'viral_load (copies/L)'].values

    df_demix['summarized'] = df_demix['summarized'].apply(lambda x: list(eval(x)))

    # Make a histogram of the number of samples per collection date
    fig = px.histogram(df_demix, x='collection_date', title='Number of Samples per Collection Date')
    fig.write_image('histogram.png')

    # Multiply summarized lineage abundances by viral load and add for each collection date
    for row in df_demix.iterrows():
        scaled_summarized = []
        for tup in df_demix.loc[row[0], 'summarized']:
            scaled_summarized.append((tup[0], tup[1] * row[1]['viral_load (copies/L)']))
        df_demix.loc[row[0], 'scaled_summarized'] = str(scaled_summarized)
  
    
    # Combine scaled abundances for each collection date
    df_demix['scaled_summarized'] = df_demix['scaled_summarized'].apply(lambda x: list(eval(x)))
    df_demix['scaled_summarized'] = df_demix['scaled_summarized'].apply(lambda x: dict(x))


    def merge_dicts(x):
        merged_dict = {}
        for d in x:
            for k in d.keys():
                # Take the average if the key is already in merged_dict
                if k in merged_dict.keys():
                    merged_dict[k] = (merged_dict[k] + d[k]) / 2

                else:
                    merged_dict[k] = d[k]
        # Set the values such that they sum to 1
        merged_dict = {k: v / sum(merged_dict.values()) for k, v in merged_dict.items()}
        return merged_dict
    
    # Merge scaled abundances for each collection date
    plot_df = df_demix.groupby('collection_date').agg({'scaled_summarized': merge_dicts}).reset_index()
    
    for row in plot_df.iterrows():
        for k, v in row[1]['scaled_summarized'].items():
            plot_df.loc[row[0], k] = v
    plot_df = plot_df.drop('scaled_summarized', axis=1)
    plot_df = plot_df.fillna(0)
    plot_df.index = plot_df['collection_date']
    plot_df = plot_df.drop('collection_date', axis=1)
    print(plot_df)
    plot_df = plot_df.rolling(7).mean()
    
    
    fig = go.Figure(layout=dict(template='plotly'))
    # Create a plotly stackplot to show the prevalence of each lineage over time
    fig = px.area(plot_df, title='Prevalence of SARS-CoV-2 Lineages in USA Wastewater Samples')

    return fig
