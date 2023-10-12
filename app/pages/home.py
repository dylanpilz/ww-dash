import dash
from dash import dcc, html, Input, Output, callback, dash_table
import pandas as pd
import plotly.express as px
import pickle
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
            html.Div([html.Img(src=dash.get_asset_url('../assets/ncbi-logo.jpeg')),html.Img(src=dash.get_asset_url('../assets/search-logo.png'),style={'height':'25%', 'width':'25%'}),
                ], style={'textAlign': 'center'}),
        ])

def aggregate_viral_load(metadata):
    metadata['collection_date'] = pd.to_datetime(metadata['collection_date'], format='mixed')
    # Drop rows where ww_surv_target_1_conc is NaN
    metadata = metadata.dropna(subset=['ww_surv_target_1_conc_unit'])
    metadata = metadata[metadata['ww_surv_target_1_conc'] != 'not provided']
    metadata = metadata[metadata['ww_surv_target_1_conc'] != 'missing']
    metadata = metadata[metadata['ww_surv_target_1_protocol'] != 'Multiplex SARS-Flu-RSV-Mpox dPCR'] 

    metadata['ww_surv_target_1_conc_unit'] = metadata['ww_surv_target_1_conc_unit'].str.lower()
    metadata['ww_surv_target_1_conc'] = metadata['ww_surv_target_1_conc'].astype(float)
    metadata.loc[metadata['ww_surv_target_1_conc_unit'].str.contains('copies/l'), 'ww_surv_target_1_conc'] /= 1000
    #metadata.loc[metadata['ww_surv_target_1_conc_unit'].str.contains('copies/ul'), 'ww_surv_target_1_conc'] *= 1000000

    print(metadata['ww_surv_target_1_conc_unit'].value_counts())
    # Create a new dataframe with the aggregated viral load for each date
    metadata['collection_date'].value_counts().to_csv('date_counts.csv')

    metadata_agg = metadata.groupby('collection_date').agg({'ww_surv_target_1_conc': 'mean'}).reset_index()

    # Print the number of samples for each date
    return metadata_agg


@callback(
    Output('lineplot', 'figure'),
    Input('lineplot', 'figure'))
def update_figure(input):
    agg_df = aggregate_viral_load(df_meta)

    fig = px.line(agg_df, x='collection_date', y='ww_surv_target_1_conc', title='Viral Load in Wastewater Samples')
    fig.update_layout(
        transition_duration=100,
        title_text='Viral Load in USA Wastewater Samples',
        title_x=0.5,
        showlegend=False,
        xaxis_title='Collection Date',
        yaxis_title='Viral Load (copies/mL)'
    )

    return fig
