import dash
import pandas as pd
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html, dash_table, State, callback_context
import plotly.express as px
import io
import datetime
import plotly.graph_objs as go
from vLab.IntegratedBioprocess.PlantwiseSimulator import PlantwiseSimulator
import time
import numpy as np
import matplotlib.pyplot as plt
import base64
from tempfile import TemporaryFile


app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

IMG_STYLE = {
    "position": "absolute",
    "top": 100,
    "left": 700,
}

sidebar = html.Div(
    [
        html.H2("NIIMBL", className="display-4"),
        html.Hr(),
        dbc.Nav(
            [
                dbc.NavLink("Home", href="/", active="exact"),
                dbc.NavLink("Select Files", href="/page-1", active="exact"),
                dbc.NavLink("Plantwise Simulator", href="/page-2", active="exact"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)

content = html.Div(id="page-content", style=CONTENT_STYLE)

app.layout = html.Div([dcc.Location(id="url"), sidebar, content])

controls = dbc.Card(
    [
            html.Div([
                "X0: ",
                dcc.Slider(0, 5, 0.1, value=0.1, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id="input-X0")
                # dcc.Input(id="input-X0", type="number", min=0, max=20, value=0.1)
            ]),
            html.Br(),

            html.Div([
                "Sg0: ",
                dcc.Slider(30, 50, 1, value=40, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-Sg0')
                # dcc.Input(id='input-Sg0', type="number", min=0, max=40, value=40)
            ]),
            html.Br(),

            html.Div([
                "Sm0: ",
                dcc.Slider(0, 5, 0.1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-Sm0')
                # dcc.Input(id='input-Sm0', type="number", min=0, max=20, value=0)
            ]),
            html.Br(),

            html.Div([
                "P10: ",
                dcc.Slider(0, 5, 0.1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-P10')
                # dcc.Input(id='input-P10', type="number", min=0, max=20, value=0)
            ]),
            html.Br(),

            html.Div([
                "P20: ",
                dcc.Slider(0, 5, 0.1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-P20')
                # dcc.Input(id='input-P20', type="number", min=0, max=20, value=0)
            ]),
            html.Br(),

            html.Div([
                "P30: ",
                dcc.Slider(0, 5, 0.1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-P30')
                # dcc.Input(id='input-P30', type="number", min=0, max=20, value=0)
            ]),
            html.Br(),

            html.Div([
                "VB0: ",
                dcc.Slider(0, 5, 0.1, value=0.5, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-VB0')
                # dcc.Input(id='input-VB0', type="number", min=0, max=20, value=0.5)
            ]),
            html.Br(),

            html.Div([
                "VH0: ",
                dcc.Input(id='input-VH0', type="text", value=1e-8)
            ]),
            html.Br(),
    ],
    body=True,
)

progress = dbc.Progress(
    [
        dbc.Progress(value=20, label="Outgrowth", color="success", bar=True),
        dbc.Progress(value=10, label="Trans", color="warning", bar=True),
        dbc.Progress(value=70, label="Production", color="danger", bar=True),
    ], style={'position': 'absolute', 'top': 650, 'left': 700, 'width': 700, 'height': '30px'}
)

img = html.Div([
    dbc.Button('Bioreactor', color="primary", id='bio', n_clicks=0, className="me-1"),
    dbc.Button('Hold Tank', color="primary", id='hold', n_clicks=0, className="me-1"),
    dbc.Button('Capture(Bind/Elute)', color="primary", id='cap', n_clicks=0, className="me-1"),
    dbc.Button('Polish(Flow Through)', color="primary", id='po', n_clicks=0, className="me-1"),
    html.Div(id='output')
], style=IMG_STYLE)

home = html.Div(
    [
        html.H4("Interactive Biopharmaceutical Simulators", style={'textAlign':'center'}),
        html.Hr(),
        html.H5('Our Goal:', style={'color': 'blue'}),
        html.P('Provide the ability to investigate different biopharmaceutical process outcomes by manipulating environment parameters and feed rates without needing the physical setup.'),
        html.Br(),
        html.H5('How?', style={'color': 'blue'}),
        html.P('Our simulators act as a digital twin of real life biopharma processes!'),
        html.P('An ideal digital twin acts perfectly like its real-life counterpart. For example, we can simulate a fermentation or antibody production cycle and extract data on the process outputs, all from behind a desktop.'),
        html.P('This allows scientists, technicians, or even curious students to see how manipulating feed rates and process parameters affect the outputs of biopharma processes.'),
        html.P('Access the simulators using the buttons on the left or continue to learn more about biopharmaceuticals!')
    ],
)

page2 = html.Div(
    [
        html.H4("Plantwise Simulator", style={'textAlign':'center'}),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(controls, md=4),
                html.Div([dbc.Col(img, md=8), dbc.Col(progress, md=8)])
            ],
            align="center",
        ),

    ],
)

page1 = html.Div(
    [
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=True
        ),
        html.Div(id='output-data-upload'),
    ],
)

def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        dash_table.DataTable(
            df.to_dict('records'),
            [{'name': i, 'id': i} for i in df.columns]
        ),

        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        html.Div('Raw Content'),
        html.Pre(contents[0:200] + '...', style={
            'whiteSpace': 'pre-wrap',
            'wordBreak': 'break-all'
        })
    ])

@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname == "/":
        return home
    elif pathname == "/page-1":
        return page1
    elif pathname == "/page-2":
        return page2

    return dbc.Jumbotron(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ]
    )

@app.callback(Output('output-data-upload', 'children'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children

@app.callback(
    Output('output', 'children'),
    Input('bio', 'n_clicks'),
    Input('input-X0', 'value'),
    Input('input-Sg0', 'value'),
    Input('input-Sm0', 'value'),
    Input('input-P10', 'value'),
    Input('input-P20', 'value'),
    Input('input-P30', 'value'),
    Input('input-VB0', 'value'),
    Input('input-VH0', 'value')
)
def plot_input(n_clicks, v1, v2, v3, v4, v5, v6, v7, v8):
    changed_id = [p['prop_id'] for p in callback_context.triggered][0]
    if 'bio' in changed_id:
        # x0 = [v1, v2, v3, v4, v5, v6, v7, v4, v5, v6, v8]
        # xC0 = [0] * (10 * 30 + 3)
        # x0 = x0 + xC0
        # outfile = TemporaryFile()
        # np.save(outfile, x0)
        # _ = outfile.seek(0)
        # # start_time = time.time()
        # solver = PlantwiseSimulator()
        # x0 = np.load(outfile, encoding='bytes', allow_pickle=True)
        # sol = solver.solve(x0)
        # elapse_time = time.time() - start_time
        # t = np.array(sol.t)
        # x = np.array(sol.x)
        # plt.plot(t, x[:, 3:6])
        # plt.axvline(solver._process_time[1], ls='--', c='k')
        # plt.axvline(solver._process_time[2], ls='--', c='k')
        # plt.title('Bioreactor', fontsize=14)
        # plt.ylabel('Concentration (mg/mL)', fontsize=14)
        # plt.xlabel('Time (h)', fontsize=14)
        # plt.legend(['Product', 'Impurity 1', 'Impurity 2'], loc='upper left')
        # buf = io.BytesIO()
        # plt.savefig(buf, format="png")
        # data = base64.b64encode(buf.getbuffer()).decode("utf8")
        # plt.close()
        # img_path = "/Users/ke.zh/vLab-0.1.0/src/vLab/1.PNG"
        # image = base64.b64encode(open(img_path, "rb").read()).decode("ascii")
        return html.Div([
            # html.Img(src="data:img/png;base64, {}".format(image))
            dbc.Spinner(color="primary"),
            dbc.Progress(value=75, striped=True)
        ])



if __name__ == "__main__":
    app.run_server(port=8000)