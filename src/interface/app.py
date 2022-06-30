import sys

import numpy as np

sys.path.append('/Users/ke.zh/vLab-0.1.0/src')

import pathlib
from collections import deque
import plotly
import dash
from dash import html, dcc, callback_context
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import dash_daq as daq
import os
import dash_bootstrap_components as dbc
from numpy import load, save
import base64
import io
from scipy.interpolate import interp1d
import subprocess
from subprocess import Popen, PIPE

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
app.config["suppress_callback_exceptions"] = True

APP_PATH = str(pathlib.Path(__file__).parent.resolve())

# app.css.append_css({
#     'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
# })

def build_banner():
    return html.Div(
        id="banner",
        className="banner",
        children=[
            html.Div(
                id="banner-text",
                children=[
                    html.H5("NIIMBL")
                ],
            ),
        ],
    )


def build_tabs():
    return html.Div(
        id="tabs",
        className="tabs",
        children=[
            dcc.Tabs(
                id="app-tabs",
                value="tab2",
                className="custom-tabs",
                children=[
                    dcc.Tab(
                        id="Specs-tab",
                        label="Variable Settings",
                        value="tab1",
                        className="custom-tab",
                        selected_className="custom-tab--selected",
                    ),
                    dcc.Tab(
                        id="Control-chart-tab",
                        label="Charts",
                        value="tab2",
                        className="custom-tab",
                        selected_className="custom-tab--selected",
                    ),
                ],
            )
        ],
    )



input_X0 = dcc.Slider(0, 5, 0.1, value=0.1, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id="input-X0")
input_Sg0 = dcc.Slider(30, 50, 1, value=40, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-Sg0')
input_Sm0 = dcc.Slider(0, 50, 1, value=10, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-Sm0')
input_La0 = dcc.Slider(0, 50, 1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-La0')
input_Am0 = dcc.Slider(0, 50, 1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-Am0')
input_P10 = dcc.Slider(0, 5, 0.1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-P10')
input_P20 = dcc.Slider(0, 5, 0.1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-P20')
input_P30 = dcc.Slider(0, 5, 0.1, value=0, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-P30')
input_VB0 = dcc.Slider(0, 5, 0.1, value=0.5, marks=None, tooltip={"placement": "bottom", "always_visible": True}, id='input-VB0')
input_VH0 = dcc.Input(id='input-VH0', type="text", value=1e-8)


def build_tab_1():
    return [
        html.Div(
            id="settings-menu",
            children=[
                html.Div(
                    id="value-setter-menu",
                    children=[
                        html.Div([
                            build_value_setter_line(
                                "value-setter-panel-header",
                                "Variables",
                                "Set value",
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial viable biomass concentration",
                                input_X0,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial glucose concentration",
                                input_Sg0,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial glutamine concentration",
                                input_Sm0,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial lactate concentration",
                                input_La0,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial ammonia concentration",
                                input_Am0,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial product conentration",
                                input_P10,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial impurity1 conentration",
                                input_P20,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial impurity2 conentration",
                                input_P30,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial bioreactor volume",
                                input_VB0,
                            ),
                            build_value_setter_line(
                                "value-setter-panel-ucl",
                                "initial hold tank volume",
                                input_VH0,
                            ),
                        ]),
                        html.Div(
                            id="button-div",
                            children=[
                                html.Button("Update", id="value-btn"),
                            ],
                        ),
                        html.Div(
                            id="value-setter-view-output", className="output-datatable"
                        ),
                    ],
                ),
            ],
        ),
    ]



def build_value_setter_line(line_num, label, col3):
    return html.Div(
        id=line_num,
        children=[
            html.Label(label, className="four columns"),
            # html.Label(value, className="four columns"),
            html.Div(col3, className="four columns"),
        ],
        className="row",
    )



def generate_section_banner(title):
    return html.Div(className="section-banner", children=title)


def build_top_panel():
    return html.Div(
        id="top-section-container",
        className="row",
        children=[
            # Metrics summary
            html.Div(
                id="metric-summary-session",
                className="three columns",
                children=[
                    generate_section_banner("Simulate"),
                    html.Div(
                        id="metric-div",
                        children=[
                            html.P("0-240h", id="period1"),
                            html.P("240-250h", id="period2"),
                            html.P("250-260h", id="period3"),
                            html.P("260-270h", id="period4"),
                            dbc.Progress(id="progress1", color="#a19d9d", striped=True),
                            dbc.Progress(id="progress2", color="#a19d9d", striped=True),
                            dbc.Progress(id="progress3", color="#a19d9d", striped=True),
                            dbc.Progress(id="progress4", color="#a19d9d", striped=True),
                            html.Button("Simulate", id="simulate_button", n_clicks=0),
                            html.Button("Restart", id="restart_button", n_clicks=0),
                            daq.StopButton(id="stop_button", size=160, n_clicks=0)
                        ],
                    ),
                ],
            ),
            # Piechart
            html.Div(
                id="ooc-piechart-outer",
                className="nine columns",
                children=[
                    generate_section_banner("Progress"),
                    build_progress_panel()
                ],
            ),
        ],
    )


def build_chart_panel():
    return html.Div(
        id="control-chart-container",
        className="twelve columns",
        children=[
            generate_section_banner("Live Charts"),

            dcc.Graph(id='live_raman_graph'),

            html.Div([
                html.Div([
                    dcc.Graph(id='live_bio_graph', animate=True),
                ], className="six columns"),
                html.Div([
                    dcc.Graph(id='live_hold_graph', animate=True),
                ], className="six columns"),
            ], className="row"),

            html.Div([
                html.Div([
                    dcc.Graph(id='live_cap_graph', animate=True),
                ], className="six columns"),
                html.Div([
                    dcc.Graph(id='live_po_graph', animate=True),
                ], className="six columns"),
            ], className="row"),


            dcc.Interval(
                id='graph-update',
                interval=1000,
                n_intervals=0
            ),
        ],
    )


app.layout = html.Div(
    id="big-app-container",
    children=[
        build_banner(),

        dcc.Interval(id="interval-component", n_intervals=0, interval=5000, disabled=True),
        dcc.Interval(id="progress_bar1", n_intervals=0, interval=500, disabled=True),
        dcc.Interval(id="progress_bar2", n_intervals=0, interval=2000, disabled=True),
        dcc.Interval(id="progress_bar3", n_intervals=0, interval=2500, disabled=True),
        dcc.Interval(id="progress_bar4", n_intervals=0, interval=3000, disabled=True),
        html.Div(
            id="app-container",
            children=[
                build_tabs(),
                html.Div(id="app-content"),
            ],
        ),
        dcc.Store(id="n-interval-stage", data=50)
    ],
)

buf = io.BytesIO()
data = base64.b64encode(buf.getbuffer()).decode("utf8")
img_path = "/Users/ke.zh/vLab-0.1.0/src/1.PNG"
image = base64.b64encode(open(img_path, "rb").read()).decode("ascii")

def build_progress_panel():
    return html.Div(
        id="progress-section-container",
        className="row",
        children=[
            html.Div(
                id="progress-summary-session",
                className="eight columns",
                children=[
                    html.Div(
                        id="metric-div",
                        children=[
                            html.Div(
                                id="progress-rows",
                                children=[
                                    html.Img(src="data:img/png;base64, {}".format(image), width=950, height=200)
                                ]
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )


@app.callback(
    [Output("app-content", "children"), Output("interval-component", "n_intervals")],
    [Input("app-tabs", "value")],
    [State("n-interval-stage", "data")],
)
def render_tab_content(tab_switch, stopped_interval):
    if tab_switch == "tab1":
        return build_tab_1(), stopped_interval
    return (
        html.Div(
            id="status-container",
            children=[
                html.Div(
                    id="graphs-container",
                    children=[build_top_panel(), build_chart_panel()],
                ),
            ],
        ),
        stopped_interval,
    )


@app.callback(
    Output("value-btn", "text"),
    [Input("value-btn", "n_clicks"),
     Input("input-X0", "value"),
     Input("input-Sg0", "value"),
     Input("input-Sm0", "value"),
     Input("input-La0", "value"),
     Input("input-Am0", "value"),
     Input("input-P10", "value"),
     Input("input-P20", "value"),
     Input("input-P30", "value"),
     Input("input-VB0", "value"),
     Input("input-VH0", "value")]
)
def set_value_setter_store(n_clicks, X0, Sg0, Sm0, La0, Am0, P10, P20, P30, VB0, VH0):
    if n_clicks != 0:
        save('value.npy', [X0, Sg0, Sm0, La0, Am0, P10, P20, P30, VB0, VH0])
        return "Saved"


# start, stop
@app.callback(
    [Output("graph-update", "disabled"), Output("stop_button", "buttonText")],
    [Input("stop_button", "n_clicks"), Input("graph-update", "n_intervals")],
    [State("graph-update", "disabled")],
)
def stop_production(n_clicks, n_intervals, current):
    changed_id = [p['prop_id'] for p in callback_context.triggered][0]
    if "stop_button" in changed_id:
        save('interval.npy', n_intervals)
        return not current, "stop" if current else "start"
    if n_clicks == 0 or n_intervals == 239 or n_intervals == 249 or n_intervals == 259 or n_intervals == 269:
        save('interval.npy', n_intervals)
        return True, "start"
    return False, "stop"


@app.callback(
    Output('progress_bar1', 'disabled'),
    [Input('simulate_button', 'n_clicks'), Input('progress1', 'value')]
)
def plot_input(n_clicks, value):
    changed_id = [p['prop_id'] for p in callback_context.triggered][0]
    if n_clicks != 0:
        if value != 100:
            return False
        else:
            return True
    return True

@app.callback(
    Output('simulate_button', 'text'),
    [Input('simulate_button', 'n_clicks')]
)
def plot_input(n_clicks):
    if n_clicks != 0:
        subprocess.Popen('python3 /Users/ke.zh/vlab/src/interface/raman_simulate.py', shell=True)
        print("DoneP")
        cmd = 'python3 /Users/ke.zh/vlab/src/interface/simulate.py'
        os.system(cmd)
        print("DoneProgress")
        return None


# @app.callback(
#     Output("progress_bar1", "n_intervals"),
#     [Input('restart_button', 'n_clicks')]
# )
# def plot_input(n_clicks):
#     if n_clicks != 0:
#         return 0


@app.callback(
    Output('progress_bar2', 'disabled'),
    [Input('progress1', 'value'), Input('progress2', 'value')]
)
def plot_input(value1, value2):
    if value1 == 100 and value2 != 100:
        save('period.npy', 1)
        return False
    return True

@app.callback(
    Output('progress_bar3', 'disabled'),
    [Input('progress2', 'value'), Input('progress3', 'value')]
)
def plot_input(value2, value3):
    if value2 == 100 and value3 != 100:
        save('period.npy', 2)
        return False
    return True

@app.callback(
    Output('progress_bar4', 'disabled'),
    [Input('progress3', 'value'), Input('progress4', 'value')]
)
def plot_input(value3, value4):
    if value3 == 100 and value4 != 100:
        save('period.npy', 3)
        return False
    if value4 == 100:
        save('period.npy', 4)
    return True

@app.callback(
    [Output("progress1", "value"), Output("progress1", "label")],
     Input("progress_bar1", "n_intervals")
)
def update_progress(n):
    progress = min(n % 110, 100)
    return progress, f"{progress} %"

@app.callback(
    [Output("progress2", "value"), Output("progress2", "label")],
     Input("progress_bar2", "n_intervals")
)
def update_progress(n):
    progress = min(n % 110, 100)
    return progress, f"{progress} %"

@app.callback(
    [Output("progress3", "value"), Output("progress3", "label")],
     Input("progress_bar3", "n_intervals")
)
def update_progress(n):
    progress = min(n % 110, 100)
    return progress, f"{progress} %"

@app.callback(
    [Output("progress4", "value"), Output("progress4", "label")],
     Input("progress_bar4", "n_intervals")
)
def update_progress(n):
    progress = min(n % 110, 100)
    return progress, f"{progress} %"


k=0
g=0
period=0
file=0
X1 = deque(maxlen=5000)
X2 = deque(maxlen=5000)

new_Y1 = deque(maxlen=5000)
new_Y2 = deque(maxlen=5000)
new_Y3 = deque(maxlen=5000)
new_Y4 = deque(maxlen=5000)
new_Y5 = deque(maxlen=5000)
new_Y6 = deque(maxlen=5000)
new_Y7 = deque(maxlen=5000)
new_Y8 = deque(maxlen=5000)
Y4 = deque(maxlen=5000)
Y5 = deque(maxlen=5000)
Y6 = deque(maxlen=5000)
Y7 = deque(maxlen=5000)
Y8 = deque(maxlen=5000)
Y9 = deque(maxlen=5000)
Y10 = deque(maxlen=5000)
Y11 = deque(maxlen=5000)
Y12 = deque(maxlen=5000)


def generate_data(period):
    if period == 0:
        data_t = load('init_data_t.npy')
        data_x = load('init_data_x.npy')
    else:
        data_t = load('data_t.npy')
        data_x = load('data_x.npy')

    data_tC = load('data_tC.npy')
    data_yplot = load('data_yplot.npy')

    new_array1 = data_x[:, 0:1]
    new_array2 = data_x[:, 1:2]
    new_array3 = data_x[:, 2:3]
    new_array4 = data_x[:, 3:4]
    new_array5 = data_x[:, 4:5]
    new_array6 = data_x[:, 5:6]
    new_array7 = data_x[:, 6:7]
    new_array8 = data_x[:, 7:8]
    array4 = data_x[:, 9:10]
    array5 = data_x[:, 10:11]
    array6 = data_x[:, 11:12]

    new_clean_array1 = [i for arr in new_array1 for i in arr]
    new_clean_array2 = [i for arr in new_array2 for i in arr]
    new_clean_array3 = [i for arr in new_array3 for i in arr]
    new_clean_array4 = [i for arr in new_array4 for i in arr]
    new_clean_array5 = [i for arr in new_array5 for i in arr]
    new_clean_array6 = [i for arr in new_array6 for i in arr]
    new_clean_array7 = [i for arr in new_array7 for i in arr]
    new_clean_array8 = [i for arr in new_array8 for i in arr]
    clean_array4 = [i for arr in array4 for i in arr]
    clean_array5 = [i for arr in array5 for i in arr]
    clean_array6 = [i for arr in array6 for i in arr]

    array7 = data_yplot[:, 0, -1]
    array8 = data_yplot[:, 1, -1]
    array9 = data_yplot[:, 2, -1]
    array10 = data_yplot[:, 4, -1]
    array11 = data_yplot[:, 5, -1]
    array12 = data_yplot[:, 6, -1]

    new_os1 = interp1d(data_t, new_clean_array1)
    new_os2 = interp1d(data_t, new_clean_array2)
    new_os3 = interp1d(data_t, new_clean_array3)
    new_os4 = interp1d(data_t, new_clean_array4)
    new_os5 = interp1d(data_t, new_clean_array5)
    new_os6 = interp1d(data_t, new_clean_array6)
    new_os7 = interp1d(data_t, new_clean_array7)
    new_os8 = interp1d(data_t, new_clean_array8)
    os4 = interp1d(data_t, clean_array4)
    os5 = interp1d(data_t, clean_array5)
    os6 = interp1d(data_t, clean_array6)
    os7 = interp1d(data_tC, array7)
    os8 = interp1d(data_tC, array8)
    os9 = interp1d(data_tC, array9)
    os10 = interp1d(data_tC, array10)
    os11 = interp1d(data_tC, array11)
    os12 = interp1d(data_tC, array12)

    if period == 1:
        data_t = np.linspace(0, 240, 241)
    if period == 2:
        data_t = np.linspace(0, 250, 251)
        data_tC = np.linspace(241, 250, 500)
    if period == 3:
        data_t = np.linspace(0, 260, 261)
        data_tC = np.linspace(241, 260, 1000)
    if period == 4:
        data_t = np.linspace(0, 270, 271)
        data_tC = np.linspace(241, 270, 1500)

    new_clean_array1 = new_os1(data_t)
    new_clean_array2 = new_os2(data_t)
    new_clean_array3 = new_os3(data_t)
    new_clean_array4 = new_os4(data_t)
    new_clean_array5 = new_os5(data_t)
    new_clean_array6 = new_os6(data_t)
    new_clean_array7 = new_os7(data_t)
    new_clean_array8 = new_os8(data_t)

    clean_array4 = os4(data_t)
    clean_array5 = os5(data_t)
    clean_array6 = os6(data_t)
    array7 = os7(data_tC)
    array8 = os8(data_tC)
    array9 = os9(data_tC)
    array10 = os10(data_tC)
    array11 = os11(data_tC)
    array12 = os12(data_tC)

    return data_t, data_tC, new_clean_array1, new_clean_array2, new_clean_array3, new_clean_array4, new_clean_array5, \
           new_clean_array6, new_clean_array7, new_clean_array8, clean_array4, clean_array5, clean_array6, array7, \
           array8, array9, array10, array11, array12

@app.callback(
    [Output('live_bio_graph', 'figure'),
     Output('live_hold_graph', 'figure'),
     Output('live_cap_graph', 'figure'),
     Output('live_po_graph', 'figure'),
     Output('live_raman_graph', 'figure')],
     Input('graph-update', 'n_intervals')
)
def update_graph_scatter(n):
    global k, period, file, g

    period = load('period.npy', allow_pickle=True)
    a, b, new_array1, new_array2, new_array3, new_array4, new_array5, new_array6, new_array7, new_array8,\
    array4, array5, array6, array7, array8, array9, array10, array11, array12 = generate_data(period)

    X1.append(a[k])

    # if k<241:
    #     X2.append(k)
        # Y7.append(0)
        # Y8.append(0)
        # Y9.append(0)
        # Y10.append(0)
        # Y11.append(0)
        # Y12.append(0)
    if k > 240:
        # X2.append(b[k-241])
        # Y7.append(array7[k-241])
        # Y8.append(array8[k-241])
        # Y9.append(array9[k-241])
        # Y10.append(array10[k-241])
        # Y11.append(array11[k-241])
        # Y12.append(array12[k-241])
        for i in range(0, 50):
            X2.append(b[g])
            Y7.append(array7[g])
            Y8.append(array8[g])
            Y9.append(array9[g])
            Y10.append(array10[g])
            Y11.append(array11[g])
            Y12.append(array12[g])
            g+=1

    new_Y1.append(new_array1[k])
    new_Y2.append(new_array2[k])
    new_Y3.append(new_array3[k])
    new_Y4.append(new_array4[k])
    new_Y5.append(new_array5[k])
    new_Y6.append(new_array6[k])
    new_Y7.append(new_array7[k])
    new_Y8.append(new_array8[k])

    Y4.append(array4[k])
    Y5.append(array5[k])
    Y6.append(array6[k])

    new_data = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y1),
        name='Cell Density',
        mode='lines+markers'
    )
    new_data2 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y2),
        name='glucose',
        mode='lines+markers'
    )
    new_data3 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y3),
        name='glutamine',
        mode='lines+markers'
    )
    new_data4 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y4),
        name='lactate',
        mode='lines+markers'
    )
    new_data5 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y5),
        name='ammonium',
        mode='lines+markers'
    )
    new_data6 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y6),
        name='Product',
        mode='lines+markers'
    )
    new_data7 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y7),
        name='Impurity1',
        mode='lines+markers'
    )
    new_data8 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(new_Y8),
        name='Impurity2',
        mode='lines+markers'
    )
    data4 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(Y4),
        name='Product',
        mode='lines+markers'
    )
    data5 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(Y5),
        name='Impurity1',
        mode='lines+markers'
    )
    data6 = plotly.graph_objs.Scatter(
        x=list(X1),
        y=list(Y6),
        name='Impurity2',
        mode='lines+markers'
    )
    data7 = plotly.graph_objs.Scatter(
        x=list(X2),
        y=list(Y7),
        name='Product',
        mode='lines+markers'
    )
    data8 = plotly.graph_objs.Scatter(
        x=list(X2),
        y=list(Y8),
        name='Impurity1',
        mode='lines+markers'
    )
    data9 = plotly.graph_objs.Scatter(
        x=list(X2),
        y=list(Y9),
        name='Impurity2',
        mode='lines+markers'
    )
    data10 = plotly.graph_objs.Scatter(
        x=list(X2),
        y=list(Y10),
        name='Product',
        mode='lines+markers'
    )
    data11 = plotly.graph_objs.Scatter(
        x=list(X2),
        y=list(Y11),
        name='Impurity1',
        mode='lines+markers'
    )
    data12 = plotly.graph_objs.Scatter(
        x=list(X2),
        y=list(Y12),
        name='Impurity2',
        mode='lines+markers'
    )

    file = min(int(n/10), 23)
    predicted = load(f'predicted:{file}.npy'.format(file=file))
    new_clean_array = [i for arr in predicted for i in arr]
    data = plotly.graph_objs.Scatter(
        x=list(np.linspace(400, 3000, 2601)),
        y=list(new_clean_array),
        name='Ramen',
        mode='lines+markers'
    )
    k+=1
    return {'data': [new_data, new_data2, new_data3, new_data4, new_data5, new_data6, new_data7, new_data8], 'layout': go.Layout(title='Bioreactor', paper_bgcolor='#ffffff', plot_bgcolor='#ffffff', xaxis=dict(range=[min(X1), max(X1)]), yaxis=dict(range=[0, 70]))},\
           {'data': [data4, data5, data6], 'layout': go.Layout(title='Hold Tank', paper_bgcolor='#ffffff', plot_bgcolor='#ffffff', xaxis=dict(range=[min(X1), max(X1)]), yaxis=dict(range=[0, 0.8]))}, \
           {'data': [data7, data8, data9], 'layout': go.Layout(title='Capture', paper_bgcolor='#ffffff', plot_bgcolor='#ffffff', xaxis=dict(range=[240, 270]), yaxis=dict(range=[0, 0.8]))},\
           {'data': [data10, data11, data12], 'layout': go.Layout(title='Polish', paper_bgcolor='#ffffff', plot_bgcolor='#ffffff', xaxis=dict(range=[240, 270]), yaxis=dict(range=[0, 0.8]))},\
           {'data': [data], 'layout': go.Layout(title='Raman:{time}h'.format(time=file*10), paper_bgcolor='#ffffff', plot_bgcolor='#ffffff', xaxis=dict(range=[400, 3000]), yaxis=dict(range=[-0.5, 1.5]))}

# @app.callback(
#     Output('live_raman_graph', 'figure'),
#     [Input('show_button', 'n_clicks')]
# )
# def plot_input(n_clicks):
#     subprocess.Popen('python3 /Users/ke.zh/vlab/src/interface/raman_simulate.py', shell=True)
#     predicted = load('predicted.npy')
#     new_clean_array = [i for arr in predicted for i in arr]
#     data = plotly.graph_objs.Scatter(
#         x=list(np.linspace(400, 3000, 2601)),
#         y=list(new_clean_array),
#         name='Impurity2',
#         mode='lines+markers'
#     )
#     print("Done")
#     return {'data': [data], 'layout': go.Layout(title='Raman', paper_bgcolor='#ffffff', plot_bgcolor='#ffffff', xaxis=dict(range=[400, 3000]), yaxis=dict(range=[-0.2, 0.8]))}


# Running the server
if __name__ == "__main__":
    save('period.npy', 0)
    app.run_server(debug=True, port=8061)
