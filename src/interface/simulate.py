import sys
sys.path.append('/Users/ke.zh/vLab-0.1.0/src')
import argparse
from vLab.IntegratedBioprocess.PlantwiseSimulator import PlantwiseSimulator
import numpy as np
from numpy import save, load
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def parse_args(value):
    parser = argparse.ArgumentParser()
    parser.add_argument('--X0', type=float, default=value[0], help='X0')
    parser.add_argument('--Sg0', type=float, default=value[1], help='Sg0')
    parser.add_argument('--Sm0', type=float, default=value[2], help='Sm0')
    parser.add_argument('--Sl0', type=float, default=value[3], help='Sl0')
    parser.add_argument('--Amm0', type=float, default=value[4], help='Amm0')
    parser.add_argument('--P10', type=float, default=value[5], help='P10')
    parser.add_argument('--P20', type=float, default=value[6], help='P20')
    parser.add_argument('--P30', type=float, default=value[7], help='P30')
    parser.add_argument('--VB0', type=float, default=value[8], help='VB0')
    parser.add_argument('--VH0', type=float, default=value[9], help='VH0')
    args = parser.parse_args()
    return args


def cal(X0, Sg0, Sm0, Sl0, Amm0, P10, P20, P30, VB0, VH0):

    x0 = [X0, Sg0, Sm0, Sl0, Amm0, P10, P20, P30, VB0, P10, P20, P30, VH0]
    xC0 = [0] * (10 * 30 + 3)
    x0 = x0 + xC0
    import time
    from vLab.IntegratedBioprocess.Util import CellCultureModel

    start_time = time.time()
    bioreactor_param = CellCultureModel()
    bioreactor_param.set_cho_cell_lines()
    t0 = 0  # initial time
    tg1 = 22 * 4  # glycerol batch period (h)
    tg2 = 10 * 4  # glycerol perfusion period (h)
    tm1 = 8 * 4  # methanol perfusion period (h)
    tm2 = 20 * 4  # methanol perfusion period (h)
    tl = 3  # load period (h)
    tw = 1  # wash period (h)
    te = 6  # elute period (h)
    rep = 3
    process_time = np.cumsum(
        [t0, tg1, tg2, tm1, tm2] + ([tl, tw, te] * rep))

    F0 = 0.5 * 60 / 1000  # typical flow rate (L/h)
    Sin_g0 = 80  # inlet glycerol concentration (g/L)
    Sin_m0 = 40  # inlet methanol concentration (g/L)
    u_Fg1 = [0, 0, 0, 0, 0, 0, 0]
    u_Cing1 = [0, 0, 0]  # glycerol batch
    u_Fg2 = [F0, 0, F0, 0, 0, 0, 0]
    u_Cing2 = [Sin_g0, 0, 0]  # glycerol perfusion to waste
    u_Fm1 = [F0, 0, F0, 0, 0, 0, 0]
    u_Cinm1 = [0, Sin_m0, 0]  # methanol perfusion to waste
    u_Fm2 = [F0, F0, 0, 0, 0, 0, 0]
    u_Cinm2 = [0, Sin_m0, 0]  # methanol perfusion to tank
    u_Fl = [F0, F0, 0, 2 * F0, 0, 0, 0]
    u_Cinl = [0, Sin_m0, 0]  # load
    u_Fw = [F0, F0, 0, 0, 2 * F0, 0, 0]
    u_Cinw = [0, Sin_m0, 0]  # wash
    u_Fe = [F0, F0, 0, 0, 0, 2 * F0, 2 * F0]
    u_Cine = [0, Sin_m0, 1]  # elute
    flow = np.array([u_Fg1, u_Fg2, u_Fm1, u_Fm2] + [u_Fl, u_Fw, u_Fe] * rep).T
    inlet = np.array([u_Cing1, u_Cing2, u_Cinm1, u_Cinm2] + [u_Cinl, u_Cinw, u_Cine] * rep).T

    solver = PlantwiseSimulator(bioreactor_param=bioreactor_param)
    solver._process_time = solver._process_time
    sol = solver.solve(x0, [0, 240], process_time=process_time, flow=flow, inlet=inlet)
    elapse_time_bioreactor = time.time() - start_time
    t = np.array(sol.t)
    x = np.array(sol.x)
    save('data_t.npy', t)
    save('data_x.npy', x)
    print("1")


    sol = solver.solve(sol.x[-1], [240, 250])
    t = np.array(sol.t)
    x = np.array(sol.x)
    tC = t[t >= solver._process_time[4]]
    nrows = len(tC)
    xC = x[t >= solver._process_time[4], 13:]
    yplot = xC[:, :(30 * 10)].reshape(nrows, 10, 30, order='F')
    old_t = load('data_t.npy')
    new_t = np.append(old_t, t)
    old_x = load('data_x.npy')
    new_x = np.array(old_x.tolist() + x.tolist())
    save('data_t.npy', new_t)
    save('data_x.npy', new_x)
    save('data_tC.npy', tC)
    save('data_yplot.npy', yplot)
    print("2")

    sol = solver.solve(sol.x[-1], [250, 260])
    t = np.array(sol.t)
    x = np.array(sol.x)
    tC = t[t >= solver._process_time[4]]
    nrows = len(tC)
    xC = x[t >= solver._process_time[4], 13:]
    yplot = xC[:, :(30 * 10)].reshape(nrows, 10, 30, order='F')
    old_t = load('data_t.npy')
    new_t = np.append(old_t, t)
    old_x = load('data_x.npy')
    new_x = np.array(old_x.tolist() + x.tolist())
    old_tC = load('data_tC.npy')
    new_tC = np.append(old_tC, tC)
    old_yplot = load('data_yplot.npy')
    new_yplot = np.array(old_yplot.tolist() + yplot.tolist())
    save('data_t.npy', new_t)
    save('data_x.npy', new_x)
    save('data_tC.npy', new_tC)
    save('data_yplot.npy', new_yplot)
    print("3")

    sol = solver.solve(sol.x[-1], [260, 270])
    t = np.array(sol.t)
    x = np.array(sol.x)
    tC = t[t >= solver._process_time[4]]
    nrows = len(tC)
    xC = x[t >= solver._process_time[4], 13:]
    yplot = xC[:, :(30 * 10)].reshape(nrows, 10, 30, order='F')
    old_t = load('data_t.npy')
    new_t = np.append(old_t, t)
    old_x = load('data_x.npy')
    new_x = np.array(old_x.tolist() + x.tolist())
    old_tC = load('data_tC.npy')
    new_tC = np.append(old_tC, tC)
    old_yplot = load('data_yplot.npy')
    new_yplot = np.array(old_yplot.tolist() + yplot.tolist())
    save('data_t.npy', new_t)
    save('data_x.npy', new_x)
    save('data_tC.npy', new_tC)
    save('data_yplot.npy', new_yplot)
    print("4")


if __name__ == '__main__':
    value = load('value.npy')
    args = parse_args(value)
    cal(args.X0, args.Sg0, args.Sm0, args.Sl0, args.Amm0, args.P10, args.P20, args.P30, args.VB0, args.VH0)
