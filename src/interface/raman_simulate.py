from vLab import Raman_Simulator
from numpy import save, load
import argparse

# def parse_args(x):
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--x', type=float, default=x, help='x')
#     args = parser.parse_args()
#     return args

def cal(a,b,c,d,e):
    rr = Raman_Simulator()
    rr.simulate(a,b,c,d,e)


if __name__ == '__main__':
    x = load('data_x.npy')
    interval = 0
    file = 0
    for interval in range(24):
        cal(x[:, 1:2][interval*10][0], x[:, 3:4][interval*10][0], x[:, 2:3][interval*10][0], x[:, 4:5][interval*10][0],file)
        file+=1