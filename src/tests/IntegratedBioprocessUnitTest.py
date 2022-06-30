import unittest
import numpy as np

from src.vLab.IntegratedBioprocess.Bioreactor import bioreactor
from src.vLab.IntegratedBioprocess.Chromatography import chromatography
from src.vLab.IntegratedBioprocess.Util import CellCultureModel, ChromatographyModel


class TestIntegratedBioprocess(unittest.TestCase):

    def test_chromatography_model(self):
        """ Set up Parameters """
        theta_C = ChromatographyModel()

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
        # Initial States
        X0 = 0.1  # initial viable biomass concentration (g/L)
        Sg0 = 40  # initial glycerol concentration (g/L)
        Sm0 = 0  # initial methanol concentration (g/L)
        P10 = 0  # initial product conentration (g/L)
        P20 = 0
        P30 = 0
        VB0 = 0.5  # initial bioreactor volume (L)
        VH0 = 1e-8  # initial hold tank volume (L)
        xC0 = [0] * (10 * theta_C.n + 3)  # np.zeros((10 * theta_C.n + 3, 1))

        x0 = [X0, Sg0, Sm0, P10, P20, P30, VB0, P10, P20, P30, VH0] + xC0

        ## Test
        # Time
        t0 = 0
        tg1 = 22  # glycerol batch period (h)
        tg2 = 10  # glycerol perfusion period (h)
        tm1 = 8  # methanol perfusion period (h)
        tm2 = 20  # methanol perfusion period (h)
        tl = 3  # load period (h)
        tw = 1  # wash period (h)
        te = 6  # elute period (h)
        rep = 3
        tran = np.cumsum([t0, tg1, tg2, tm1, tm2] + ([tl, tw, te] * rep))
        u_F = np.array([u_Fg1, u_Fg2, u_Fm1, u_Fm2] + [u_Fl, u_Fw, u_Fe] * rep).T
        u_Cin = np.array([u_Cing1, u_Cing2, u_Cinm1, u_Cinm2] + [u_Cinl, u_Cinw, u_Cine] * rep).T

        '''Test Case 1'''
        t = t0
        x = x0
        u_F = u_F[:, 10]
        u_Cin = u_Cin[:, 10]
        u_FC = u_F[3] + u_F[4] + u_F[5]
        u_Cin = np.concatenate([np.array(x[7:10]) * u_F[3], [u_Cin[2] * u_F[5]]]) / u_FC
        u_Cin[np.isnan(u_Cin)] = 0
        u_Cin = u_Cin + 1
        dxCdt = chromatography(t, x[11:], theta_C, np.array([u_FC, u_F[6]]), u_Cin)
        ground_truth_from_MIT_matlab_code = ([611.4649681528662] * 4) + ([0] * 299)
        self.assertTrue(
            np.allclose(dxCdt, ground_truth_from_MIT_matlab_code, rtol=5e-10, atol=5e-8),
            msg="chromatography calculation error"
        )

        '''Test Case 2'''
        u_Cin = np.array([1, 2, 3, 4])
        dxCdt = chromatography(t, x[11:], theta_C, np.array([u_FC, 0.1]), u_Cin)
        ground_truth_from_MIT_matlab_code = ([611.4649681528662,
                                              1.222929936305732e3,
                                              1.834394904458599e3,
                                              2.445859872611465e3]) + ([0] * 299)
        self.assertTrue(
            np.allclose(dxCdt, ground_truth_from_MIT_matlab_code, rtol=5e-10, atol=5e-8),
            msg="chromatography calculation error"
        )

        '''Test Case 4'''
        ground_truth_from_MIT_matlab_code = np.loadtxt('./data/chromagraphy_test_data.csv')
        x = [i + 1 for i in x]
        dxCdt = chromatography(t, x[11:], theta_C, np.array([u_FC, 0.1]), u_Cin)
        self.assertTrue(
            np.allclose(dxCdt, ground_truth_from_MIT_matlab_code, rtol=5e-10, atol=5e-8),
            msg="chromatography calculation error"
        )
        '''Test Case 5'''
        x_data = np.loadtxt('./data/x.txt', delimiter=',')
        x = x_data[1500, :]
        _u_F = np.array([0.03, 0.03, 0., 0., 0., 0.06, 0.06])
        _u_FC = 0.06
        _u_Cin = np.array([0, 0, 0, 1])
        ground_truth_from_MIT_matlab_code = np.loadtxt('./data/dxCdt.txt')
        dxCdt = chromatography(t, x[11:], theta_C, np.array([_u_FC, _u_F[6]]), _u_Cin)
        self.assertTrue(
            np.allclose(dxCdt, ground_truth_from_MIT_matlab_code, rtol=5e-10, atol=5e-8),
            msg="chromatography calculation error"
        )
        '''Test Case 6'''
        x = x_data[2000, :]
        _u_F = np.array([0.03, 0.03, 0.3, 0.5, 0.1, 0.06, 0.06])
        _u_FC = 0.1
        _u_Cin = np.array([0.5, 1, 0.1, 1])
        ground_truth_from_MIT_matlab_code = np.loadtxt('./data/dxCdt2.txt')
        dxCdt = chromatography(t, x[11:], theta_C, np.array([_u_FC, _u_F[6]]), _u_Cin)
        self.assertTrue(
            np.allclose(dxCdt, ground_truth_from_MIT_matlab_code, rtol=5e-10, atol=5e-8),
            msg="chromatography calculation error"
        )

    def test_bioreactor_model(self):
        theta = CellCultureModel()
        theta_C = ChromatographyModel()
        F0 = 0.5 * 60 / 1000  # typical flow rate (L/h)
        Sin_g0 = 80  # inlet glycerol concentration (g/L)
        Sin_m0 = 40  # inlet methanol concentration (g/L)

        # CPP
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

        # Initial States
        X0 = 0.1  # initial viable biomass concentration (g/L)
        Sg0 = 40  # initial glycerol concentration (g/L)
        Sm0 = 0  # initial methanol concentration (g/L)
        P10 = 0  # initial product conentration (g/L)
        P20 = 0
        P30 = 0
        VB0 = 0.5  # initial bioreactor volume (L)
        VH0 = 1e-8  # initial hold tank volume (L)
        xC0 = [0] * (10 * theta_C.n + 3)  # np.zeros((10 * theta_C.n + 3, 1))

        x0 = [X0, Sg0, Sm0, P10, P20, P30, VB0, P10, P20, P30, VH0] + xC0

        ## Test
        # Time
        t0 = 0
        tg1 = 22  # glycerol batch period (h)
        tg2 = 10  # glycerol perfusion period (h)
        tm1 = 8  # methanol perfusion period (h)
        tm2 = 20  # methanol perfusion period (h)
        tl = 3  # load period (h)
        tw = 1  # wash period (h)
        te = 6  # elute period (h)
        rep = 3
        tran = np.cumsum([t0, tg1, tg2, tm1, tm2] + ([tl, tw, te] * rep))
        u_F = np.array([u_Fg1, u_Fg2, u_Fm1, u_Fm2] + [u_Fl, u_Fw, u_Fe] * rep).T
        u_Cin = np.array([u_Cing1, u_Cing2, u_Cinm1, u_Cinm2] + [u_Cinl, u_Cinw, u_Cine] * rep).T

        t = t0
        x = x0
        u_F = u_F[:, 1]
        u_Cin = u_Cin[:, 1]
        dxBdt = bioreactor(t, x[:7], theta, np.array([u_F[0], u_F[1] + u_F[2]]), u_Cin[:2])
        ground_truth_from_MIT_matlab_code = [0.0258, 2.3631, 0, 0, 0.0000, 0.0003, 0]
        '''Test Case 1'''
        self.assertTrue(
            np.allclose(dxBdt, ground_truth_from_MIT_matlab_code, rtol=5e-5, atol=5e-4),
            msg="bioreactpr gradient calculation error"
        )

        '''Test Case 2'''
        x = [i + 5 for i in x]
        dxBdt = bioreactor(t, x[:7], theta, np.array([u_F[0], u_F[1] + u_F[2]]), u_Cin[:2])
        ground_truth_from_MIT_matlab_code = [2.3201, -1.6919, -2.8773, 0.0729, -0.0250, -0.0041, 0]
        self.assertTrue(
            np.allclose(dxBdt, ground_truth_from_MIT_matlab_code, rtol=5e-5, atol=5e-4),
            msg="bioreactpr gradient calculation error"
        )


if __name__ == '__main__':
    unittest.main()
