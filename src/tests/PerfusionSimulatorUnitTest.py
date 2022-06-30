import unittest
import numpy as np

from vLab.GlycosylationModelBase.GlycosylationDerivative import steady_state_inner_derivative
from vLab.GlycosylationModelBase.GlycosylationModelParams import GlycosylationModelParamClass, CellCultureVariables
from vLab.GlycosylationModelBase.GlycosylationNetwork import GlycosylationNetwork
from src.vLab.PerfusionSimulator.GlycosylationODESolver import ODESolver


def compute_derivative(y, z, x, p, fp):
    deriv = steady_state_inner_derivative(z, y, x, p, fp)
    return np.round(deriv, 4)


class TestPerfusionSimulator(unittest.TestCase):

    def test_glycosylation_inside_Golgi_apparatus(self):
        fp = GlycosylationNetwork(network_data_path='data/Network Description.csv')
        p = GlycosylationModelParamClass()

        x = CellCultureVariables(1.5, 0.01, 0.1198, 66.3856,
                                 np.array([0.490 + 1.452, 0.117 + 0.379, 0.058 + 0.190]) * 1e3,
                                 np.array([1.62, 0.043, 0.1158, 0.040]) * 1e3)
        # set up initial conditions
        ic = np.zeros((fp.nos + fp.nns + fp.nn))
        ic[0] = x.mabtiter  # umol / L
        ic[fp.nos: (fp.nos + fp.nns)] = x.nscyt * 40  # nucleotide sugar concentrations in umol / L.third entry is mystery
        ic[fp.nos + 3] = x.udpgalcyt * 1e3 * 40  # updating with correct UDP-Gal concentration
        ic[(fp.nos + fp.nns):] = x.ncyt  # sum of nucleotide concentrations in umol / L

        # z = 0.1
        ground_truth_from_MIT_matlab_code = [-657.8735, 657.8735, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.7306, 0.0111, 0.0000, 0.0000,
                                             -2.7306, -0.0111, -0.0000]
        deriv = compute_derivative(ic, 0.1, x, p, fp)
        self.assertTrue(
            np.allclose(deriv, ground_truth_from_MIT_matlab_code),
            msg="derivative calculation error"
        )

        # z = 0.1
        ground_truth_from_MIT_matlab_code = [-35.4457, 35.4457, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0035, 0.0000, 0.0000, 0.0000,
                                             -0.0035, -0.0000, -0.0000]

        deriv = compute_derivative(ic, 0.01, x, p, fp)
        self.assertTrue(
            np.allclose(deriv, ground_truth_from_MIT_matlab_code),
            msg="derivative calculation error"
        )

    def test_glycosylation_inside_Golgi_apparatus_cell_culture_variable(self):
        fp = GlycosylationNetwork(network_data_path='data/Network Description.csv')
        p = GlycosylationModelParamClass()

        # x = 30
        x = CellCultureVariables(1.5, 0.01, 0.1198, 30,
                                 np.array([0.490 + 1.452, 0.117 + 0.379, 0.058 + 0.190]) * 1e3,
                                 np.array([1.62, 0.043, 0.1158, 0.040]) * 1e3)
        ic = np.zeros((fp.nos + fp.nns + fp.nn))
        ic[0] = x.mabtiter  # umol / L
        ic[fp.nos: (fp.nos + fp.nns)] = x.nscyt * 40  # nucleotide sugar concentrations in umol / L.third entry is mystery
        ic[fp.nos + 3] = x.udpgalcyt * 1e3 * 40  # updating with correct UDP-Gal concentration
        ic[(fp.nos + fp.nns):] = x.ncyt  # sum of nucleotide concentrations in umol / L

        ground_truth_from_MIT_matlab_code = [-416.8248, 416.8248, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.7306, 0.0111, 0.0000, 0.0000,
                                             -2.7306, -0.0111, -0.0000]

        deriv = compute_derivative(ic, 0.1, x, p, fp)
        self.assertTrue(
            np.allclose(deriv, ground_truth_from_MIT_matlab_code),
            msg="derivative calculation error"
        )

    def test_ode_solver_for_glycosylation(self):
        fp = GlycosylationNetwork(network_data_path='data/Network Description.csv')
        p = GlycosylationModelParamClass()
        x = CellCultureVariables(1.5, 0.01, 0.1198, 66.3856,
                                 np.array([0.490 + 1.452, 0.117 + 0.379, 0.058 + 0.190]) * 1e3,
                                 np.array([1.62, 0.043, 0.1158, 0.040]) * 1e3)
        # compute boundary conditions
        ic = np.zeros((fp.nos + fp.nns + fp.nn))
        ic[0] = x.mabtiter  # umol / L
        ic[fp.nos:(fp.nos + fp.nns)] = x.nscyt * 40  # nucleotide sugar concentrations in umol / L.third entry is mystery
        ic[fp.nos + 3] = x.udpgalcyt * 1e3 * 40  # updating with correct UDP-Gal concentration
        ic[(fp.nos + fp.nns):] = x.ncyt  # sum of nucleotide concentrations in umol / L

        t = [0, 1] # np.linspace(0, 1, 1001)
        ode_solver = ODESolver(t, ic, x, p, fp)
        HM, FA1G1, FA2G0, FA2G1, FA2G2, SIA = ode_solver.solve()
        os = ode_solver.os[:, -1]

        ground_truth_from_MIT_matlab_code = [0.0000, 0.0004, 0.0002, 0.0012, 0.0312, 0.0009, 0.0003, 0.1478, 0.0000,
                                             0.1180, 0.0052, 0.0392, 0.0000, 25.7738, 0.0081, 0.0258, 0.0258, 0.0000,
                                             16.5879, 16.5879, 0.0001, 0.0102, 0.0003, 0.0003, 6.4873, 0.2133, 0.2133,
                                             0.0001, 0.0001, 0.0532, 0.0532, 0.0000, 0.0007]
        self.assertTrue(
            np.allclose(os, ground_truth_from_MIT_matlab_code, rtol=5e-1, atol=5e-1),
            msg="glycosylation calculation error"
        )


if __name__ == '__main__':
    unittest.main()
