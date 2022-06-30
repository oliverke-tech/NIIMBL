import unittest
import numpy as np
from scipy.integrate import odeint

from vLab.DynamicGlycosylationSimulator.OuterFunction import outerfunction
from vLab.GlycosylationModelBase.GlycosylationNetwork import GlycosylationNetwork
from vLab.GlycosylationModelBase.GlycosylationModelParams import GlycosylationModelParamClass, label_matrix, \
    default_initial_states


class TestDynamicGlycosylaationSimulator(unittest.TestCase):

    def test_outerfunction(self):
        inputs = np.loadtxt('data/inputs.csv', delimiter='\t')
        u_start = [1.25, 0.7, 0.09, 5]
        ninputs = 4
        size_t = 1200
        num_macro_states = 6
        T = 20  # sampling time [min]
        tott = 792000  # total time [min]
        nstep = round(tott / T)  # number of sampling instances

        ## Inputs (vary all of them during identification)
        nsplit = inputs.shape[0] - 1  # the last set of inputs is not simulated to obtain outputs

        q_SS = inputs[:, 1:]
        nu = inputs[:, 1:].shape[1]  # number of control inputs
        u_F = []
        for ii in range(inputs.shape[0] - 1):  # 2:inputs.shape[0]-1:
            u_F.append(np.tile(q_SS[ii, :], (int(nstep / nsplit) + 1, 1)))
        u_F = np.vstack(u_F)

        # Time steps
        dt = 20  # sampling time [min]
        T = 2160
        tsp = np.arange(0, T + 0.1, dt).astype(int)  # total MPC run time
        init = default_initial_states[:6]
        init_rest = default_initial_states[6:39]
        jj = 1

        fp = GlycosylationNetwork(network_data_path='data/Network Description.csv')
        p = GlycosylationModelParamClass(is_dynamic=True)
        ig = label_matrix[:, 1].astype('float')
        yout_macro = None
        for ii in range(len(tsp)):
            u_in = u_F[ii, :]
            feed_cond = np.hstack([np.reshape([0, 20], (2, 1)), np.tile(u_in, (2, 1))])
            yout_macro = outerfunction(ig, p, fp, dt, feed_cond, init, init_rest)
            ic_states = yout_macro[:6]
            init = ic_states  # Initialize next iteration
            init_rest = yout_macro[6:39]

        ground_truth_from_MIT_matlab_code = [1.80000000000000,
                                             63.6673873643804,
                                             0.161999999999992,
                                             4.76737310570132,
                                             0.409071217509809,
                                             8.99999999999949,
                                             3.59316705672924e-09,
                                             0.000218585701331633,
                                             8.50070321225160e-05,
                                             0.000621159720709343,
                                             0.00870003935363342,
                                             0.00985256318906216,
                                             0.00499745400444354,
                                             0.0886145000230394,
                                             0.000628571397580493,
                                             0.0749589251728266,
                                             0.0131830336600201,
                                             0.275036490586520,
                                             0.00741222006606196,
                                             8.78050573638822,
                                             0.0528003466483801,
                                             0.296705540927756,
                                             0.296705540927756,
                                             0.000132942801060332,
                                             5.92159363807140,
                                             5.92159363807140,
                                             0.000819895171492308,
                                             0.149966436372544,
                                             0.00403815591870216,
                                             0.00403815591870216,
                                             2.45772737199435,
                                             0.0676984330256695,
                                             0.0676984330256687,
                                             0.00122047929469295,
                                             0.00122047929469296,
                                             0.0180466641500125,
                                             0.0180466641500126,
                                             1.44852991572927e-05,
                                             0.000196722024564630]
        self.assertTrue(
            np.allclose(yout_macro, ground_truth_from_MIT_matlab_code, rtol=1e-4, atol=1e-4),
            msg="dynamic glycosylation calculation error"
        )


if __name__ == '__main__':
    unittest.main()
