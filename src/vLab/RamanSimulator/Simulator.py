import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import pickle
from numpy import load, save


class Raman_Simulator:
    """
    Performs Raman Spectra data simulation based on concentrations of molecules.
    """

    def __init__(self):
        self.concen_base = pd.Series(np.array([6.190, 0.000, 6.120, 0.200, 0.341]),
                                     index=['Glucose', 'Lactate', 'Glutamine', 'Ammonia', 'Glutamate'])
        # self.Raman_base = pd.read_csv('data/RamanBaseline', index_col=(0), header=None).squeeze()
        self.Raman_base = pd.read_csv('/Users/ke.zh/vlab/data/RamanBaseline', index_col=(0), header=None).squeeze()

    def _convert_unit(self, Glc, Lac, Gln, NH3, Glu):
        return Glc * 180.156, \
               Lac * 90.08, \
               Gln * 146.14, \
               NH3 * 17.031, \
               Glu * 147.13

    def simulate(self, Glc, Lac, Gln, NH3, file, Glu=-100.0, show_figure=False, converted=True):
        """
        :param float Glc: concentrations (g/l) of Glucose
        :param float Lac: concentrations (g/l) of Lactate
        :param float Gln: concentrations (g/l) of Glutamine
        :param float NH3: concentrations (g/l) of Ammonia
        :param float Glu: concentrations (g/l) of Glutamate
        :param bool show_figure: whether or not show the figure; default false.
        :return: Simulated Raman Spectra Data
        """
        if not converted:
            Glc, Lac, Gln, NH3, Glu = self._convert_unit(Glc, Lac, Gln, NH3, Glu)

        if Glu > 0:
            concen = pd.Series(np.array([Glc, Lac, Gln, NH3, Glu]),
                               index=['Glucose', 'Lactate', 'Glutamine', 'Ammonia', 'Glutamate'])
            # reg_model = pickle.load(open('src/vLab/RamanSimulator/simulator_wi_Glu.sav', 'rb'))
            reg_model = pickle.load(open('/Users/ke.zh/vlab/src/vLab/RamanSimulator/simulator_wi_Glu.sav', 'rb'))
            predicted_diff = pd.DataFrame(reg_model.predict(np.array(concen).reshape(1, -1)))
            predicted_diff.columns = self.Raman_base.index
            predicted = predicted_diff + self.Raman_base

        else:
            concen = pd.Series(np.array([Glc, Lac, Gln, NH3]), index=['Glucose', 'Lactate', 'Glutamine', 'Ammonia'])
            # reg_model = pickle.load(open('src/vLab/RamanSimulator/simulator_wo_Glu.sav', 'rb'))
            reg_model = pickle.load(open('/Users/ke.zh/vlab/src/vLab/RamanSimulator/simulator_wo_Glu.sav', 'rb'))
            predicted_diff = pd.DataFrame(reg_model.predict(np.array(concen).reshape(1, -1)))
            predicted_diff.columns = self.Raman_base.index
            predicted = predicted_diff + self.Raman_base

        save(f'predicted:{file}.npy'.format(file=file), np.transpose(predicted))

        if show_figure:
            import matplotlib.pyplot as plt
            plt.figure()
            np.transpose(predicted).plot(legend=False, title="Simulated Raman Spectra Data")
            plt.xlabel('RamanShift(cm-1)')
            plt.ylabel('Counts')
            plt.show()

        simulated_Raman = predicted.squeeze()
        return simulated_Raman
