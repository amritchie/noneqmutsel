# File: KimuraMutSelMatrix.py

# File: pyvolve_sim_seqs.py
#
# Simple class to allow use of different fixation rate formulas in Stephanie Spielman's pyvolve.

import pyvolve
import numpy as np

ZERO = 1e-8

class KimuraMutSelMatrix(pyvolve.MutSel_Matrix):

    def __init__(self, *args):
        super(KimuraMutSelMatrix, self).__init__(*args)
        try:
            self._N = self.params['popsize']
        except:
            self._N = 1.;
            
    def _calc_fixrate_fitness(self, source, target, parameters):
        sij = (parameters['fitness'][target]/parameters['fitness'][source]) - 1.
        if abs(sij) <= ZERO:
            fixation_rate = 1.0/self._N
        else:
            fixation_rate = (1. - np.exp(-sij))/(1 - np.exp(-2.*self._N*sij))
            if fixation_rate <= ZERO:
                fixation_rate = ZERO
        return fixation_rate
