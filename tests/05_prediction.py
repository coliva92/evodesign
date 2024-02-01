from evodesign.Prediction.ESMFoldRemoteApi import ESMFoldRemoteApi
from Bio.PDB.Atom import Atom
from evodesign.Chain import Chain
import numpy as np
import os

sequence = 'GREETINGSTRAVELER'
filename = 'tests/esmfold_prediction.pdb'
predictor = ESMFoldRemoteApi()
backbone, plddt = predictor(sequence, filename)
test_passed = type(plddt) == np.float64 and \
              0 <= plddt and plddt <= 1. and \
              type(backbone) == list and \
              len(backbone) == 4 * len(sequence) and \
              sum([ isinstance(atom, Atom) for atom in backbone ]) == len(backbone) and \
              sum([ atom.get_name() in Chain.BACKBONE_ATOMS for atom in backbone ]) == len(backbone) and \
              os.path.isfile(filename)
print('PASSED' if test_passed else 'FAILED')
os.remove(filename)
