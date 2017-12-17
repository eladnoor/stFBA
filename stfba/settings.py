# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 16:36:02 2016

@author: noore
"""

import os
import pandas as pd
import inspect
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
BASE_DIR = os.path.join(*os.path.split(SCRIPT_DIR)[0:-1])
DATA_DIR = os.path.join(BASE_DIR, 'data')
MODEL_DIR = os.path.join(BASE_DIR, 'model')
RESULT_DIR = os.path.join(BASE_DIR, 'res')
IJO1366_JSON_FNAME = os.path.join(MODEL_DIR, 'iJO1366.json')
CORE_SBML_FNAME = os.path.join(MODEL_DIR, 'e_coli_core.xml.gz')

M = 1e6
eps = 1e-6

R = 8.31e-3   # kJ/(K*mol)
DEFAULT_TEMP = 298.15  # K
DEFAULT_IONIC_STRENGTH = 0.1  # mM
DEFAULT_PH = 7.0
DEFAULT_PMG = 14.0
DEFAULT_PHASE = 'aqueous'
RT = R * DEFAULT_TEMP
RTlog10 = RT * np.log(10)


def get_stoichiometry_from_model(model):
    sparse = []
    for rxn in model.reactions:
        for met, coeff in rxn.metabolites.items():
            sparse.append([rxn.id, met.id, coeff])

    sparse = pd.DataFrame(sparse, columns=['bigg.reaction',
                                           'bigg.metabolite', 'stoichiometry'])
    S = sparse.pivot(index='bigg.metabolite', columns='bigg.reaction',
                     values='stoichiometry')
    S.fillna(0, inplace=True)
    return S

