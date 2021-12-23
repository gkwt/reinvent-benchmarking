import numpy as np

import rdkit.Chem as Chem
from rdkit.Chem import Draw, Descriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

def fitness_function(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        log_P = Descriptors.MolLogP(mol)
        return log_P
    except:
        return None

class custom_score():

    def __init__(self):
        pass

    def __call__(self, smi):
        return fitness_function(smi)
