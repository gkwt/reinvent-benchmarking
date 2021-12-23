import rdkit.Chem as Chem
from rdkit.Chem import Draw, Descriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

class fitness_function():

    def __init__(self):
        pass

    def __call__(self, smi):
        try:
            mol = Chem.MolFromSmiles(smi)
            log_P = Descriptors.MolLogP(mol)
            return log_P
        except:
            return -100.0
