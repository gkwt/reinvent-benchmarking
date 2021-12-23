import numpy as np
import copy

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

class EarlyStopping():
    ''' Class that checks criteria for early stopping. Saves the best model weights.
    '''
    def __init__(self, patience, min_delta, mode='minimize'):
        self.patience = patience
        self.best_weights = None
        self.checkpoint = 0
        self.best_epoch = 0
        if mode == 'maximize':
            self.monitor_fn = lambda a, b: np.greater(a - min_delta, b)
            self.best_val = -np.inf
        elif mode == 'minimize':
            self.monitor_fn = lambda a, b: np.less(a + min_delta, b)
            self.best_val = np.inf
        else:
            raise ValueError(f'Mode should be either minimize or maximize.')

    def check_criteria(self, net, epoch, new_val):
        ''' Compare with value in memory. If there is an improvement, reset the checkpoint and
        save the model weights.
        Return True if stopping criteria is met (checkpoint is exceeds patience), otherwise, return False.
        '''
        if self.monitor_fn(new_val, self.best_val):
            self.best_val = new_val
            self.checkpoint = 0
            self.best_weights = copy.deepcopy(net.state_dict())
            self.best_epoch = epoch
        else:
            self.checkpoint += 1
        
        return self.checkpoint > self.patience

    def restore_best(self, net, verbose=True):
        if verbose: print(f'        Early stopping at epoch: {self.best_epoch}       loss: {self.best_val}')
        net.load_state_dict(self.best_weights)
        return net