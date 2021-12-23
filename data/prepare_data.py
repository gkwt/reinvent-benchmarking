import pandas as pd

data_path = 'hce.txt'
sep = ' '
header = None
smile_name = 0

# dataset load
data = pd.read_csv(data_path, sep=sep, header=header)
smiles = data[smile_name]

# create smile file
with open('data.smi', 'w') as f:
    for smi in smiles:
        f.write(smi+'\n')
