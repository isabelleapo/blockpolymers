from contextlib import contextmanager
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess as sp

def factors(x):
    factorlist = []
    for i in range(1,x+1):
        if x%i == 0:
            factorlist.append(i)
    return factorlist


def generate_xyz(mol):
    print(Chem.MolToMolBlock(mol), file=open('mol','w+'))

    p = sp.Popen(['babel', 'mol', 'xyz'], stdout=sp.PIPE)
    o, e = p.communicate()

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def run_calc(calc_params):
    p = sp.Popen(calc_params, stdout=sp.PIPE, encoding='utf8')
    output, _ = p.communicate()
    return output
