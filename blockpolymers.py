import stk
import rdkit, rdkit.Chem as rdkit
import subprocess as sp
import itertools
from .utilities import *
from joblib import Parallel, delayed
import re
import os

class BlockPolymers:
    """
    Parameters
    ----------
    A_smiles_path: :class:'str'
        Path to text file containing list of smiles for monomer A

    B_smiles_path: :class:'str'
        Path to text file containing list of smiles for monomer B

    length: :class:'int'
        Number of monomer units in polymer

    nconfs: :class:'int'
        Number of conformers to embed within conformer search

    solvent: :class:'str'
        Solvent to be used in xTB (available solvents: toluene,thf,
        methanol, h2o, ether, chcl3, acetonitrile, acetone, cs2)

    min_osc_strength:  :class:'float'
        Minimum oscillator strength of the excitation selected to be
        the optical gap

    Methods
    -------
    screen:
        Builds binary co-polymers out of all possible combinations of monomers
        in library A and library B. Calculates IP, EA, optical gap
        and oscillator strength

    conformer_search:
        Finds lowest energy conformers, returns the .mol files and saves them
        into 'lowest_conformers' folder

    Returns
    -------
    .mol file containing lowest energy conformer
    Text file containing calculated properties


    """
    def __init__(self,
                 A_smiles_path,
                 B_smiles_path,
                 length = 12,
                 nconfs = 30,
                 solvent = None,
                 min_osc_strength = 0
                 ):

        self.A_smiles_path = A_smiles_path
        self.B_smiles_path = B_smiles_path
        self.length = length
        self.nconfs = nconfs
        self.solvent = solvent
        self.min_osc_strength = min_osc_strength

    def _create_smiles_list(self):
        with open(self.A_smiles_path, 'r') as f:
            f2 = f.read()
            A_list = f2.split()

        with open(self.B_smiles_path, 'r') as f3:
            f4 = f3.read()
            B_list = f4.split()

        return A_list, B_list
        print(A_list)

    def get_sequences(self):
        sets = int(self.length/2)
        setlist = []
        for n in factors(sets):
            y, sett = n, int(sets/n)
            pattern = ['A']*y + ['B']*y
            pattern2 = ''.join(tuple(pattern))
            setlist.append(pattern2)
        return setlist

    def build_polymers(self, monomer_pair):
        a, b = monomer_pair
        A_id, A_smiles = a
        B_id, B_smiles = b
        Amol = rdkit.MolFromSmiles(A_smiles)
        Bmol = rdkit.MolFromSmiles(B_smiles)
        A = stk.StructUnit2.rdkit_init(Amol, 'bromine')
        B = stk.StructUnit2.rdkit_init(Bmol, 'bromine')
        polymers = []
        try:
            for x in self.get_sequences():
                print(self.get_sequences())
                seqstring = ''
                for l  in x:
                    if l == 'A':
                        seqstring += str(0)
                    if l == 'B':
                        seqstring += str(1)
                id = f'{A_id}_{B_id}_{seqstring}'
                N = int(self.length/len(x))
                polymer = stk.Polymer([A,B],
                          stk.Linear(x, [0]*self.length, n=N))

                polymer_smiles = rdkit.MolToSmiles(
                                 rdkit.RemoveHs(polymer.mol), canonical= True)
                pol_id = (id,A_smiles,B_smiles,polymer_smiles, polymer)
                polymers.append(pol_id)
        except Exception as e:
            print(f'Error in building polymer: {e}')

        return polymers

    def conformer_search(self, polymers):
        for poly in list(polymers):
            id, _, _, _, polymer = poly
            polymer = polymer.mol
            confs = rdkit.AllChem.EmbedMultipleConfs(
                              polymer, self.nconfs, rdkit.AllChem.ETKDG())
            rdkit.SanitizeMol(polymer)

            lowest_energy = 10**10
            for conf in confs:
                ff = rdkit.AllChem.MMFFGetMoleculeForceField(
                polymer, rdkit.AllChem.MMFFGetMoleculeProperties(polymer), confId=conf)
                ff.Initialize()
                energy = ff.CalcEnergy()

                if energy < lowest_energy:
                    lowest_energy = energy
                    lowest_conf = conf

            rdkit.MolToMolFile(polymer, f'lowest_conformers/{id}_lowest-conformer.mol',
                           confId=lowest_conf)


    def _create_prop_files(self):
        try:
            os.mkdir('lowest_conformers')
        except Exception as e:
            print(e)

        with open('properties.txt', 'w') as f:
            f.write(('ID\txTB_IP\txTB_EA\txTB_Opticalgap\txTB_oscillator_strength\tSmiles\tA\tB\n'))

    def properties(self, polymers):
        for poly in list(polymers):
            id, A_smiles, B_smiles, polymer_smiles, polymer = poly
            with open(f'lowest_conformers/{id}_lowest-conformer.mol', 'r') as f:
                molfile = f.read()
            molfile = rdkit.AddHs(rdkit.MolFromMolBlock(molfile))
            rdkit.AllChem.EmbedMolecule(molfile, rdkit.AllChem.ETKDG())
            dirname = str(id)+'_properties'
            try:
                os.mkdir(dirname)
            except Exception as e:
                pass
            with cd(dirname):
                os.system('pwd')
                xyz = generate_xyz(molfile)

                if self.solvent != None:
                    output = run_calc(['xtb','xyz','-opt','-gbsa', self.solvent])

                else:
                    output = run_calc(['xtb','xyz','-opt','-gbsa', self.solvent])

                if self.solvent != None:
                    output = run_calc(['xtb','xtbopt.xyz', '-vip', '-gbsa', self.solvent])

                else:
                    output = run_calc(['xtb','xtbopt.xyz', '-vip', '-gbsa'])


                with open('temp.txt', 'w', encoding='utf8') as f:
                    f.write(output)

                with open('temp.txt', 'r', encoding='utf8') as f:
                    temp = f.read()
                    temp = str(temp)

                pattern = re.compile('(?<=delta SCC IP [(]eV[)].).*\d\.\d{4}')
                ip = pattern.findall(temp)
                print(ip)

                rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
                o, e = rt.communicate()

                if self.solvent != None:
                    output = run_calc(['xtb','xtbopt.xyz', '-vea', '-gbsa', self.solvent])

                else:
                    output = run_calc(['xtb','xtbopt.xyz', '-vea', '-gbsa'])

                with open('temp.txt', 'w', encoding='utf8' ) as f:
                    f.write(output)

                with open('temp.txt', 'r', encoding='utf8') as f:
                    temp = f.read()
                    temp = str(temp)

                pattern = re.compile('(?<=delta SCC EA [(]eV[)].).*\d\.\d{4}')
                ea = pattern.findall(temp)

                rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
                o, e = rt.communicate()

                output = run_calc(['xtb','xtbopt.xyz'])

                output = run_calc(['stda','-xtb', '-e', '8'])

                with open('temp.txt', 'w', encoding='utf8') as f:
                    f.write(output)

                with open('temp.txt', 'r', encoding='utf8') as f:
                    temp = f.read()
                    temp = str(temp)

                pattern = re.compile(r'\n\s+\d+\s+(\d\.\d+)\s+[\d.-]+\s+(\d+\.\d+)[>\(\)\s\d.-]')
                og = pattern.findall(temp)

                rt = sp.Popen(['rm','temp.txt'], stdout=sp.PIPE)
                o, e = rt.communicate()

                polymer_smiles = rdkit.MolToSmiles(rdkit.RemoveHs(molfile), canonical=True)

                string = str(id)+'\t'
                for match in ip:
                    match = match.strip()
                    string += match+'\t'
                for match in ea:
                    match = match.strip()
                    string += match+'\t'
                opgap = []
                osc = []
                for match in og:
                    if float(match[1]) > self.min_osc_strength:
                        opgap.append(match[0])
                        osc.append(match[1])
                        break
                string += opgap[0]+'\t'
                string += osc[0]+'\t'
                string += polymer_smiles+'\t'
                string += A_smiles+'\t'
                string += B_smiles+'\n'

                with open('props.txt', 'w') as f:
                    f.write(string)
            os.system('pwd')




    def _write_to_props(self):
        dirs = str(os.listdir('.'))
        pattern = re.compile('\d+_\d+_\d+_properties')
        match = pattern.findall(dirs)
        for m in match:
            try:
                with open(m+'/props.txt', 'r') as p:
                    prop = p.read()
                with open('properties.txt', 'a') as f:
                    f.write(prop)
            except Exception as e:
                print(e, 'Error in writing to file')

    def _poly_screen(self, item):
        try:
            polymers = self.build_polymers(item)
            self.conformer_search(polymers)
        except Exception as e:
            print(e, 'Error in polymer formation or conformer search')
        try:
            self.properties(polymers)
        except Exception as e:
            print(e, 'Error in property calculation')
            os.chdir('../')

    def screen(self, nprocs=1):
        self._create_prop_files()
        A_list, B_list = self._create_smiles_list()
        results = Parallel(n_jobs=nprocs)(delayed(self._poly_screen)(item) for item in itertools.product(enumerate(A_list), enumerate(B_list)))
        self._write_to_props()
