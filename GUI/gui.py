#Packages
pip install git+https://github.com/NissrineH/Docking_Interface/blob/main/environment.yml
import streamlit as st
from pymol import cmd
import py3Dmol
import prolif as plf
import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from openbabel import pybel
from vina import Vina
from meeko import MoleculePreparation
from meeko import obutils
import prolif as plf
from prolif.plotting.network import LigNetwork
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser
import sys, os, random
sys.path.insert(1, '/home/nissrine/Jupyter_Dock/GUI/utilities/')
from utils import fix_protein, getbox, generate_ledock_file, pdbqt_to_sdf, dok_to_sdf
import warnings
warnings.filterwarnings("ignore")
import subprocess
import pandas as pd

#Working folder
os.chdir('/home/nissrine/Jupyter_Dock/GUI/test/')

bashCommand = "mkdir -p ledock_inputfiles/"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
bashCommand = "mkdir -p ledock_outfiles/"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)



#Classes
class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:          
            return 0
            
class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


#################################################################################################################################

#Functions

#Get pdb format of the receptor 
def load_receptor(receptor):
    pdb_lists = [receptor] 
    for x in pdb_lists:
        cmd.fetch(x, type='pdb')
        cmd.select(x)
        cmd.delete(x)
        
# Get chain of interest from the receptor 
def select_chains(ch):
    structure = PDBParser().get_structure(receptor, receptor+'.pdb')    
    for chain in ch:
        pdb_chain_file = '4giz_chain_{}.pdb'.format(ch)                                 
        io_w_no_h = PDBIO()               
        io_w_no_h.set_structure(structure)
        io_w_no_h.save('{}'.format(pdb_chain_file), ChainSelect(chain))

# Remove water molecules and ions
def remove_HOH(ch):
    pdb = PDBParser().get_structure('4giz_chain_{}'.format(ch), '4giz_chain_{}.pdb'.format(ch))
    io = PDBIO()
    io.set_structure(pdb)
    io.save('4giz_chain_{}_non_het.pdb'.format(ch, NonHetSelect()))

# Receptor preparation: addH and charges
def receptor_preparation(ch):
    bashCommand = "../bin/lepro_linux_x86 {'4giz_chain_{}_non_het.pdb'.format(ch)}"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.rename('pro.pdb','4giz_chain_{}_non_het_H.pdb'.format(ch))
              
#Get ligands
def load_ligands(ligands, smiles):
    for ligand in ligands:
        name=ligand.name
        lines = open(name).read().splitlines()
        for line in range(int(round((len(lines)*ratio)/100))):
            lines[line]=lines[line].split(' ')[0]
            smiles.append(lines[line])
    return(smiles)

# Ligands preparation: addH and charges
def ligands_preparation(smiles):
    for index,smi in enumerate(smiles):
        mol=pybel.readstring(string=smi,format='smiles')
        mol.title='mol_'+str(index)
        mol.make3D('mmff94s')
        mol.localopt(forcefield='mmff94s', steps=500)
        out=pybel.Outputfile(filename='ledock_inputfiles/'+'mol_'+str(index)+'.mol2',format='mol2',overwrite=True)
        out.write(mol)
        out.close()
    
# Docking box definition   
def box(ch):    
    cmd.load(filename='4giz_chain_{}_non_het.pdb'.format(ch),format='pdb',object='prot') 
    X,Y,Z=getbox(selection='prot',extending=6.0,software='ledock')
    cmd.delete('all')
    #print(X,'\n',Y,'\n',Z)
    return (X, Y, Z)

# Molecular Docking
def ledock(ch, X, Y, Z):
    l_list=[]
    for file in os.listdir('ledock_inputfiles/'):
        if 'mol2' in file:
            l_list.append('ledock_inputfiles/'+file+'\n')
            
    generate_ledock_file(receptor='4giz_chain_{}_non_het.pdb'.format(ch),l_list=l_list,
                         l_list_outfile='ligand.list',
                         x=[X['minX'],X['maxX']],
                         y=[Y['minY'],Y['maxY']],
                         z=[Z['minZ'],Z['maxZ']],
                         n_poses=10,
                         rmsd=1.0,
                         out='dock.in')

    
def dock():  
    bashCommand = "../bin/ledock_linux_x86 dock.in"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print(output, error)

#DOK results file conversion to SDF
def results_file():
    for file in os.listdir('ledock_inputfiles/'):
        if '.dok' in file:
            os.rename('ledock_inputfiles/'+file,'ledock_outfiles/'+file)
            dok_to_sdf(dok_file='ledock_outfiles/'+file,output='ledock_outfiles/'+file.replace('dok','sdf'))




def docking(receptor, ligands, ratio):
    ch='C'
    smiles=[]
    load_receptor(receptor)
    select_chains(ch)
    remove_HOH(ch)
    rec_clean=receptor_preparation(ch)
    sm=load_ligands(ligands, smiles)
    ligands_preparation(sm)
    X, Y, Z=box(ch)
    ledock(ch, X, Y, Z)
    dock()
    results_file()

    #return (error)

#################################################################################################################################
#Interface
col1, col2, col3 = st.columns(3)
with col2:
    st.subheader('Molecular Docking')
receptor = st.text_input('ID receptor')
ligands = st.file_uploader('Choose ligands path', accept_multiple_files=True)
for ligand in ligands:
    bytes_data = ligand.read()
    #st.write("filename:", ligand.name)
    #st.write(bytes_data)
ratio = st.slider('Choose a ratio', min_value=0.00, max_value=100.00,value=0.01, step=0.001)
st.write("The selected ratio is ", ratio, '%')
col1, col2, col3 = st.columns(3)
with col2:
    if st.button('Get Docking results'):
        result = docking(receptor, ligands, ratio)
        for file in os.listdir('ledock_outfiles/'):
            if 'sdf' in file:
                pose=Chem.SDMolSupplier('ledock_outfiles/'+file,False)[0]
                p=Chem.MolToMolBlock(pose)
                #print('Name: {} | Pose: {} | Score: {}'.format(file.split('.')[0],pose.GetProp('Pose'),pose.GetProp('Score')))

                st.write('{} ------- Score: {}'.format(file.split('.')[0],pose.GetProp('Score')))



