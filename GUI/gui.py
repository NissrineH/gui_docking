#Packages
import streamlit as st

def load_ligands(ligands, smiles):
    for ligand in ligands:
        name=ligand.name
        lines = open(name).read().splitlines()
        for line in range(int(round((len(lines)*ratio)/100))):
            lines[line]=lines[line].split(' ')[0]
            smiles.append(lines[line])
    return(smiles)
        
def docking(receptor, ligands, ratio):
    ch='C'
    smiles=[]
    load_ligands(ligands, smiles)


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


