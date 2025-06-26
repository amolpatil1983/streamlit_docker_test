import streamlit as st
from utils.molecule_utils import smiles_to_mol_image, calculate_mol_weight, get_rdk_mol_from_sml, get_atomic_contributions_to_logp 

st.title("RDKit + Sklearn Streamlit App")

# User input for SMILES string
smiles = st.text_input("Enter SMILES string", "CC(=O)OC1=CC=CC=C1C(=O)O")

if smiles:
    # Generate molecule image
    img_path = smiles_to_mol_image(smiles)
    if img_path:
        st.image(img_path, caption="Molecule Visualization")
    
    # Calculate molecular weight
    mol_weight = calculate_mol_weight(smiles)
    if mol_weight:
        st.write(f"Molecular Weight: {mol_weight:.2f} g/mol")
    
    rdkmolobj = get_rdk_mol_from_sml(smiles)
    image = get_atomic_contributions_to_logp(rdkmolobj)
            if image:
                st.image(image, caption="Atomic Contributions to LogP", use_column_width=True)
            else:
                st.warning("Failed to generate visualization. Please check the SMILES string.")

        
