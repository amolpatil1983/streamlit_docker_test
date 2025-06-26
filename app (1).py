import streamlit as st
from utils.molecule_utils import smiles_to_mol_image, calculate_mol_weight, predict_solubility

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
    
    # Predict solubility using sklearn model
    solubility = predict_solubility(smiles)
    if solubility:
        st.write(f"Predicted Solubility (logS): {solubility:.2f}")