import streamlit as st
from utils.molecule_utils import smiles_to_mol_image, calculate_mol_weight, get_rdk_mol_from_sml, get_atomic_contributions_to_logp 

st.title(":rainbow[Chemiformatics Assistant To Experimental Chemist]")
st.write(
    """
    Enrich your research with quick cheminformatics analysis.
    e.g. View atomic contributions to logP values in your molecule.
    """
)

# User input for SMILES string
smiles = st.text_input("Enter SMILES string, e.g O=C(O)C[C@H](O)C[C@H](O)CCn2c(c(c(c2c1ccc(F)cc1)c3ccccc3)C(=O)Nc4ccccc4)C(C)C for atorvastatin")

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
        st.image(image, caption="Atomic Contributions to LogP") #, use_container_width=True)
    #else:
    #    st.warning("Failed to generate visualization. Please check the SMILES string.")

        
