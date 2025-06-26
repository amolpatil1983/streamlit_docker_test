import streamlit as st
from utils.molecule_utils import smiles_to_mol_image, calculate_mol_weight, get_rdk_mol_from_sml, get_atomic_contributions_to_logp

st.set_page_config(layout="wide")
st.title(":rainbow[Cheminformatics Assistant]")

st.info("""
### 🧪 What It Does  
Generates a **static 2D map** visualizing atomic contributions to a molecule's **Crippen LogP** value, highlighting how each atom affects **lipophilicity**.
---
### 💡 Why It Matters  
Mapping LogP contributions helps optimize molecular structures for improved **ADME** properties:  
*Absorption, Distribution, Metabolism, Excretion* – crucial in drug development.  
📖 [Lipinski, 2004](https://www.sciencedirect.com/science/article/abs/pii/S0169409X00001290?via%3Dihub)
---
### 🔬 Scientific Basis  
Based on **Crippen's LogP method**  
📖 [Wildman & Crippen, 1999](https://doi.org/10.1021/ci990307l)  
Calculates **hydrophobicity**, a key factor in drug **solubility**, **membrane permeability**, and **bioavailability**.
---
### 🧬 How to Use  
1. Input a **valid SMILES string** of your molecule in the field below.  
2. Press **Enter** to generate the visualization.  
🔧 Need a SMILES string? Try the [PubChem Sketcher](https://pubchem.ncbi.nlm.nih.gov/#input=draw&draw=true) to draw and convert your molecule.
---
### 🎯 Analyzing Results  
- **Green atoms**: Contribute to **hydrophobicity**  
- **Red atoms**: Contribute to **hydrophilicity**  
Use this insight to identify **key features** influencing lipophilicity and guide your **molecular design**.
---
### 🚀 Start mapping your molecule’s potential today!
""")

# User input for SMILES string
smiles = st.text_input("Enter SMILES string", "O=C(O)C[C@H](O)C[C@H](O)CCn2c(c(c(c2c1ccc(F)cc1)c3ccccc3)C(=O)Nc4ccccc4)C(C)C")

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

        
