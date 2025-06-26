import streamlit as st
from utils.molecule_utils import smiles_to_mol_image, calculate_mol_weight, get_rdk_mol_from_sml, get_atomic_contributions_to_logp 

st.set_page_config(page_title=":rainbow[Chemiformatics Assistant To Experimental Chemist]", layout="centered")
st.write(
    """
    :rainbow[Enrich your research with quick cheminformatics analysis. e.g.\n View atomic contributions to logP values in your molecule.]
    """
)
st.title("CrippenLogP 2D Map Tool")
st.markdown("""
    <style>
    .main-title {
        font-size: 2.5em;
        color: #2E86C1;
        text-align: center;
        font-weight: bold;
        margin-bottom: 20px;
    }
    .section-header {
        font-size: 1.8em;
        color: #1B4F72;
        margin-top: 20px;
        font-weight: bold;
    }
    .content-text {
        font-size: 1.1em;
        color: #34495E;
        line-height: 1.6;
        margin-bottom: 15px;
    }
    .link {
        color: #2874A6;
        text-decoration: none;
        font-weight: bold;
    }
    .link:hover {
        text-decoration: underline;
    }
    .input-instruction {
        font-size: 1.2em;
        color: #17202A;
        font-style: italic;
        text-align: center;
        margin-top: 20px;
    }
    </style>
""", unsafe_allow_html=True)

st.markdown("""
    <div class="main-title">Explore Molecular Lipophilicity with Our 2D CrippenLogP Map Tool</div>

    <div class="section-header">What It Does</div>
    <div class="content-text">
        Generates a static 2D map visualizing atomic contributions to a molecule's CrippenLogP value, 
        highlighting how each atom affects lipophilicity.
    </div>

    <div class="section-header">Scientific Basis</div>
    <div class="content-text">
        Based on Crippen's LogP method 
        (<a href="https://doi.org/10.1021/ci990307l" class="link">Wildman & Crippen, 1999</a>), 
        this tool calculates hydrophobicity, a critical factor in drug solubility, membrane permeability, 
        and bioavailability.
    </div>

    <div class="section-header">Why It Matters</div>
    <div class="content-text">
        Mapping LogP contributions helps optimize molecular structures for improved ADME 
        (absorption, distribution, metabolism, excretion) properties in drug development 
        (<a href="https://doi.org/10.1016/S1359-6446(03)02671-2" class="link">Lipinski, 2004</a>).
    </div>

    <div class="section-header">How to Use</div>
    <div class="content-text">
        Input a valid SMILES string of your molecule in the input field below and press "Enter." 
        Need a SMILES string? Use an online structure editor like 
        <a href="https://www.chemdoodle.com/" class="link">ChemDoodle</a> to draw your molecule 
        and generate its SMILES.
    </div>

    <div class="section-header">Analyzing Results</div>
    <div class="content-text">
        The static map color-codes atoms by LogP contribution (positive: hydrophobic, negative: hydrophilic). 
        Review the legend for value ranges. Identify key structural features influencing lipophilicity 
        to guide molecular design.
    </div>

    <div class="input-instruction">Start mapping your molecule’s potential today!</div>
""", unsafe_allow_html=True)

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

        
