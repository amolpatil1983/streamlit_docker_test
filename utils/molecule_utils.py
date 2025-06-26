from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem import rdEHTTools
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import SimilarityMaps
from sklearn.linear_model import LinearRegression
import numpy as np
import os

def smiles_to_mol_image(smiles, output_path="mol_image.png"):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        Draw.MolToFile(mol, output_path, size=(300, 300))
        return output_path
    except:
        return None

def calculate_mol_weight(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Descriptors.MolWt(mol)
    except:
        return None

def predict_solubility(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Calculate simple descriptors for solubility prediction
        features = [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol)
        ]
        # Dummy sklearn model (replace with trained model for production)
        model = LinearRegression()
        # Mock training data (for demo purposes)
        X_train = np.array([[100, 1.0, 1, 2], [200, 2.0, 0, 3], [300, 0.5, 2, 1]])
        y_train = np.array([-2.0, -1.5, -3.0])
        model.fit(X_train, y_train)
        # Predict solubility
        return model.predict([features])[0]
    except:
        return None

# MolFromSmiles
def get_rdk_mol_from_sml(smiles):
    try:
        rdk_mol = Chem.MolFromSmiles(smiles)
        return rdk_mol
    except:
        return None
        
# Function to show atomic contributions to LogP values
def get_atomic_contributions_to_logp(rdk_mol):
    try:
        # Calculate Crippen (logP) contributions for each atom
        at_contribs = rdMolDescriptors._CalcCrippenContribs(rdk_mol)
        # Create a 2D drawing canvas
        d = Draw.MolDraw2DCairo(400, 400)
        # Generate similarity map from logP contributions
        SimilarityMaps.GetSimilarityMapFromWeights(rdk_mol, [x[0] for x in at_contribs], draw2d=d)
        d.FinishDrawing()
        # Convert drawing to PNG bytes
        png_data = d.GetDrawingText()
        # Convert PNG bytes to an image for Streamlit
        image = Image.open(io.BytesIO(png_data))
        return image
    except:
        return None
        
def molsim(smiles):
    try:
        mh = Chem.AddHs(smiles)
        rdDistGeom.EmbedMolecule(mh)
        _,res = rdEHTTools.RunMol(mh)
        static_chgs = res.GetAtomicCharges()[:atorvastatin.GetNumAtoms()]
        d = Draw.MolDraw2DCairo(400, 400)
        SimilarityMaps.GetSimilarityMapFromWeights(atorvastatin,list(static_chgs),draw2d=d)
        d.FinishDrawing()
        # show_png(d.GetDrawingText())
        return d
    except:
        return None
