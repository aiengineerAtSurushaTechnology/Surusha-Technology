# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 13:18:36 2022

@author: Acer
"""

#%%
import streamlit as st
from stmol import showmol
import py3Dmol

from rdkit import Chem
from rdkit.Chem import AllChem

#st.title("""RDkit + Py3MOL """)

def run_search():
    def makeblock(smi):
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mblock = Chem.MolToMolBlock(mol)
        return mblock
    
    def render_mol(xyz):
        xyzview = py3Dmol.view()
        xyzview.addModel(xyz,'mol')
        xyzview.setStyle({'stick':{}})
        xyzview.setBackgroundColor('white')
        xyzview.zoomTo()
        showmol(xyzview,height=500, width=700)
        
    compound_smiles = st.text_input('SMILES please', 'CC')
    blk = makeblock(compound_smiles)
    render_mol(blk)