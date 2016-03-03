#! /usr/bin/env python
#
# Copyright (C) 2016 Jing Lu <ajingnk@gmail.com>
# License: Apache

# -*- coding: utf-8 -*-

# pylint: disable=too-few-public-methods
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs

from .model import SMILE_COLUMNNAME, POTENCY


class FingerPrintType:
    """
    representing fingerprint types
    """
    def __init__(self, name, fp_func, metadata):
        """
        Initialize the fingerprint type

        :param name: name of fingerprint
        :param fp_func: the function to generate a fingerprint
        :param meta: string for fingerprint document
        :return:
        """
        self.name = name
        self.fp_func = fp_func
        self.metadata = metadata

    def to_dict(self):
        """
        Show the information for this fingerprint

        :return: dictionary with basic info
        """
        return {"name": self.name, "metadata": self.metadata}


class PropertyType:
    """
    representing biological or chemical properties
    """
    def __init__(self, name, metadata, transfunc=None, colname=None):
        """
        :param name: a string name of property
        :param metadata: a string with property meaning
        :param transfunc: a function to generate the property
        :return:
        """
        self.name = name
        self.metadata = metadata
        self.transfunc = transfunc  # the signiture for transfunc is transfunc(a_value, a_mol_dict)
        self.colname = colname

    def set_col_name(self, col_name):
        """
        Set the property name from the input file

        :param col_name: original column name in the input file
        :return:
        """
        self.colname = col_name

    def gen_property(self, mol_dict = None):
        """
        generate value for this property type

        :param prop_name: the name of the property
        :param mol_dict: other information about the molecule
        :return: a generated value for this property
        """
        # if self.colname is None:
        #     raise Exception("Please set the column name for this property")

        if self.colname in mol_dict:
            return mol_dict[self.colname]

        if self.name in mol_dict:
            return mol_dict[self.name]

        if not self.transfunc is None and not mol_dict is None:
            return self.transfunc(mol_dict)
        else:
            raise Exception("please provide the transformation function and molecule information for " + str(self))


    def to_dict(self):
        """
        Show the information

        :return: dictionary with basic info
        """
        return {"name": self.name, "metadata": self.metadata}

    def __str__(self):
        return self.name

# property functions
def _lig_eff(mol_dict):
    smile = mol_dict[SMILE_COLUMNNAME]
    m = Chem.MolFromSmiles(smile)
    num_heavy = m.GetNumHeavyAtoms()
    if POTENCY in mol_dict:
        ic50  = mol_dict[POTENCY]
        return round(1.37 * (9 - math.log10(ic50)) / num_heavy, 5)
    elif "pIC50" in mol_dict:
        pic50 = mol_dict["pIC50"]
        return round(1.37 * pic50 / num_heavy, 5)
    else:
        raise Exception("Cannot calculate ligand efficiency, please change your input file column name to IC50 or pIC50.")


def _slogp(mol_dict):
    smile = mol_dict[SMILE_COLUMNNAME]
    m = Chem.MolFromSmiles(smile)
    return round(Chem.rdMolDescriptors.CalcCrippenDescriptors(m)[0], 5)

def _pic50(mol_dict):
    ic50 = mol_dict[POTENCY]
    return round(9 - math.log10(float(ic50)), 5)

ecfp6 = FingerPrintType(name = "ECFP6", fp_func= lambda mol: AllChem.GetMorganFingerprint(mol, 3), metadata = "Extended Connectivity fingerprint, implemented in <a href=\"http://www.rdkit.org\">RDKit</a>. <br/>Parameters used: Radius = 3")

atom_pair = FingerPrintType(name = "AtomPair", fp_func= Pairs.GetAtomPairFingerprint, metadata = "Atom Pairs as Molecular Features, describe in  R.E. Carhart, D.H. Smith, R. Venkataraghavan; \"Atom Pairs as Molecular Features in Structure-Activity Studies: Definition and Applications\" JCICS 25, 64-73 (1985).implemented in <a href=\"http://www.rdkit.org\">RDKit</a>. <br/>")

lig_eff = PropertyType(name = "Lig_Eff", metadata = "Ligand efficiency. The value is calculated by the function 1.37 * pIC50 / a_heavy", transfunc = _lig_eff)

slogp   = PropertyType(name = "SLogP", metadata = "SLogP, the coefficients are a measure of the difference in solubility of the compound in water and octanol. describe in    S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999) R.E. Carhart, D.H. Smith, R. Venkataraghavan; \"Atom Pairs as Molecular Features in Structure-Activity Studies:", transfunc = _slogp)

ic50   = PropertyType(name = "IC50", colname = "IC50", metadata = "The half maximal inhibitory concentration (IC50) is a measure of the effectiveness of a substance in inhibiting a specific biological or biochemical function.")

pic50 = PropertyType(name = "pIC50", metadata = "This number assumes IC50 in nM unit, so it is calculated by 9 - log(IC50). Please change your data or the code to make it appropriate.", transfunc = _pic50)

bindingdb = {"name": "BindingDB", "link": "https://www.bindingdb.org/bind/chemsearch/marvin/MolStructure.jsp?monomerid="}
chebi = {"name": "CHEBI", "link": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId="}
pubchem = {"name": "PubChem", "link": "https://pubchem.ncbi.nlm.nih.gov/compound/"}

DEFAULT_FINGERPRINT_TYPES = [ecfp6, atom_pair]
DEFAULT_ACTIVITY_TYPES = [ic50, pic50]
DEFAULT_PROPERTY_TYPES = [lig_eff, slogp]
DEFAULT_EXTERNAL = [bindingdb, pubchem]