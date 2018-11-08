#!/usr/bin/env python2
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

m = Chem.SDMolSupplier("../zinc_537755.sdf").next()
print m.GetNumAtoms()
print rdMolDescriptors.CalcNumHBD(m)
print rdMolDescriptors.CalcNumHBA(m)
print rdMolDescriptors.CalcNumRotatableBonds(m)
print rdMolDescriptors.CalcExactMolWt(m)
print rdMolDescriptors.CalcTPSA(m)
print rdMolDescriptors.CalcCrippenDescriptors(m)[0]
