#!/usr/bin/env python2
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys

suppl = Chem.SDMolSupplier('../cdk2.sdf')
ms = [suppl.next(), suppl.next()]
fps = [MACCSkeys.GenMACCSKeys(x) for x in ms]
print DataStructs.FingerprintSimilarity(fps[0], fps[1])
