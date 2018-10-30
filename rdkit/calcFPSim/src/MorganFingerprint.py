#!/usr/bin/env python2
from rdkit import Chem
from rdkit import DataStructs
#from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem

m1 = Chem.MolFromSmiles('Cc1ccccc1')
fp1 = AllChem.GetMorganFingerprint(m1,2)
#fp1 = AllChem.GetMorganFingerprintAsBitVect(m1,2,nBits=1024)
m2 = Chem.MolFromSmiles('Cc1ncccc1')
fp2 = AllChem.GetMorganFingerprint(m2,2)
#fp2 = AllChem.GetMorganFingerprintAsBitVect(m2,2,nBits=1024)
#print DataStructs.DiceSimilarity(fp1,fp2)
print DataStructs.TanimotoSimilarity(fp1,fp2)
