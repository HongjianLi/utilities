#!/usr/bin/env python2
import sys
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdGeometry

suppl = Chem.SDMolSupplier(sys.argv[1]) # Path to a input SDF file.
mols = [mol for mol in suppl if mol is not None]
#for mol in mols:
#	AllChem.Compute2DCoords(mol)
#	Draw.MolToFile(mol, mol.GetProp("_Name") + "." + sys.argv[2]) # Supported file formats are pdf, svg, ps, png.
mcs = rdFMCS.FindMCS(mols)
ref = Chem.MolFromSmarts(mcs.smartsString)
AllChem.Compute2DCoords(ref)
refCnf = ref.GetConformer()
molsPerRow = 3
subImgSize = (400, 400)
nRows = 1+(len(mols)-1)//molsPerRow
img = Image.new("RGBA", (molsPerRow * subImgSize[0], nRows * subImgSize[1]), (255,255,255,0))
for k, mol in enumerate(mols):
#	AllChem.GenerateDepictionMatching2DStructure(mol, ref)
	coordMap = {}
	for i, idx in enumerate(mol.GetSubstructMatch(ref)):
		pt3 = refCnf.GetAtomPosition(i)
		coordMap[idx] = rdGeometry.Point2D(pt3.x, pt3.y) # pt3.z == 0.0
	AllChem.Compute2DCoords(mol, clearConfs=True, canonOrient=False, coordMap=coordMap)
	row = k//molsPerRow
	col = k-molsPerRow*row
	img.paste(Draw.MolToImage(mol, subImgSize, legend=mol.GetProp("_Name")), (col * subImgSize[0], row * subImgSize[1]))
img.save(sys.argv[2]) # Path to a output PNG or JPG file.
