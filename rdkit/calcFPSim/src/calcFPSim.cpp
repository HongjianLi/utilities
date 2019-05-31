#include <iostream>
#include <fstream>
#include <GraphMol/FileParsers/MolSupplier.h>
//#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
//#include <GraphMol/Fingerprints/MACCS.h>
//#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/BitOps.h>
using namespace std;
using namespace RDKit;

int main(const int argc, const char* argv[])
{
	if (argc < 3) {
//		cout << argv[0] << " Cc1ccccc1 Cc1ncccc1" << endl;
		cout << argv[0] << " mol0.sdf mol1.sdf" << endl;
		return 1;
	}
	const array<unique_ptr<ROMol>, 2> mols {
		unique_ptr<ROMol>(SDMolSupplier(argv[1]).next()),
		unique_ptr<ROMol>(SDMolSupplier(argv[2]).next()),
//		unique_ptr<ROMol>(static_cast<ROMol*>(SmilesToMol(argv[1]))),
//		unique_ptr<ROMol>(static_cast<ROMol*>(SmilesToMol(argv[2]))),
	};
	// Get fingerprints. Use Morgan fingerprints by default. Other fingerprints include MACCS, RDK, Layered, Pattern.
	const array<unique_ptr<ExplicitBitVect>, 2> fps {
		unique_ptr<ExplicitBitVect>(MorganFingerprints::getFingerprintAsBitVect(*mols[0], 2, 2048)), // mol, radius, nBits, ... https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MorganFingerprints.html
		unique_ptr<ExplicitBitVect>(MorganFingerprints::getFingerprintAsBitVect(*mols[1], 2, 2048)),
//		unique_ptr<ExplicitBitVect>(MACCSFingerprints::getFingerprintAsBitVect(*mols[0])),
//		unique_ptr<ExplicitBitVect>(MACCSFingerprints::getFingerprintAsBitVect(*mols[1])),
	};
//		MACCSFingerprints::getFingerprintAsBitVect(mol); // https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MACCSFingerprints.html
//		RDKFingerprintMol(mol);     // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
//		LayeredFingerprintMol(mol); // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
//		PatternFingerprintMol(mol); // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
	// Calculate the similarity of two fingerprints. Use Tanimoto similarity by default. Other similarities include Dice, Cosine, Tversky, etc. https://www.rdkit.org/docs/cppapi/BitOps_8h.html
	cout << TanimotoSimilarity(*fps[0], *fps[1]) << endl;
}
