#include <iostream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <DataStructs/BitOps.h>
using namespace std;
using namespace RDKit;

int main(const int argc, const char* argv[])
{
	if (argc < 3) {
		cout << argv[0] << " Cc1ccccc1 Cc1ncccc1" << endl;
		return 1;
	}
	const array<unique_ptr<ROMol>, 2> mols {
		unique_ptr<ROMol>(static_cast<ROMol*>(SmilesToMol(argv[1]))),
		unique_ptr<ROMol>(static_cast<ROMol*>(SmilesToMol(argv[2]))),
	};
	// Get fingerprints. Use Morgan fingerprints by default. Other fingerprints include MACCS, RDK, Layered, Pattern.
	const array<unique_ptr<ExplicitBitVect>, 2> fps {
		unique_ptr<ExplicitBitVect>(MorganFingerprints::getFingerprintAsBitVect(*mols[0], 2, 2048)), // mol, radius, nBits, ... https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MorganFingerprints.html
		unique_ptr<ExplicitBitVect>(MorganFingerprints::getFingerprintAsBitVect(*mols[1], 2, 2048)),
	};
	// Calculate the similarity of two fingerprints. Use Tanimoto similarity by default. Other similarities include Dice, Cosine, Tversky, etc. https://www.rdkit.org/docs/cppapi/BitOps_8h.html
	cout << TanimotoSimilarity(*fps[0], *fps[1]) << endl;
}
