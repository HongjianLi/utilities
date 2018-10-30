#include <iostream>
#include <fstream>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
//#include <GraphMol/Fingerprints/MACCS.h>
//#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <DataStructs/BitOps.h>
using namespace std;
using namespace RDKit;

int main(const int argc, const char* argv[])
{
	if (argc < 2) {
		cout << argv[0] << " input.smi" << endl;
		return 1;
	}
	SmilesMolSupplier supplier(argv[1], "\t", 0, 1, false); // Default arguments: const std::string &fileName, const std::string &delimiter=" \t", int smilesColumn=0, int nameColumn=1, bool titleLine=true, bool sanitize=true)
	while (!supplier.atEnd()) {
		// Obtain a pointer to the current molecule with heavy atoms only.
		const unique_ptr<ROMol> mol_ptr(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		// Obtain a reference to the molecule to avoid writing *mol_ptr.
		auto& mol = *mol_ptr;
		// Obtain the molecule name from the smi file
		string name;
		mol.getProp("_Name", name);
		// Get fingerprints
		const unique_ptr<ExplicitBitVect> fp(MorganFingerprints::getFingerprintAsBitVect(mol, 2, 2048)); // mol, radius, nBits, ... https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MorganFingerprints.html#a1bf757d66784abf5a4ebf2c869be8261;
//		MACCSFingerprints::getFingerprintAsBitVect(mol); // https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MACCSFingerprints.html#a763a8cf42c312fead7cdb351e20771ef
//		RDKFingerprintMol(mol);     // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
//		LayeredFingerprintMol(mol); // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
//		PatternFingerprintMol(mol); // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
		// Write the fingerprint to *.fp files
		ofstream ofs(name + ".fp");
		ofs << fp->toString();
	}
//	cout << TanimotoSimilarity(*fp0, *fp1) << endl;
//	cout << DiceSimilarity(*fp0, *fp1) << endl;
}
