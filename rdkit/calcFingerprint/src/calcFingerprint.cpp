#include <iostream>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
//#include <GraphMol/Fingerprints/MACCS.h>
//#include <GraphMol/Fingerprints/Fingerprints.h>
using namespace std;
using namespace RDKit;

int main(const int argc, const char* argv[])
{
	if (argc < 2) {
		cout << argv[0] << " input.sdf" << endl;
		return 1;
	}
	SDMolSupplier supplier(argv[1], true, false, true); // Default arguments: const std::string &fileName, bool sanitize=true, bool removeHs=true, strictParsing=true
	while (!supplier.atEnd()) {
		// Obtain a pointer to the current molecule with heavy atoms only.
		const unique_ptr<ROMol> mol_ptr(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		// Obtain a reference to the molecule to avoid writing *mol_ptr.
		auto& mol = *mol_ptr;
		// Obtain the molecule name from the smi file
		const string name = mol.getProp<string>("_Name"); // mol.getPropList() https://www.rdkit.org/docs/cppapi/classRDKit_1_1RDProps.html#ad63e121bf0725c67b9dd689b4c889bd5
		// Get fingerprints
		const unique_ptr<ExplicitBitVect> fp(MorganFingerprints::getFingerprintAsBitVect(mol, 2, 2048)); // mol, radius, nBits, ... https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MorganFingerprints.html#a1bf757d66784abf5a4ebf2c869be8261;
//		const unique_ptr<SparseIntVect<uint32_t>> fpInt(MorganFingerprints::getFingerprint(mol, 2)); // mol, radius, ... This fingerprint in SparseIntVect type will serialize to a much larger string than using ExplicitBitVect.
//		MACCSFingerprints::getFingerprintAsBitVect(mol); // https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MACCSFingerprints.html#a763a8cf42c312fead7cdb351e20771ef
//		RDKFingerprintMol(mol);     // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
//		LayeredFingerprintMol(mol); // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
//		PatternFingerprintMol(mol); // https://www.rdkit.org/docs/cppapi/Fingerprints_8h.html
		// Write the fingerprint
		const auto fpString = fp->toString();
		cout << name << '\t' << fpString.size() << '\t' << fpString << endl;
	}
}
