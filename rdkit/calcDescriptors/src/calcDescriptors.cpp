#include <iostream>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
using namespace std;
using namespace RDKit;
using namespace RDKit::Descriptors;

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << argv[0] << " input.smi" << endl;
		return 1;
	}
	SmilesMolSupplier supplier(argv[1], "\t", 1, 0, false); // Default parameters are: const std::string &fileName, const std::string &delimiter=" \t", int smilesColumn=0, int nameColumn=1, bool titleLine=true, bool sanitize=true)
	cout << "ID	canonicalSMILES	molFormula	numAtoms	numHBD	numHBA	numRotatableBonds	numRings	exactMW	tPSA	clogP" << endl;
	while (!supplier.atEnd()) {
		// Obtain a pointer to the current molecule with heavy atoms only.
		const unique_ptr<ROMol> mol_ptr(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		// Obtain a reference to the molecule to avoid writing *mol_ptr.
		const auto& mol = *mol_ptr;
		// Obtain the molecule name from the smi file
		const string name = mol.getProp<string>("_Name"); // mol.getPropList() https://www.rdkit.org/docs/cppapi/classRDKit_1_1RDProps.html#ad63e121bf0725c67b9dd689b4c889bd5
		cout
			<< name << '\t'
			<< MolToSmiles(mol) << '\t' // Default parameters are: const ROMol& mol, bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1, bool canonical = true, bool allBondsExplicit = false, bool allHsExplicit = false, bool doRandom = false. https://www.rdkit.org/docs/cppapi/namespaceRDKit.html#a3636828cca83a233d7816f3652a9eb6b
			<< calcMolFormula(mol) << '\t'
			<< mol.getNumAtoms() << '\t'
//			<< mol.getNumHeavyAtoms() << '\t'
			<< calcNumHBD(mol) << '\t'
			<< calcNumHBA(mol) << '\t'
			<< calcNumRotatableBonds(mol) << '\t'
			<< calcNumRings(mol) << '\t'
			<< calcExactMW(mol) << '\t'
			<< calcTPSA(mol) << '\t'
			<< calcClogP(mol) << endl;
	}
}
