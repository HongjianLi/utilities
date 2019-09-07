#include <iostream>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
//#include <GraphMol/FileParsers/MolWriters.h>
using namespace std;
using namespace RDKit;

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << argv[0] << " input.sdf" << endl;
		return 1;
	}
//	SmilesWriter writer(argv[2], "\t", "Name", false); // This writer class does not allow to customize column ordering. It always outputs SMILES in the 1st column and name in the 2nd column. Default parameters are: const std::string &fileName, const std::string &delimiter=" ", const std::string &nameHeader="Name", bool includeHeader=true, bool isomericSmiles=true, bool kekuleSmiles=false)
	SDMolSupplier supplier(argv[1]); // Default arguments: const std::string &fileName, bool sanitize=true, bool removeHs=true, strictParsing=true
	while (!supplier.atEnd()) {
		// Obtain a pointer to the current molecule with heavy atoms only.
		const unique_ptr<ROMol> mol_ptr(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		// Obtain a reference to the molecule to avoid writing *mol_ptr.
		auto& mol = *mol_ptr;
		// Obtain the molecule name from the smi file
		const string name = mol.getProp<string>("_Name"); // mol.getPropList() https://www.rdkit.org/docs/cppapi/classRDKit_1_1RDProps.html#ad63e121bf0725c67b9dd689b4c889bd5
		cout
			<< name << '\t'
			<< MolToSmiles(mol) << '\n' // Default parameters are: const ROMol& mol, bool doIsomericSmiles = true, bool doKekule = false, int rootedAtAtom = -1, bool canonical = true, bool allBondsExplicit = false, bool allHsExplicit = false, bool doRandom = false. https://www.rdkit.org/docs/cppapi/namespaceRDKit.html#a3636828cca83a233d7816f3652a9eb6b
		;
//		writer.write(mol);
	}
}
