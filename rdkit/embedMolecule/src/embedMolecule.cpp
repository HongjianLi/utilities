#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/FileParsers/MolWriters.h>
using namespace std;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << argv[0] << " input.smi" << endl;
		return 1;
	}
	SmilesMolSupplier supplier(argv[1], "\t", 1, 0, false); // Default arguments are: const std::string &fileName, const std::string &delimiter=" \t", int smilesColumn=0, int nameColumn=1, bool titleLine=true, bool sanitize=true)
	while (!supplier.atEnd()) {
		// Obtain a pointer to the current molecule with heavy atoms only.
		const unique_ptr<ROMol> smi_ptr(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		// Add hydrogens to the molecule. To get good conformations, it's almost always a good idea to add hydrogens to the molecule first. http://rdkit.org/docs/GettingStartedInPython.html#writing-molecules
		const unique_ptr<ROMol> mol_ptr(addHs(*smi_ptr));
		// Obtain a reference to the molecule to avoid writing *mol_ptr.
		auto& mol = *mol_ptr;
		// Obtain the molecule name from the smi file
		const string name = mol.getProp<string>("_Name");
		// Generate conformers with knowledge.
		EmbedParameters params(srETKDGv3);
		params.randomSeed = 209;
		const auto confId = EmbedMolecule(mol, params);
		// Check if conformers are generated.
		cout << name << '\t' << confId + 1 << endl;
		if (confId == -1) {
			continue;
		}
		// Create output streams.
		SDWriter writer(name + ".sdf");
		// Write the generated conformer.
		writer.write(mol);
	}
}
