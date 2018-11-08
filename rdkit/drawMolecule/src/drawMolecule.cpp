#include <fstream>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
using namespace std;
using namespace RDKit;
using namespace RDDepict;

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << argv[0] << " input.smi" << endl;
		return 1;
	}
	SmilesMolSupplier supplier(argv[1], "\t", 1, 0, false); // Default parameters are: const std::string &fileName, const std::string &delimiter=" \t", int smilesColumn=0, int nameColumn=1, bool titleLine=true, bool sanitize=true)
	while (!supplier.atEnd()) {
		// Obtain a pointer to the current molecule with heavy atoms only.
		const unique_ptr<ROMol> mol_ptr(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		// Obtain a reference to the molecule to avoid writing *mol_ptr.
		auto& mol = *mol_ptr;
		// Obtain the molecule name from the smi file
		const string name = mol.getProp<string>("_Name");
		// Compute a 2D conformer
		compute2DCoords(mol);
		// Draw to SVG
		ofstream ofs(name + ".svg");
		MolDraw2DSVG drawer(600, 600, ofs); // width, height, output
		drawer.drawMolecule(mol);
		drawer.finishDrawing();
	}
}
