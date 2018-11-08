#include <iostream>
#include <boost/shared_ptr.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
using namespace std;
using namespace boost;
using namespace RDKit;
using namespace RDDepict;

int main(const int argc, const char* argv[]) {
	if (argc < 2) {
		cout << argv[0] << " input.smi" << endl;
		return 1;
	}
	vector<boost::shared_ptr<ROMol>> mols;
	SmilesMolSupplier supplier(argv[1], "\t", 1, 0, false); // Default parameters are: const std::string &fileName, const std::string &delimiter=" \t", int smilesColumn=0, int nameColumn=1, bool titleLine=true, bool sanitize=true)
	while (!supplier.atEnd()) {
		mols.emplace_back(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
	}
	const auto mcs = findMCS(mols); // Find the maximum common substructure
	unique_ptr<ROMol> ref(SmartsToMol(mcs.SmartsString));
	compute2DCoords(*ref);
	const auto& refCnf = ref->getConformer();
	const size_t panelWidth = 300;
	const size_t panelHeight = 300;
	MolDraw2DSVG drawer(panelWidth * mols.size(), panelHeight, panelWidth, panelHeight); // Horizontally place the molecules. https://www.rdkit.org/docs/cppapi/classRDKit_1_1MolDraw2DSVG.html
	vector<ROMol*> vmols;
	vector<string> legends;
	vector<vector<int>> vhighlight_atoms;
	for (const auto& mol : mols) {
		vector<pair<int, int>> matchVect;
		SubstructMatch(*mol, *ref, matchVect);
		map<int, Point2D> coordMap;
		vector<int> highlight_atoms;
		for (const auto& match : matchVect)
		{
			const auto& pt3 = refCnf.getAtomPos(match.first);
			coordMap[match.second] = Point2D(pt3.x, pt3.y); // pt3.z == 0.0
			highlight_atoms.push_back(match.second);
		}
		compute2DCoords(*mol, &coordMap);
		vmols.push_back(mol.get());
		legends.emplace_back(mol->getProp<string>("_Name"));
		vhighlight_atoms.emplace_back(std::move(highlight_atoms));
	}
	drawer.drawMolecules(vmols, &legends, &vhighlight_atoms); // https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/MolDraw2D/MolDraw2D.h
	drawer.finishDrawing();
	cout << drawer.getDrawingText();
}
