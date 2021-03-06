#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <limits>
#include <filesystem>
#include <fstream>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
using namespace std;
using namespace std::filesystem;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace RDGeom;

template<typename T>
auto dist2(const T& p0, const T& p1)
{
	const auto d0 = p0[0] - p1[0];
	const auto d1 = p0[1] - p1[1];
	const auto d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << argv[0] << " input.smi" << endl;
		return 1;
	}

	// Initialize constants.
	const size_t num_references = 4;
	const size_t num_subsets = 5;
	const array<string, num_subsets> SubsetSMARTS
	{{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	}};

	// Wrap SMARTS strings to ROMol objects.
	array<unique_ptr<ROMol>, num_subsets> SubsetMols;
	for (size_t k = 0; k < num_subsets; ++k)
	{
		SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
	}

	// Initialize variables.
	array<vector<int>, num_subsets> subsets;
	array<Point3D, num_references> references;
	array<vector<double>, num_references> dista;

	// Create root directories.
	const path sdfRootPath = "sdf";
	const path usrRootPath = "usrcat";
	create_directory(sdfRootPath);
	create_directory(usrRootPath);

	// Loop over the input SMILES.
	SmilesMolSupplier supplier(argv[1], "\t", 1, 0, false); // Default arguments are: const std::string &fileName, const std::string &delimiter=" \t", int smilesColumn=0, int nameColumn=1, bool titleLine=true, bool sanitize=true)
	while (!supplier.atEnd()) {
		// Obtain a pointer to the current molecule with heavy atoms only.
		const unique_ptr<ROMol> smi_ptr(supplier.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		// Add hydrogens to the molecule. To get good conformations, it's almost always a good idea to add hydrogens to the molecule first. http://rdkit.org/docs/GettingStartedInPython.html#writing-molecules
		const unique_ptr<ROMol> mol_ptr(addHs(*smi_ptr));
		// Obtain a reference to the molecule to avoid writing *mol_ptr.
		auto& mol = *mol_ptr;
		// Obtain the molecule name from the smi file
		const auto name = mol.getProp<string>("_Name");

		// Generate conformers with knowledge.
		EmbedParameters params(ETKDGv2); // https://github.com/rdkit/rdkit/pull/1597
		params.randomSeed = 209;
		const auto confIds = EmbedMultipleConfs(mol, 4, params);

		// Check if conformers are generated.
		cout << name << '\t' << confIds.size() << endl;
		if (confIds.empty()) {
			continue;
		}

		// Get the number of points, excluding hydrogens.
		const auto num_points = mol.getNumHeavyAtoms();
		// Get the number of atoms, including hydrogens.
		const auto num_atoms = mol.getNumAtoms();

		// Categorize atoms into pharmacophoric subsets.
		for (size_t k = 0; k < num_subsets; ++k)
		{
			vector<vector<pair<int, int>>> matchVect;
			SubstructMatch(mol, *SubsetMols[k], matchVect);
			const auto num_matches = matchVect.size();
			auto& subset = subsets[k];
			subset.resize(num_matches);
			for (size_t i = 0; i < num_matches; ++i)
			{
				subset[i] = matchVect[i].front().second;
			}
		}
		const auto& subset0 = subsets.front();
		assert(subset0.size() == num_points);

		// Create molecule directories
		const auto sdfMolPath = sdfRootPath / name;
		const auto usrMolPath = usrRootPath / name;
		create_directory(sdfMolPath);
		create_directory(usrMolPath);

		// Loop over the generated conformers.
		for (const auto confId : confIds) {
			// Create output streams and write as sdf.
			const string filestem = to_string(confId + 1); // confId is zero-based. Adding by one to make [0, 1, 2, 3] to become [1, 2, 3, 4].
			SDWriter writer(sdfMolPath / (filestem + ".sdf"));
			writer.write(mol, confId);

			// Determine the reference points.
			const auto& conf = mol.getConformer(confId);
			for (auto& ref : references)
			{
				ref.x = ref.y = ref.z = 0;
			}
			auto& ctd = references[0];
			auto& cst = references[1];
			auto& fct = references[2];
			auto& ftf = references[3];
			for (const auto i : subset0)
			{
				const auto& a = conf.getAtomPos(i);
				ctd += a;
			}
			ctd /= num_points;
			auto cst_dist = numeric_limits<double>::max();
			auto fct_dist = numeric_limits<double>::lowest();
			auto ftf_dist = numeric_limits<double>::lowest();
			for (const auto i : subset0)
			{
				const auto& a = conf.getAtomPos(i);
				const auto this_dist = dist2(a, ctd);
				if (this_dist < cst_dist)
				{
					cst = a;
					cst_dist = this_dist;
				}
				if (this_dist > fct_dist)
				{
					fct = a;
					fct_dist = this_dist;
				}
			}
			for (const auto i : subset0)
			{
				const auto& a = conf.getAtomPos(i);
				const auto this_dist = dist2(a, fct);
				if (this_dist > ftf_dist)
				{
					ftf = a;
					ftf_dist = this_dist;
				}
			}

			// Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
			for (size_t k = 0; k < num_references; ++k)
			{
				const auto& reference = references[k];
				auto& distp = dista[k];
				distp.reserve(num_atoms); // Here use the number of all atoms instead of just heavy atoms because subset0[0] does not necessarily start from 0. In other words, hydrogens could appear before heavy atoms in the sdf file.
				for (size_t i = 0; i < num_points; ++i)
				{
					distp[subset0[i]] = sqrt(dist2(conf.getAtomPos(subset0[i]), reference));
				}
			}

			// Loop over pharmacophoric subsets and reference points.
			ofstream ofs(usrMolPath / (filestem + ".txt"));
			ofs << setprecision(6);
			for (const auto& subset : subsets)
			{
				const auto n = subset.size();
				for (size_t k = 0; k < num_references; ++k)
				{
					// Load distances from precalculated ones.
					const auto& distp = dista[k];
					vector<double> dists(n);
					for (size_t i = 0; i < n; ++i)
					{
						dists[i] = distp[subset[i]];
					}

					// Compute moments.
					array<double, 3> m{};
					if (n > 2)
					{
						const auto v = 1.0 / n;
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i];
							m[0] += d;
						}
						m[0] *= v;
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i] - m[0];
							m[1] += d * d;
						}
						m[1] = sqrt(m[1] * v);
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i] - m[0];
							m[2] += d * d * d;
						}
						m[2] = cbrt(m[2] * v);
					}
					else if (n == 2)
					{
						m[0] = 0.5 *     (dists[0] + dists[1]);
						m[1] = 0.5 * fabs(dists[0] - dists[1]);
					}
					else if (n == 1)
					{
						m[0] = dists[0];
					}

					// Write moments to file.
					ofs << m[0] << endl << m[1] << endl << m[2] << endl;
				}
			}
		}
	}
}
