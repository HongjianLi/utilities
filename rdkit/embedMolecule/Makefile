CC=g++

all: bin/embedMolecule bin/embedMultipleConfs bin/embedEncode

bin/embedMolecule: obj/embedMolecule.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitDistGeomHelpers -lRDKitFileParsers -lRDKitGraphMol -lRDKitRDGeneral

bin/embedMultipleConfs: obj/embedMultipleConfs.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitDistGeomHelpers -lRDKitFileParsers -lRDKitGraphMol -lRDKitRDGeneral

bin/embedEncode: obj/embedEncode.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitMolAlign -lRDKitDistGeomHelpers -lRDKitForceFieldHelpers -lRDKitSubstructMatch -lRDKitFileParsers -lRDKitSmilesParse -lRDKitDepictor -lRDKitGraphMol -lRDKitDistGeometry -lRDKitRDGeneral

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
