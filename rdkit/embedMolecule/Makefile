CC=g++

bin/embedMolecule: obj/embedMolecule.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitDistGeomHelpers -lRDKitFileParsers -lRDKitGraphMol

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
