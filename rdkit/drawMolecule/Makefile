CC=g++

bin/drawMolecule: obj/drawMolecule.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitMolDraw2D -lRDKitChemReactions -lRDKitFileParsers -lRDKitSmilesParse -lRDKitDepictor -lRDKitGraphMol -lRDKitRDGeneral

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
