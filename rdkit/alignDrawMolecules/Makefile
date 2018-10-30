CC=g++

bin/alignDrawMolecules: obj/alignDrawMolecules.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitMolDraw2D -lRDKitGraphMol -lRDKitFileParsers -lRDKitDepictor -lRDKitRDGeneral -lRDKitFMCS -lRDKitSmilesParse -lRDKitSubstructMatch

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit -I${BOOST_ROOT}/include
