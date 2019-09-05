CC=g++

bin/calcDescriptors: obj/calcDescriptors.o
	${CC} -o $@ $< -L${RDKIT_ROOT}/lib -lRDKitMolAlign -lRDKitForceFieldHelpers -lRDKitDescriptors -lRDKitSubstructMatch -lRDKitFileParsers -lRDKitSmilesParse -lRDKitDepictor -lRDKitGraphMol -lRDKitRDGeneral

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
