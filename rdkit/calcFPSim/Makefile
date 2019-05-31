CC=g++

bin/calcFPSim: obj/calcFPSim.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFingerprints -lRDKitFileParsers -lRDKitSmilesParse -lRDKitGraphMol -lRDKitDataStructs -lRDKitRDGeneral

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
