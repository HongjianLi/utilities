CC=g++

all: bin/calcFPSimSmi bin/calcFPSimSdf

bin/calcFPSimSmi: obj/calcFPSimSmi.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFingerprints -lRDKitSmilesParse -lRDKitGraphMol -lRDKitDataStructs -lRDKitRDGeneral

bin/calcFPSimSdf: obj/calcFPSimSdf.o
	${CC} -o $@ $^ -L${RDKIT_ROOT}/lib -lRDKitFingerprints -lRDKitFileParsers -lRDKitGraphMol -lRDKitDataStructs -lRDKitRDGeneral

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall -I${RDKIT_ROOT}/include/rdkit -I${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
