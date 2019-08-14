CC=g++

splitMolecule: bin/splitSdf bin/splitMol2 bin/splitPdbqt

bin/splitSdf: obj/splitSdf.o
	${CC} -o $@ $<

bin/splitMol2: obj/splitMol2.o
	${CC} -o $@ $<

bin/splitPdbqt: obj/splitPdbqt.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
