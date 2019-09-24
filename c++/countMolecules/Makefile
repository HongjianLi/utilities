CC=g++

all: bin/countMoleculesSdf bin/countMoleculesMol2

bin/countMoleculesSdf: obj/countMoleculesSdf.o
	${CC} -o $@ $<

bin/countMoleculesMol2: obj/countMoleculesMol2.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
