CC=g++

all: bin/countMoleculesSdf

bin/countMoleculesSdf: obj/countMoleculesSdf.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
