CC=g++

bin/calcRmsdPdbqt: obj/calcRmsdPdbqt.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall

clean:
	rm -f bin/* obj/*
