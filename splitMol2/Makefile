CC=g++

bin/splitMol2: obj/splitMol2.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -O2 -Wall

clean:
	rm -f bin/* obj/*
