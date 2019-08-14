CC=g++

calcFooter: bin/calcFooterTxt bin/calcFooterSdf bin/calcFooterPdbqt

bin/calcFooterTxt: obj/calcFooterTxt.o
	${CC} -o $@ $<

bin/calcFooterSdf: obj/calcFooterSdf.o
	${CC} -o $@ $<

bin/calcFooterPdbqt: obj/calcFooterPdbqt.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
