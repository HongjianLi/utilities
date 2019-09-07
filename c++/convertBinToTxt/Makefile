CC=g++

convertBinToTxt: bin/convertBinToTxt-f32 bin/convertBinToTxt-u32 bin/convertBinToTxt-u64 bin/convertBinToTxt-u16

bin/convertBinToTxt-f32: obj/convertBinToTxt-f32.o
	${CC} -o $@ $<

bin/convertBinToTxt-u32: obj/convertBinToTxt-u32.o
	${CC} -o $@ $<

bin/convertBinToTxt-u64: obj/convertBinToTxt-u64.o
	${CC} -o $@ $<

bin/convertBinToTxt-u16: obj/convertBinToTxt-u16.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
