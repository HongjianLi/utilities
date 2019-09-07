CC=g++

convertTxtToBin: bin/convertTxtToBin-f32 bin/convertTxtToBin-u32 bin/convertTxtToBin-u64 bin/convertTxtToBin-u16

bin/convertTxtToBin-f32: obj/convertTxtToBin-f32.o
	${CC} -o $@ $<

bin/convertTxtToBin-u32: obj/convertTxtToBin-u32.o
	${CC} -o $@ $<

bin/convertTxtToBin-u64: obj/convertTxtToBin-u64.o
	${CC} -o $@ $<

bin/convertTxtToBin-u16: obj/convertTxtToBin-u16.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
