CC=g++

decodeColumn: bin/decodeColumn-f32 bin/decodeColumn-u32 bin/decodeColumn-u64 bin/decodeColumn-u8

bin/decodeColumn-f32: obj/decodeColumn-f32.o
	${CC} -o $@ $<

bin/decodeColumn-u32: obj/decodeColumn-u32.o
	${CC} -o $@ $<

bin/decodeColumn-u64: obj/decodeColumn-u64.o
	${CC} -o $@ $<

bin/decodeColumn-u8: obj/decodeColumn-u8.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
