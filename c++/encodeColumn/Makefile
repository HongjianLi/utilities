CC=g++

encodeColumn: bin/encodeColumn-f32 bin/encodeColumn-u32 bin/encodeColumn-u64 bin/encodeColumn-u8

bin/encodeColumn-f32: obj/encodeColumn-f32.o
	${CC} -o $@ $<

bin/encodeColumn-u32: obj/encodeColumn-u32.o
	${CC} -o $@ $<

bin/encodeColumn-u64: obj/encodeColumn-u64.o
	${CC} -o $@ $<

bin/encodeColumn-u8: obj/encodeColumn-u8.o
	${CC} -o $@ $<

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++2a -O2 -Wall

clean:
	rm -f bin/* obj/*
