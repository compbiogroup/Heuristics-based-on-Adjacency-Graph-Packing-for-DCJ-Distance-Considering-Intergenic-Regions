.PHONY: clean test debug all

all:
	cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
	cmake --build build
	mv build/dec .

clean:
	cd build && make clean
