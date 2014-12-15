

all:
	g++ main_jit.cpp -g -I ../include -L ../bin -lHalide `libpng-config --cflags --ldflags` -lpthread -ldl -o bin/main_jit

jit_test:
	LD_LIBRARY_PATH=../bin ./bin/main_jit images/cones.png
