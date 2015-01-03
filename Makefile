
all: aot

aot_run:
	g++ main_aot_run.cpp halide_haze_removal.o `libpng-config --cflags --ldflags` -lpthread -ldl -o bin/main_aot

aot:
	g++ main_aot_generate.cpp -g -I ../include -L ../bin -lHalide -lpthread -ldl -lz -o bin/main_aot_generate

aot_generate:
	LD_LIBRARY_PATH=../bin ./bin/main_aot_generate

jit:
	g++ main_jit.cpp -g -I ../include -L ../bin -lHalide `libpng-config --cflags --ldflags` -lpthread -ldl -o bin/main_jit

jit_test:
	LD_LIBRARY_PATH=../bin ./bin/main_jit images/cones.png
