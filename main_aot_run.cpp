
#include <stdio.h>

#include "halide_haze_removal.h"
#include "static_image.h"
#include "image_io.h"

int main(int argc, char **argv){

    Image<uint8_t> input =  load<uint8_t>(argv[1]);
    Image<uint8_t> output(input.width(), input.height(), input.channels());


    halide_haze_removal(input, output);

    save(output, argv[2]);

    return 0;
}  // int main
