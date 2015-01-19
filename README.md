NOTE
====
*Makefile* is very, very naive style. But the most importance to note is path to halide library.

Binaries should be easily built if we have correct path to halide libary(`-I` and `-L` value of g++).

Build and Usage
===============
* `make aot`

   compile *bin/main_aot_generate* for `make aot_generate`
* `make aot_generate`

   generate *halide_haze_removal.o* and *halide_haze_removal.h* for `make aot_run`
* `make aot_run`

   compile AOT binary out. We will have *bin/main_aot*. Its usage: `main_aot <input_png> <output_png>`
* `make jit`

   compile JIT binary out. We will have *bin/main_jit*. It's usage is as `make jit_test`
