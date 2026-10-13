#!/bin/sh
# Compile the state exchange library inside the GLM runtime image.
# Arguments: path to the folder with glm.h and glm_types.h (GLM3 commit 8453b1d), output path.
set -eu
gcc -O2 -Wall -shared -fPIC -I"$1" native/native_exchange.c -o "$2" -ldl -lnetcdf -lm
