#!/bin/sh
if [ ! -d "bin" ]; then
  mkdir bin
fi

if [ ! -d "out" ]; then
  mkdir out
fi

echo "Building gaps!"
make -C gaps

echo "Building fetregister and fetbenchmark!"
make 
