rm -rf build/
mkdir build
cd build

ifort -c ../main.f90 -o f.o
icpc -c ../temp.cpp -o c.o
ifort f.o c.o -o out.o -lc -lstdc++

./out.o