swig -c++ -python -shadow lib.i
g++  -std=c++11 -O2 -fPIC -fopenmp -c lib.cpp 
g++  -std=c++11 -O2 -fPIC -c lib_wrap.cxx -I /usr/include/python2.7 
g++  -shared lib_wrap.o lib.o -o _lib.so 
