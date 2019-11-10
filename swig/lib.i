%module lib

%{
#define SWIG_FILE_WITH_INIT
#include "lib.h"
%}

%include "std_vector.i"

namespace std {
    %template(DoubleVector)  vector<double>;
}

%callback("%s_cb");
void eulerIntegrator(vector<double>, ode_sys);
%nocallback;


%include "lib.h"
