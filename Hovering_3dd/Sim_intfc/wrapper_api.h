#include "uavNavigationRL2D_2022a_icra.h"
#include <iostream>
#include <string>


extern "C" double* run_a_step(double* inputs,int num_of_inputs);

extern "C" void* set_model_parameters(double parameters, int idx);

extern "C" double* get_model_parameters();

extern "C" void* LoadModel();
extern "C" void* InitModel();
extern "C" void* TermModel();
