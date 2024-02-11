#include <iostream>
#include <dlfcn.h>

using namespace std;

extern "C" typedef double* (*run_a_step_t)(double* inputs, int num_of_inputs);

extern "C" typedef void* (*set_model_parameters_t)(double parameters, int num_of_parameters);

extern "C" typedef double* (*get_model_parameters_t)();

extern "C" typedef void* (*LoadModel_t)();
extern "C" typedef void* (*InitModel_t)();
extern "C" typedef void* (*TermModel_t)();

int main() {

	const int sim_len = 10000;
	double test_Output_log[sim_len];
	void* lib = dlopen("./libinterface_simulink.so", RTLD_LAZY | RTLD_GLOBAL);
	if (!lib)
	{
		fprintf(stderr, "dlopen failed: %s\n", dlerror());
	};

	// Load functions
	run_a_step_t run_a_step = (run_a_step_t)dlsym(lib, "run_a_step");
	set_model_parameters_t set_model_parameters=(set_model_parameters_t)dlsym(lib, "set_model_parameters");
	get_model_parameters_t get_model_parameters=(get_model_parameters_t)dlsym(lib, "get_model_parameters");
	LoadModel_t LoadModel = (LoadModel_t)dlsym(lib, "LoadModel");
	InitModel_t InitModel = (InitModel_t)dlsym(lib, "InitModel");
	TermModel_t TermModel = (TermModel_t)dlsym(lib, "TermModel");

	//Load and init model
	LoadModel();
	InitModel();

	double* a;
	double* res;
	double input_mock=1;
	a = &input_mock;
	double const_para;
	const_para = 0.5;
	set_model_parameters(const_para, 2);



	for (int i = 0; i < sim_len; i++) {
		res=run_a_step(a,1);
		test_Output_log[i] = *res;
	}

	//Conclude
	TermModel();

	string s;
	cout << "This is the caller lib." << endl;
	cin >> s;
	dlclose(lib);
	return 0;
}