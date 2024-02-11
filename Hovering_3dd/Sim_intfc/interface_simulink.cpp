// interface_simulink.cpp : Defines the entry point for the application.
//
/*==================================*
 * Global data local to this module *
 *==================================*/


#include "test.h"
#include "interface_simulink.h"

using namespace std;

/* Golbal Definitions*/

#ifndef MODEL
# error Must specify a model name.  Define MODEL=name.
#else
/* create generic macros that work with any model */
# define EXPAND_CONCAT(name1,name2) name1 ## name2
# define CONCAT(name1,name2) EXPAND_CONCAT(name1,name2)
# define MODEL_INITIALIZE CONCAT(MODEL,_initialize)
# define MODEL_STEP       CONCAT(MODEL,_step)
# define MODEL_TERMINATE  CONCAT(MODEL,_terminate)
# define RT_MDL_TYPE      CONCAT(MODEL,_M_TYPE)
#endif

#ifndef MULTITASKING
static boolean_T OverrunFlags[1];    /* ISR overrun flags */
static boolean_T eventFlags[1];      /* necessary for overlapping preemption */
#else
static boolean_T OverrunFlags[NUMST];
static boolean_T eventFlags[NUMST];
#endif

const char* RT_MEMORY_ALLOCATION_ERROR = "memory allocation error";


/* Function: rt_InitModel ====================================================
 *
 * Abstract:
 *   Initializes the model and prints the error status
 *
 */

static void rt_InitModel(RT_MDL_TYPE* S)
{
#if defined(MULTITASKING)
    int i;
    for (i = 0; i < NUMST; i++) {
        OverrunFlags[i] = 0;
        eventFlags[i] = 0;
    }
#else
    OverrunFlags[0] = 0;
    eventFlags[0] = 0;
#endif

    /************************
     * Initialize the model *
     ************************/
    MODEL_INITIALIZE(S);
}

/* Function: rt_TermModel ====================================================
 *
 * Abstract:
 *   Terminates the model and prints the error status
 *
 */
static int_T rt_TermModel(RT_MDL_TYPE* S)
{
    const char_T* errStatus = (const char_T*)(rtmGetErrorStatus(S));
    int_T i = 0;

    if (errStatus != NULL && strcmp(errStatus, "Simulation finished")) {
        (void)printf("%s\n", errStatus);
#if defined(MULTITASKING)
        for (i = 0; i < NUMST; i++) {
            if (OverrunFlags[i]) {
                (void)printf("ISR overrun - sampling rate too"
                    "fast for sample time index %d.\n", i);
            }
        }
#else
        if (OverrunFlags[i]) {
            (void)printf("ISR overrun - base sampling rate too fast.\n");
        }
#endif
        MODEL_TERMINATE(S);
        return(1);
    }

    MODEL_TERMINATE(S);
    return(0);
}

int main()
{
    const int sim_len = 10000;
    real_T test_Output_log[sim_len];


    RT_MDL_TYPE* S;
    const char_T* errmsg;
    int_T ret;

    /* External mode */
    //rtParseArgsForExtMode(argc, argv);

    /*******************************************
     * warn if the model will run indefinitely *
     *******************************************/
#if MAT_FILE==0 && EXT_MODE==0
    printf("warning: the simulation will run with no stop time; "
        "to change this behavior select the 'MAT-file logging' option\n");
    fflush(NULL);
#endif

    /************************
     * Initialize the model *
     ************************/

    S = MODEL();
    if (S == NULL) {
        (void)fprintf(stderr, "Memory allocation error during model "
            "registration");
        return(1);
    }
    errmsg = (const char_T*)(rtmGetErrorStatus(S));
    if (errmsg != NULL) {
        (void)fprintf(stderr, "Error during model registration: %s\n", errmsg);
        MODEL_TERMINATE(S);
        return(1);
    }

    if (errmsg != NULL) {
        (void)fprintf(stderr, "Error starting data logging: %s\n", errmsg);
        MODEL_TERMINATE(S);
        return(1);
    }


    (void)printf("\n** starting the model **\n");

    rt_InitModel(S);
    test_Input = 1.0;
    test_P.Constant_Value = 2.0;
    for (int i = 0; i < sim_len; i++) {
        MODEL_STEP(S);
        test_Output_log[i] = test_Output;
        if (i > 400) {
            test_P.Constant_Value = 10.0;
        }
    }

    ret = rt_TermModel(S);

	string s;
	cout << "Hello CMake." << endl;
	cin >>s;
	return 0;
}
