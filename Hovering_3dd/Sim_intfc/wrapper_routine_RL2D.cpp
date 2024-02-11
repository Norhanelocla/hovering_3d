#include "wrapper_api.h"

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

/* Local Definitions*/

RT_MDL_TYPE* S;
const char_T* errmsg;
int_T ret;
double* out_res;



double* run_a_step(double* inputs, int num_of_inputs) {
    //num_of_steps++;
    //num_of_steps_el = &num_of_steps;
    for (int i = 0; i < num_of_inputs; i++) {
        uavNavigationRL2_velocity_cmd_x = inputs[0];
        uavNavigationRL2_velocity_cmd_y = inputs[1];
        uavNavigationRL2_velocity_cmd_z= inputs[2];
        uavNavigationRL2D_2022a_u_cmd_z=inputs[3];
    }
    MODEL_STEP(S);
    out_res = &uavNavigationRL2D_states_output[0];
    return out_res;
}
/*
void* set_model_parameters(double parameters, int idx) {
    if (idx == 0) {
        uavNavigationRL2D_2021b_P.upper_limit_y_Value = parameters;
    }
    else if (idx == 1) {
        uavNavigationRL2D_2021b_P.upper_limit_x_Value = parameters;
    }
    else if (idx == 2) {
        uavNavigationRL2D_2021b_P.lower_limit_y_Value = parameters;
    }
    else if (idx == 3) {
        uavNavigationRL2D_2021b_P.lower_limit_x_Value = parameters;
    }
    else if (idx == 4) {
        uavNavigationRL2D_2021b_P.uav_initial_x_position = parameters;
    }
    else if (idx == 5) {
        uavNavigationRL2D_2021b_P.uav_initial_y_position = parameters;
    }
    else if (idx == 6) {
        uavNavigationRL2D_2021b_P.uav_initial_z_position = parameters;
    }
    else if (idx == 7) {
        uavNavigationRL2D_2021b_P.T_motor = parameters;
    }
    else if (idx == 8) {
        uavNavigationRL2D_2021b_P.K_roll = parameters;
    }
    else if (idx == 9) {
        uavNavigationRL2D_2021b_P.T_roll = parameters;
    }
    else if (idx == 10) {
        uavNavigationRL2D_2021b_P.tau_roll = parameters;
    }
    else if (idx == 11) {
        uavNavigationRL2D_2021b_P.K_pitch = parameters;
    }
    else if (idx == 12) {
        uavNavigationRL2D_2021b_P.T_pitch = parameters;
    }
    else if (idx == 13) {

        uavNavigationRL2D_2021b_P.tau_pitch = parameters;
    }
    else if (idx == 14) {
        uavNavigationRL2D_2021b_P.K_z = parameters;
    }
    else if (idx == 15) {
        uavNavigationRL2D_2021b_P.T_aero = parameters;
    }
    else if (idx == 16) {
        uavNavigationRL2D_2021b_P.tau_z = parameters;
    }
    else if (idx == 17) {
        uavNavigationRL2D_2021b_P.tau_x = parameters;

    }
    else if (idx == 18) {
        uavNavigationRL2D_2021b_P.tau_y = parameters;
    }
    else if (idx == 19) {
        uavNavigationRL2D_2021b_P.kp_roll = parameters;
    }
    else if (idx == 20) {
        uavNavigationRL2D_2021b_P.kd_roll = parameters;
    }
    else if (idx == 21) {
        uavNavigationRL2D_2021b_P.kp_pitch = parameters;
    }
    else if (idx == 22) {
        uavNavigationRL2D_2021b_P.kd_pitch = parameters;
    }
    else if (idx == 23) {
        uavNavigationRL2D_2021b_P.kp_yaw = parameters;
    }
    else if (idx == 24) {
        uavNavigationRL2D_2021b_P.kd_yaw = parameters;
    }
    else if (idx == 25) {
        uavNavigationRL2D_2021b_P.kp_x = parameters;
    }
    else if (idx == 26) {
        uavNavigationRL2D_2021b_P.kd_x = parameters;
    }
    else if (idx == 27) {
        uavNavigationRL2D_2021b_P.kp_y = parameters;
    }
    else if (idx == 28) {
        uavNavigationRL2D_2021b_P.kd_y = parameters;
    }
    else if (idx == 29) {
        uavNavigationRL2D_2021b_P.kp_z = parameters;
    }
    else if (idx == 30) {
        uavNavigationRL2D_2021b_P.kd_z = parameters;
    }

    else {
        printf("warning: out of bounds index \n");
    } 
}
*/

double* get_model_parameters() {

}

/* Function: rt_InitModel ====================================================
 *
 * Abstract:
 *   Initializes the model and prints the error status
 *
 */

void rt_InitModel(RT_MDL_TYPE* S){
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
int_T rt_TermModel(RT_MDL_TYPE* S)
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

void* InitModel() {
    //num_of_steps = 0;
    rt_InitModel(S);
}

void* TermModel() {
    //num_of_steps = 0;
    int_T res = rt_TermModel(S);
}

void* LoadModel() {

#if MAT_FILE==0 && EXT_MODE==0
    printf("warning: the simulation will run with no stop time; "
        "to change this behavior select the 'MAT-file logging' option\n");
    fflush(NULL);
#endif
    S = MODEL();
    if (S == NULL) {
        (void)fprintf(stderr, "Memory allocation error during model "
            "registration");
    }
    errmsg = (const char_T*)(rtmGetErrorStatus(S));
    if (errmsg != NULL) {
        (void)fprintf(stderr, "Error during model registration: %s\n", errmsg);
        MODEL_TERMINATE(S);
    }

    if (errmsg != NULL) {
        (void)fprintf(stderr, "Error starting data logging: %s\n", errmsg);
        MODEL_TERMINATE(S);
    }
    (void)printf("\n** starting the model **\n");
}