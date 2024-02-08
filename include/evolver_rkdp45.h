#ifndef __RK__
#define __RK__
#include "common.h"
/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int evolver_rkdp45(
       int (*derivs)(double x,double * y,double * dy,
         void * parameters_and_workspace, ErrorMsg error_message),
       double t_ini,
       double t_final,
       double * y_inout,
       int* used_in_output,
       int neq,
       void * parameters_and_workspace_for_derivs,
       double rtol,
       double minimum_variation,
       int (*timescale_and_approximation)(double x,
                  void * parameters_and_workspace,
                  double * timescales,
                  ErrorMsg error_message),
       double timestep_over_timescale,
       double* t_vec,
       int tres,
       int (*output)(double x,double y[],double dy[],int index_x,void * parameters_and_workspace,
         ErrorMsg error_message),
       int (*print_variables)(double x, double y[], double dy[], void *parameters_and_workspace,
            ErrorMsg error_message),
       ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
