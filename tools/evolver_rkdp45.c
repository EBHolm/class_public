#include "evolver_rkdp45.h"

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
       ErrorMsg error_message){
  
  /** Handle options: */
  
  double abstol = 1e-15; // Same as in ndf15
  int stats[3] = {0, 0, 0};
  int verbose = 0;
  int output_return;

  double *dy,*err,*ynew,*ytemp, *ki;
  double h,absh,hmax,errmax,errtemp,hmin,hnew;
  double t,tnew;
  int tdir, i, j, k, idx=0;
  int s=7;
  double ci[s];
  double bi[s];
  double bi_diff[s];
  double ai[s][s];
  double bi_vec_y[s],bi_vec_dy[s];
  double threshold = abstol/rtol;
  int nofailed;
  int done = _FALSE_;
  double pow_grow=0.2;
  double rh,maxtmp;
  //Interpolation variables:
  double i01,i02,i03;
  double ixx[5][3];
  double ti,ss1,ss2,ss3,ss4,*yinterp,*dyinterp;
  dy = (double *) malloc(sizeof(double)*neq);
  ytemp = (double *) malloc(sizeof(double)*neq);
  yinterp = (double *) malloc(sizeof(double)*neq);
  dyinterp = (double *) malloc(sizeof(double)*neq);
  ynew = (double *) malloc(sizeof(double)*neq);
  err = (double *) malloc(sizeof(double)*neq);
  pow_grow = 0.2;
  ki = (double *) malloc(sizeof(double)*s*neq);
  /** Set method parameters for Runge-Kutta Dormand-Prince method:
      ------------------------------------------------------------------
  */
  ci[0] = 0.0; ci[1] = 0.2; ci[2] = 0.3; ci[3] = 0.8;
  ci[4] = 8.0/9.0; ci[5] = 1.0; ci[6] = 1.0;
  bi[0] = 35.0/384.0; bi[1] = 0.0; bi[2] = 500.0/1113.0;
  bi[3] = 125.0/192.0; bi[4] = -2187.0/6784.0; bi[5] = 11.0/84.0; bi[6] = 0.0;
  bi_diff[0] = 71.0/57600.0; bi_diff[1] = 0.0; bi_diff[2] = -71.0/16695.0;
  bi_diff[3] = 71.0/1920.0; bi_diff[4] = -17253.0/339200.0;
  bi_diff[5] = 22.0/525.0; bi_diff[6] = -1.0/40.0;
  ai[1][0] = 0.2;
  ai[2][0] = 3.0/40.0; ai[2][1] = 9.0/40.0;
  ai[3][0] = 44.0/45.0; ai[3][1] = -56.0/15.0; ai[3][2] = 32.0/9.0;
  ai[4][0] = 19372.0/6561.0; ai[4][1] = -25360.0/2187.0;
  ai[4][2] = 64448.0/6561.0; ai[4][3] = -212.0/729.0;
  ai[5][0] = 9017.0/3168.0; ai[5][1] = -355.0/33.0; ai[5][2] = 46732.0/5247.0;
  ai[5][3] = 49.0/176.0; ai[5][4] = -5103.0/18656.0;
  ai[6][0] = 35.0/384.0; ai[6][1] = 0.0; ai[6][2] = 500.0/1113.0;
  ai[6][3] = 125.0/192.0; ai[6][4] = -2187.0/6784.0; ai[6][5] = 11.0/84.0;
  ixx[0][0] = 1500.0/371.0; ixx[0][1] = -1000.0/159.0; ixx[0][2] = 1000.0/371.0;
  ixx[1][0] = -125.0/32.0; ixx[1][1] = 125.0/12.0; ixx[1][2] = -375.0/64.0;
  ixx[2][0] = 9477.0/3392.0; ixx[2][1] = -729.0/106.0; ixx[2][2] = 25515.0/6784.0;
  ixx[3][0] = -11.0/7.0; ixx[3][1] = 11.0/3.0; ixx[3][2] = -55.0/28.0;
  ixx[4][0] = 1.5; ixx[4][1] = -4.0; ixx[4][2] = 2.5;
  i01 = -183.0/64.0; i02=37.0/12.0; i03 = -145.0/128.0;
  /**
     Done setting method parameters.
  */

  t = t_ini;
  //initialise ki
  class_call((*derivs)(t,y_inout,ki,parameters_and_workspace_for_derivs,error_message),error_message,error_message);
  hmin = 100.0*DBL_MIN*fabs(t);
  hmax = fabs(t_final-t_ini)/10.0;
  if (t_vec!=NULL)
    absh = MIN(hmax,fabs(t_vec[1]-t_vec[0]));
  else
    absh = hmax;
  if (absh==0.0)
    absh = hmax;
  //Compute h_initial:
  for (k=0,rh=0.0; k<neq; k++){
    maxtmp = MAX(fabs(y_inout[k]),threshold);
    rh = MAX(rh,fabs(ki[k])/maxtmp);
  }
  rh /= 0.8*pow(rtol,pow_grow);
  if (absh*rh>1.0)
    absh = 1.0/rh;
  if (t_final>t_ini)
    tdir = 1;
  else
    tdir = -1;

  hnew = absh*tdir;
  h = hnew;
  //Find current index:
  if (t_vec != NULL){
    //Output at specified points.
    for(idx = 0; (t_vec[idx]-t)*tdir<0.0; idx++);
  }
  nofailed = _TRUE_;
  
  while ((t-t_final)*tdir<0.0){
    h = hnew;
    hmin = 100.0*DBL_MIN*fabs(t);
    if (fabs(h)<hmin){
      printf("h = %e",fabs(h));
      return _FAILURE_;
    }
    if (fabs(h)>0.9*fabs(t_final-t))
      h = t_final-t;
    //Try to take a step h.
    for (k=0; k<neq; k++){
      //Set ynew and err to the value after i=0 iteration:
      ynew[k] = y_inout[k]+h*bi[0]*ki[k];
      err[k] = h*bi_diff[0]*ki[k];
    }
    for (i=1; i<s; i++){
      //Note: starting at i=1 since method has FSAL property.
      //Reset ytemp.
      for (k=0; k<neq; k++)
  ytemp[k] = y_inout[k];
      for (j=0; j<i; j++){
  for (k=0; k<neq; k++){
    ytemp[k] += h*ai[i][j]*ki[j*neq+k];
  }
      }
      // Calculate k_i:
      if (verbose>2){
  printf("Evaluating ODE at t=%g, ci[i] = %g. h=%g\n",
         t+ci[i]*h,ci[i],h);
  printf("y_inout = [%g,%g]. ytemp = [%g,%g]\n",
         y_inout[0],y_inout[1],ytemp[0],ytemp[1]);
      }
      class_call((*derivs)(t+ci[i]*h,ytemp,ki+i*neq,parameters_and_workspace_for_derivs,error_message),error_message,error_message);
      // Update ynew and err:
      for (k=0; k<neq; k++){
  ynew[k] += h*bi[i]*ki[i*neq+k];
  err[k] += h*bi_diff[i]*ki[i*neq+k];
      }
    }
    if (verbose>3)
      printf("Finished loop over i, new y has been found.\n");
    // Got new y and error estimate.
    for (k=0,errmax = 0.0; k<neq; k++){
      errtemp = fabs(err[k]/MAX(threshold,fabs(ynew[k])));
      if (errtemp>errmax){
  errmax = errtemp;
      }
    }
    if (verbose>3)
      printf("h: %g, errmax = %g\n",h,errmax);
    if (errmax>rtol){
      if (verbose>4)
  printf("Step rejected..\n");
      stats[1]++;
      //Step rejected.
      if (nofailed == _TRUE_){
  hnew = tdir*MAX(hmin, fabs(h) * MAX(0.1, 0.8*pow(rtol/errmax,pow_grow)));
      }
      else{
  hnew = tdir*MAX(hmin,0.5*fabs(h));
  nofailed = _FALSE_;
      }
    }
    else{
      //Step accepted.
      stats[0]++;
      if (print_variables!=NULL){
        class_call((*print_variables)(t+h,ynew,ki+6*neq,
          parameters_and_workspace_for_derivs,error_message),
         error_message,error_message);
      }
      if (verbose>1)
  printf("Step accepted. t=%g, h=%g\n",t,h);
      nofailed = _TRUE_;
      
      hnew = tdir*MAX(hmin, fabs(h) * MAX(0.1, 0.8*pow(rtol/errmax,pow_grow)));
      tnew = t+h;
      //Do we need to write output?
      if (t_vec==NULL){
  //Refined output
  for (idx=1; idx<tres; idx++){
    ti = t+idx/((double) tres)*h;
    ss1 = (ti-t)/h; ss2=ss1*ss1; ss3=ss2*ss1; ss4=ss2*ss2;
    bi_vec_y[0] = ss1+i01*ss2+i02*ss3+i03*ss4;
    bi_vec_dy[0] = 1.0+i01*2.0*ss1+i02*3.0*ss2+i03*4.0*ss3;
    //bi_vec_y[1] = 0.0; bi_vec_dy[1] = 0.0;
    for (i=2; i<7; i++){
      bi_vec_y[i] = ixx[i-2][0]*ss2+ixx[i-2][1]*ss3+ixx[i-2][2]*ss4;
      bi_vec_dy[i] = ixx[i-2][0]*2*ss1+ixx[i-2][1]*3*ss2+ixx[i-2][2]*4*ss3;
    }
    for (k=0; k<neq; k++){
      if (used_in_output[k] == _TRUE_){
        //Interpolate
        yinterp[k] = y_inout[k];
        dyinterp[k] = 0.0;
        for (i=0; i<7; i++){
    if (i!=1){
      yinterp[k] +=h*bi_vec_y[i]*ki[i*neq+k];
      dyinterp[k] +=bi_vec_dy[i]*ki[i*neq+k];
    }
        }
      }
    }
    class_evolver_output((*output)(ti,yinterp,dyinterp,idx,parameters_and_workspace_for_derivs,error_message),
         error_message,error_message);
  }
  class_evolver_output((*output)(tnew,ynew,ki+6*neq,tres,parameters_and_workspace_for_derivs,error_message),
       error_message,error_message);
      }
      else{
  for(; (idx<tres)&&((tnew-t_vec[idx])*tdir>=0.0); idx++){
    if (tnew==t_vec[idx]){
      //We have hit the point exactly. Use ynew and dy=ki+6*neq
      class_evolver_output((*output)(tnew,ynew,ki+6*neq,idx,parameters_and_workspace_for_derivs,error_message),
           error_message,error_message);
    }
    else{
      //Interpolate to get output using the information in the ki-matrix:
      ti = t_vec[idx];
      ss1 = (ti-t)/h; ss2=ss1*ss1; ss3=ss2*ss1; ss4=ss2*ss2;
      bi_vec_y[0] = ss1+i01*ss2+i02*ss3+i03*ss4;
      bi_vec_dy[0] = 1.0+i01*2.0*ss1+i02*3.0*ss2+i03*4.0*ss3;
      //bi_vec_y[1] = 0.0; bi_vec_dy[1] = 0.0;
      for (i=2; i<7; i++){
        bi_vec_y[i] = ixx[i-2][0]*ss2+ixx[i-2][1]*ss3+ixx[i-2][2]*ss4;
        bi_vec_dy[i] = ixx[i-2][0]*2*ss1+ixx[i-2][1]*3*ss2+ixx[i-2][2]*4*ss3;
      }
      for (k=0; k<neq; k++){
        if (used_in_output[k] == _TRUE_){
    //Interpolate
    yinterp[k] = y_inout[k];
    dyinterp[k] = 0.0;
    for (i=0; i<7; i++){
      if (i!=1){
        yinterp[k] +=h*bi_vec_y[i]*ki[i*neq+k];
        dyinterp[k] +=bi_vec_dy[i]*ki[i*neq+k];
      }
    }
        }
      }
      class_evolver_output((*output)(ti,yinterp,dyinterp,idx,parameters_and_workspace_for_derivs,error_message),
           error_message,error_message);
      if (done == _TRUE_) {
        break;
      }
    }
  }
      }
      for (k=0; k<neq; k++){
  //Update y:
  y_inout[k] = ynew[k];
  //Update k0:
  ki[k] = ki[6*neq+k];
      }
      t = tnew;
    }
    if (done == _TRUE_) {
      break;
    }
  }
  if (verbose>0)
    printf(" Successful steps: %d\n Failed steps: %d\n Function evaluations: %d\n",
     stats[0],stats[1],stats[2]);
  free(dy);
  free(ynew);
  free(ytemp);
  free(ki);
  free(err);
  return _SUCCESS_;
}
