/**
 * module with tools for manipulating arrays
 * Julien Lesgourgues, 18.04.2010
 */

#include "tools_arrays.h"

/**
 * Called by thermodynamics_init(); perturb_sources().
 */
int array_derive(
		 double * array,
		 int n_columns,
		 int n_lines,
		 int index_x,   /** from 0 to (n_columns-1) */
		 int index_y,
		 int index_dydx,
		 char * errmsg) {
  
  int i;

  double dx1,dx2,dy1,dy2,weight1,weight2;

  if ((index_dydx == index_x) || (index_dydx == index_y)) {
    sprintf(errmsg,"%s(L:%d) Output column %d must differ from input columns %d and %d",__func__,__LINE__,index_dydx,index_x,index_y);
    return _FAILURE_;
  }

  dx2=array[1*n_columns+index_x]-array[0*n_columns+index_x];
  dy2=array[1*n_columns+index_y]-array[0*n_columns+index_y];

  for (i=1; i<n_lines-1; i++) {

    dx1 = dx2;
    dy1 = dy2;
    dx2 = array[(i+1)*n_columns+index_x]-array[i*n_columns+index_x];
    dy2 = array[(i+1)*n_columns+index_y]-array[i*n_columns+index_y];
    if ((dx1 == 0) || (dx2 == 0)) {
      sprintf(errmsg,"%s(L:%d) dx1=0 or dx2=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }
    weight1 = dx2*dx2;
    weight2 = dx1*dx1;
    array[i*n_columns+index_dydx] = (weight1*dy1+weight2*dy2) / (weight1*dx1+weight2*dx2);
    
    if (i == 1)
      array[(i-1)*n_columns+index_dydx] = 2.*dy1/dx1 - array[i*n_columns+index_dydx];

    if (i == n_lines-2)
      array[(i+1)*n_columns+index_dydx] = 2.*dy2/dx2 - array[i*n_columns+index_dydx];
  } 
  
  return _SUCCESS_;
}

int array_derive_spline(
		 double * x_array,
		 int n_lines,
		 double * array,
		 double * array_splined,
		 int n_columns,
		 int index_y,
		 int index_dydx,
		 char * errmsg) {
  
  int i;

  double h;

  if (index_dydx == index_y) {
    sprintf(errmsg,"%s(L:%d) Output column %d must differ from input columns %d",__func__,__LINE__,index_dydx,index_y);
    return _FAILURE_;
  }

  for (i=0; i<n_lines-1; i++) {

    h = x_array[i+1] - x_array[i];
    if (h == 0) {
      sprintf(errmsg,"%s(L:%d) h=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    array[i*n_columns+index_dydx] = 
      (array[(i+1)*n_columns+index_y] - array[i*n_columns+index_y])/h
      - h / 6. * (array_splined[(i+1)*n_columns+index_y] + 2. * array_splined[i*n_columns+index_y]);

  } 
  
  array[(n_lines-1)*n_columns+index_dydx] = 
    (array[(n_lines-1)*n_columns+index_y] - array[(n_lines-2)*n_columns+index_y])/h
    + h / 6. * (2. * array_splined[(n_lines-1)*n_columns+index_y] + array_splined[(n_lines-2)*n_columns+index_y]);
  
  return _SUCCESS_;
}

int array_derive_spline_table_line_to_line(
		 double * x_array,
		 int n_lines,
		 double * array,
		 int n_columns,
		 int index_y,
		 int index_ddy,
		 int index_dy,
		 char * errmsg) {
  
  int i;

  double h;

  for (i=0; i<n_lines-1; i++) {

    h = x_array[i+1] - x_array[i];
    if (h == 0) {
      sprintf(errmsg,"%s(L:%d) h=0, stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    array[i*n_columns+index_dy] = 
      (array[(i+1)*n_columns+index_y] - array[i*n_columns+index_y])/h
      - h / 6. * (array[(i+1)*n_columns+index_ddy] + 2. * array[i*n_columns+index_ddy]);

  } 
  
  array[(n_lines-1)*n_columns+index_dy] = 
    (array[(n_lines-1)*n_columns+index_y] - array[(n_lines-2)*n_columns+index_y])/h
    + h / 6. * (2. * array[(n_lines-1)*n_columns+index_ddy] + array[(n_lines-2)*n_columns+index_ddy]);
  
  return _SUCCESS_;
}

int array_derive1_order2_table_line_to_line(
				       double * x_array,
				       int n_lines,
				       double * array,
				       int n_columns,
				       int index_y,
				       int index_dy,
				       char * errmsg) {

  int i;
  double dxp,dxm,dyp,dym;

  i=1;
  dxp = x_array[2] - x_array[1];
  dxm = x_array[0] - x_array[1];
  dyp = *(array+2*n_columns+index_y) - *(array+1*n_columns+index_y);
  dym = *(array+0*n_columns+index_y) - *(array+1*n_columns+index_y);

  if ((dxp*dxm*(dxm-dxp)) == 0.) {
    sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
    return _FAILURE_;
  }

  *(array+1*n_columns+index_dy) = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));

  *(array+0*n_columns+index_dy) = *(array+1*n_columns+index_dy) 
    - (x_array[1] - x_array[0]) * 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm));

  for (i=2; i<n_lines-1; i++) {

    dxp = x_array[i+1] - x_array[i];
    dxm = x_array[i-1] - x_array[i];
    dyp = *(array+(i+1)*n_columns+index_y) - *(array+i*n_columns+index_y);
    dym = *(array+(i-1)*n_columns+index_y) - *(array+i*n_columns+index_y);

    if ((dxp*dxm*(dxm-dxp)) == 0.) {
      sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    *(array+i*n_columns+index_dy) = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));

  }

  *(array+(n_lines-1)*n_columns+index_dy) = *(array+(n_lines-2)*n_columns+index_dy) 
    + (x_array[n_lines-1] - x_array[n_lines-2]) * 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm)); 

  return _SUCCESS_;

}

int array_derive2_order2_table_line_to_line(
				       double * x_array,
				       int n_lines,
				       double * array,
				       int n_columns,
				       int index_y,
				       int index_dy,
				       int index_ddy,
				       char * errmsg) {

  int i;
  double dxp,dxm,dyp,dym;

  for (i=1; i<n_lines-1; i++) {

    dxp = x_array[i+1] - x_array[i];
    dxm = x_array[i-1] - x_array[i];
    dyp = *(array+(i+1)*n_columns+index_y) - *(array+i*n_columns+index_y);
    dym = *(array+(i-1)*n_columns+index_y) - *(array+i*n_columns+index_y);

    if ((dxp*dxm*(dxm-dxp)) == 0.) {
      sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    *(array+i*n_columns+index_dy) = (dyp*dxm*dxm-dym*dxp*dxp)/(dxp*dxm*(dxm-dxp));
    *(array+i*n_columns+index_ddy) = 2.*(dyp*dxm-dym*dxp)/(dxp*dxm*(dxp-dxm));

  }

  *(array+0*n_columns+index_dy) = *(array+1*n_columns+index_dy) 
    - (x_array[1] - x_array[0]) * *(array+1*n_columns+index_ddy);
  *(array+0*n_columns+index_ddy) = *(array+1*n_columns+index_ddy);

  *(array+(n_lines-1)*n_columns+index_dy) = *(array+(n_lines-2)*n_columns+index_dy) 
    + (x_array[n_lines-1] - x_array[n_lines-2]) * *(array+(n_lines-2)*n_columns+index_ddy); 
  *(array+(n_lines-1)*n_columns+index_ddy) = *(array+(n_lines-2)*n_columns+index_ddy);

  return _SUCCESS_;

}

int array_integrate_spline_table_line_to_line(
		 double * x_array,
		 int n_lines,
		 double * array,
		 int n_columns,
		 int index_y,
		 int index_ddy,
		 int index_inty,
		 char * errmsg) {
  
  int i;

  double h;

  *(array+0*n_columns+index_inty)  = 0.; 

  for (i=0; i < n_lines-1; i++) {

    h = (x_array[i+1]-x_array[i]);

    *(array+(i+1)*n_columns+index_inty) = *(array+i*n_columns+index_inty) +
      (array[i*n_columns+index_y]+array[(i+1)*n_columns+index_y])*h/2.+
      (array[i*n_columns+index_ddy]+array[(i+1)*n_columns+index_ddy])*h*h*h/24.;

  }
  
  return _SUCCESS_;
}


 /**
 * Not called.
 */
int array_derive_two(
		     double * array,
		     int n_columns,
		     int n_lines,
		     int index_x,   /** from 0 to (n_columns-1) */
		     int index_y,
		     int index_dydx,
		     int index_ddydxdx,
		     char * errmsg) {
  
  int i;

  double dx1,dx2,dy1,dy2,weight1,weight2;

  if ((index_dydx == index_x) || (index_dydx == index_y)) {
    sprintf(errmsg,"%s(L:%d) : Output column &d must differ from input columns %d and %d",__func__,__LINE__,index_dydx,index_x,index_y);
    return _FAILURE_;
  }

  dx2=*(array+1*n_columns+index_x)-*(array+0*n_columns+index_x);
  dy2=*(array+1*n_columns+index_y)-*(array+0*n_columns+index_y);

  for (i=1; i<n_lines-1; i++) {

    dx1 = dx2;
    dy1 = dy2;
    dx2 = *(array+(i+1)*n_columns+index_x)-*(array+i*n_columns+index_x);
    dy2 = *(array+(i+1)*n_columns+index_y)-*(array+i*n_columns+index_y);
    weight1 = dx2*dx2;
    weight2 = dx1*dx1;

    if ((dx1 == 0.) && (dx2 == 0.)) {
      sprintf(errmsg,"%s(L:%d) stop to avoid division by zero",__func__,__LINE__);
      return _FAILURE_;
    }

    *(array+i*n_columns+index_dydx) = (weight1*dy1+weight2*dy2) / (weight1*dx1+weight2*dx2);
    *(array+i*n_columns+index_ddydxdx) = (dx2*dy1-dx1*dy2) / (weight1*dx1+weight2*dx2);
    
    if (i == 1) {
      *(array+(i-1)*n_columns+index_dydx) = 2.*dy1/dx1 - *(array+i*n_columns+index_dydx);
      *(array+(i-1)*n_columns+index_ddydxdx) = *(array+i*n_columns+index_ddydxdx);
    }

    if (i == n_lines-2) {
      *(array+(i+1)*n_columns+index_dydx) = 2.*dy2/dx2 - *(array+i*n_columns+index_dydx);
      *(array+(i+1)*n_columns+index_dydx) = *(array+i*n_columns+index_ddydxdx);
    }
  } 
  
  return _SUCCESS_;
}

int array_spline(
		  double * array,
		  int n_columns,
		  int n_lines,
		  int index_x,   /** from 0 to (n_columns-1) */
		  int index_y,
		  int index_ddydx2,
		  short spline_mode,
		  char * errmsg) {

  int i,k;
  double p,qn,sig,un;
  double * u;
  double dy_first;
  double dy_last;

  if (n_lines < 3) {
    sprintf(errmsg,"%s(L:%d) n_lines=%d, while routine needs n_lines >= 3",__func__,__LINE__,n_lines);
    return _FAILURE_;
  }

  u = malloc((n_lines-1) * sizeof(double));
  if (u == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate u",__func__,__LINE__);
    return _FAILURE_;
  }

  if (spline_mode == _SPLINE_NATURAL_) {
    *(array+0*n_columns+index_ddydx2) = u[0] = 0.0;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_first = 
	((*(array+2*n_columns+index_x)-*(array+0*n_columns+index_x))*
	 (*(array+2*n_columns+index_x)-*(array+0*n_columns+index_x))*
	 (*(array+1*n_columns+index_y)-*(array+0*n_columns+index_y))-
	 (*(array+1*n_columns+index_x)-*(array+0*n_columns+index_x))*
	 (*(array+1*n_columns+index_x)-*(array+0*n_columns+index_x))*
	 (*(array+2*n_columns+index_y)-*(array+0*n_columns+index_y)))/
	((*(array+2*n_columns+index_x)-*(array+0*n_columns+index_x))*
	 (*(array+1*n_columns+index_x)-*(array+0*n_columns+index_x))*
	 (*(array+2*n_columns+index_x)-*(array+1*n_columns+index_x)));

      *(array+0*n_columns+index_ddydx2) = -0.5;

      u[0] = 
	(3./(*(array+1*n_columns+index_x) -  *(array+0*n_columns+index_x)))*
	((*(array+1*n_columns+index_y) -  *(array+0*n_columns+index_y))/
	 (*(array+1*n_columns+index_x) -  *(array+0*n_columns+index_x))
	 -dy_first);
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }

  for (i=1; i < n_lines-1; i++) {

      sig = (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x))
	/ (*(array+(i+1)*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));

      p = sig * *(array+(i-1)*n_columns+index_ddydx2) + 2.0;
      *(array+i*n_columns+index_ddydx2) = (sig-1.0)/p;
      u[i] = (*(array+(i+1)*n_columns+index_y) - *(array+i*n_columns+index_y))
	/ (*(array+(i+1)*n_columns+index_x) - *(array+i*n_columns+index_x))
	- (*(array+i*n_columns+index_y) - *(array+(i-1)*n_columns+index_y))
	/ (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));
      u[i]= (6.0 * u[i] / 
	     (*(array+(i+1)*n_columns+index_x) - *(array+(i-1)*n_columns+index_x))
	     - sig * u[i-1]) / p;

    }

  if (spline_mode == _SPLINE_NATURAL_) {
    qn=0.;
    un=0.;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_last = 
	((*(array+(n_lines-3)*n_columns+index_x)-*(array+(n_lines-1)*n_columns+index_x))*
	 (*(array+(n_lines-3)*n_columns+index_x)-*(array+(n_lines-1)*n_columns+index_x))*
	 (*(array+(n_lines-2)*n_columns+index_y)-*(array+(n_lines-1)*n_columns+index_y))-
	 (*(array+(n_lines-2)*n_columns+index_x)-*(array+(n_lines-1)*n_columns+index_x))*
	 (*(array+(n_lines-2)*n_columns+index_x)-*(array+(n_lines-1)*n_columns+index_x))*
	 (*(array+(n_lines-3)*n_columns+index_y)-*(array+(n_lines-1)*n_columns+index_y)))/
	((*(array+(n_lines-3)*n_columns+index_x)-*(array+(n_lines-1)*n_columns+index_x))*
	 (*(array+(n_lines-2)*n_columns+index_x)-*(array+(n_lines-1)*n_columns+index_x))*
	 (*(array+(n_lines-3)*n_columns+index_x)-*(array+(n_lines-2)*n_columns+index_x)));

      qn=0.5;
      un =
	(3./(*(array+(n_lines-1)*n_columns+index_x) -  *(array+(n_lines-2)*n_columns+index_x)))*
	(dy_last-(*(array+(n_lines-1)*n_columns+index_y) -  *(array+(n_lines-2)*n_columns+index_y))/
	 (*(array+(n_lines-1)*n_columns+index_x) -  *(array+(n_lines-2)*n_columns+index_x)));
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }

  *(array+(n_lines-1)*n_columns+index_ddydx2) = 
    (un-qn*u[n_lines-2])/(qn* *(array+(n_lines-2)*n_columns+index_ddydx2)+1.0);

  for (k=n_lines-2; k>=0; k--)
    *(array+k*n_columns+index_ddydx2) = *(array+k*n_columns+index_ddydx2) *
      *(array+(k+1)*n_columns+index_ddydx2) + u[k];
  
  free(u);

  return _SUCCESS_;
 }

int array_spline_table_line_to_line(
				    double * x, /* vector of size x_size */
				    int n_lines,
				    double * array,
				    int n_columns,
				    int index_y,
				    int index_ddydx2,
				    short spline_mode,
				    char * errmsg) {
  
  int i,k;
  double p,qn,sig,un;
  double * u;
  double dy_first;
  double dy_last;

  u = malloc((n_lines-1) * sizeof(double));
  if (u == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate u",__func__,__LINE__);
    return _FAILURE_;
  }

  if (spline_mode == _SPLINE_NATURAL_) {
    *(array+0*n_columns+index_ddydx2) = u[0] = 0.0;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_first = 
	((x[2]-x[0])*(x[2]-x[0])*
	 (*(array+1*n_columns+index_y)-*(array+0*n_columns+index_y))-
	 (x[1]-x[0])*(x[1]-x[0])*
	 (*(array+2*n_columns+index_y)-*(array+0*n_columns+index_y)))/
	((x[2]-x[0])*(x[1]-x[0])*(x[2]-x[1]));
      *(array+0*n_columns+index_ddydx2) = -0.5;
      u[0] = 
	(3./(x[1] -  x[0]))*
	((*(array+1*n_columns+index_y) -  *(array+0*n_columns+index_y))/
	 (x[1] - x[0])-dy_first);
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }

  for (i=1; i < n_lines-1; i++) {

      sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);

      p = sig * *(array+(i-1)*n_columns+index_ddydx2) + 2.0;
      *(array+i*n_columns+index_ddydx2) = (sig-1.0)/p;
      u[i] = (*(array+(i+1)*n_columns+index_y) - *(array+i*n_columns+index_y))
	/ (x[i+1] - x[i])
	- (*(array+i*n_columns+index_y) - *(array+(i-1)*n_columns+index_y))
	/ (x[i] - x[i-1]);
      u[i]= (6.0 * u[i] / 
	     (x[i+1] - x[i-1])
	     - sig * u[i-1]) / p;

  }

  if (spline_mode == _SPLINE_NATURAL_) {
    qn=0.;
    un=0.;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_last = 
	((x[n_lines-3]-x[n_lines-1])*(x[n_lines-3]-x[n_lines-1])*
	 (*(array+(n_lines-2)*n_columns+index_y)-*(array+(n_lines-1)*n_columns+index_y))-
	 (x[n_lines-2]-x[n_lines-1])*(x[n_lines-2]-x[n_lines-1])*
	 (*(array+(n_lines-3)*n_columns+index_y)-*(array+(n_lines-1)*n_columns+index_y)))/
	((x[n_lines-3]-x[n_lines-1])*(x[n_lines-2]-x[n_lines-1])*(x[n_lines-3]-x[n_lines-2]));
      qn=0.5;
      un =
	(3./(x[n_lines-1] - x[n_lines-2]))*
	(dy_last-(*(array+(n_lines-1)*n_columns+index_y) -  *(array+(n_lines-2)*n_columns+index_y))/
	 (x[n_lines-1] - x[n_lines-2]));
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }

  *(array+(n_lines-1)*n_columns+index_ddydx2) = 
    (un-qn*u[n_lines-2])/(qn* *(array+(n_lines-2)*n_columns+index_ddydx2)+1.0);

  for (k=n_lines-2; k>=0; k--)
    *(array+k*n_columns+index_ddydx2) = *(array+k*n_columns+index_ddydx2) *
      *(array+(k+1)*n_columns+index_ddydx2) + u[k];
  
  free(u);

  return _SUCCESS_;
 }

int array_spline_table_lines(
			     double * x, /* vector of size x_size */
			     int x_size,
			     double * y_array, /* array of size x_size*y_size with elements 
						  y_array[index_x*y_size+index_y] */
			     int y_size,   
			     double * ddy_array, /* array of size x_size*y_size */
			     short spline_mode,
			     char * errmsg
			     ) {

  double * p;
  double * qn;
  double * un; 
  double * u;
  double sig;
  int index_x;
  int index_y;
  double dy_first;
  double dy_last;

  u = malloc((x_size-1) * y_size * sizeof(double));
  p = malloc(y_size * sizeof(double));
  qn = malloc(y_size * sizeof(double));
  un = malloc(y_size * sizeof(double));
  if (u == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate u",__func__,__LINE__);
    return _FAILURE_;
  }
  if (p == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate p",__func__,__LINE__);
    return _FAILURE_;
  }
  if (qn == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate qn",__func__,__LINE__);
    return _FAILURE_;
  }
  if (un == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate un",__func__,__LINE__);
    return _FAILURE_;
  }
  

  index_x=0;

  if (spline_mode == _SPLINE_NATURAL_) {
    for (index_y=0; index_y < y_size; index_y++) {
      ddy_array[index_x*y_size+index_y] = u[index_x*y_size+index_y] = 0.0;
    }
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {

      for (index_y=0; index_y < y_size; index_y++) {

	dy_first = 
	  ((x[2]-x[0])*(x[2]-x[0])*
	   (y_array[1*y_size+index_y]-y_array[0*y_size+index_y])-
	   (x[1]-x[0])*(x[1]-x[0])*
	   (y_array[2*y_size+index_y]-y_array[0*y_size+index_y]))/
	  ((x[2]-x[0])*(x[1]-x[0])*(x[2]-x[1]));
	
	ddy_array[index_x*y_size+index_y] = -0.5;
	
	u[index_x*y_size+index_y] =
	  (3./(x[1] -  x[0]))*
	  ((y_array[1*y_size+index_y]-y_array[0*y_size+index_y])/
	   (x[1] - x[0])-dy_first);
	
      }
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }
    

  for (index_x=1; index_x < x_size-1; index_x++) {

    sig = (x[index_x] - x[index_x-1])/(x[index_x+1] - x[index_x-1]);

    for (index_y=0; index_y < y_size; index_y++) {

      p[index_y] = sig * ddy_array[(index_x-1)*y_size+index_y] + 2.0;

      ddy_array[index_x*y_size+index_y] = (sig-1.0)/p[index_y];

      u[index_x*y_size+index_y] = 
	(y_array[(index_x+1)*y_size+index_y] - y_array[index_x*y_size+index_y])
	/ (x[index_x+1] - x[index_x])
	- (y_array[index_x*y_size+index_y] - y_array[(index_x-1)*y_size+index_y])
	/ (x[index_x] - x[index_x-1]);

      u[index_x*y_size+index_y] = (6.0 * u[index_x*y_size+index_y] /
				   (x[index_x+1] - x[index_x-1]) 
				   - sig * u[(index_x-1)*y_size+index_y]) / p[index_y];
    }

  }

  if (spline_mode == _SPLINE_NATURAL_) {

    for (index_y=0; index_y < y_size; index_y++) {
      qn[index_y]=un[index_y]=0.0;
    }

  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {

      for (index_y=0; index_y < y_size; index_y++) {

	dy_last = 
	  ((x[x_size-3]-x[x_size-1])*(x[x_size-3]-x[x_size-1])*
	   (y_array[(x_size-2)*y_size+index_y]-y_array[(x_size-1)*y_size+index_y])-
	   (x[x_size-2]-x[x_size-1])*(x[x_size-2]-x[x_size-1])*
	   (y_array[(x_size-3)*y_size+index_y]-y_array[(x_size-1)*y_size+index_y]))/
	  ((x[x_size-3]-x[x_size-1])*(x[x_size-2]-x[x_size-1])*(x[x_size-3]-x[x_size-2]));

	qn[index_y]=0.5;

	un[index_y]=
	  (3./(x[x_size-1] - x[x_size-2]))*
	  (dy_last-(y_array[(x_size-1)*y_size+index_y] - y_array[(x_size-2)*y_size+index_y])/
	   (x[x_size-1] - x[x_size-2]));	

      }
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }
    
  index_x=x_size-1;


  for (index_y=0; index_y < y_size; index_y++) {
    ddy_array[index_x*y_size+index_y] = 
      (un[index_y] - qn[index_y] * u[(index_x-1)*y_size+index_y]) /
      (qn[index_y] * ddy_array[(index_x-1)*y_size+index_y] + 1.0);
  }

  for (index_x=x_size-2; index_x >= 0; index_x--) {
    for (index_y=0; index_y < y_size; index_y++) {

      ddy_array[index_x*y_size+index_y] = ddy_array[index_x*y_size+index_y] *
	ddy_array[(index_x+1)*y_size+index_y] + u[index_x*y_size+index_y];

    }
  }

  free(qn);
  free(un);
  free(p);
  free(u);

  return _SUCCESS_;
 }

int array_spline_table_columns(
		       double * x, /* vector of size x_size */
		       int x_size,
		       double * y_array, /* array of size x_size*y_size with elements 
					  y_array[index_y*x_size+index_x] */
		       int y_size,    
		       double * ddy_array, /* array of size x_size*y_size */
		       short spline_mode,
		       char * errmsg
		       ) {

  double * p;
  double * qn;
  double * un; 
  double * u;
  double sig;
  int index_x;
  int index_y;
  double dy_first;
  double dy_last;

  u = malloc((x_size-1) * y_size * sizeof(double));
  p = malloc(y_size * sizeof(double));
  qn = malloc(y_size * sizeof(double));
  un = malloc(y_size * sizeof(double));
  if (u == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate u",__func__,__LINE__);
    return _FAILURE_;
  }
  if (p == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate p",__func__,__LINE__);
    return _FAILURE_;
  }
  if (qn == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate qn",__func__,__LINE__);
    return _FAILURE_;
  }
  if (un == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate un",__func__,__LINE__);
    return _FAILURE_;
  }

  index_x=0;

  if (spline_mode == _SPLINE_NATURAL_) {
    for (index_y=0; index_y < y_size; index_y++) {
      ddy_array[index_y*x_size+index_x] = 0.0;
      u[index_x*y_size+index_y] = 0.0;
    }
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {

      for (index_y=0; index_y < y_size; index_y++) {

	dy_first = 
	  ((x[2]-x[0])*(x[2]-x[0])*
	   (y_array[index_y*x_size+1]-y_array[index_y*x_size+0])-
	   (x[1]-x[0])*(x[1]-x[0])*
	   (y_array[index_y*x_size+2]-y_array[index_y*x_size+0]))/
	  ((x[2]-x[0])*(x[1]-x[0])*(x[2]-x[1]));
	
	ddy_array[index_y*x_size+index_x] = -0.5;
	
	u[index_x*y_size+index_y] =
	  (3./(x[1] -  x[0]))*
	  ((y_array[index_y*x_size+1]-y_array[index_y*x_size+0])/
	   (x[1] - x[0])-dy_first);
	
      }
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }
  
  for (index_x=1; index_x < x_size-1; index_x++) {

    sig = (x[index_x] - x[index_x-1])/(x[index_x+1] - x[index_x-1]);

    for (index_y=0; index_y < y_size; index_y++) {

      p[index_y] = sig * ddy_array[index_y*x_size+(index_x-1)] + 2.0;

      ddy_array[index_y*x_size+index_x] = (sig-1.0)/p[index_y];

      u[index_x*y_size+index_y] =	
	(y_array[index_y*x_size+(index_x+1)] - y_array[index_y*x_size+index_x])
	/ (x[index_x+1] - x[index_x])
	- (y_array[index_y*x_size+index_x] - y_array[index_y*x_size+(index_x-1)])
	/ (x[index_x] - x[index_x-1]);

      u[index_x*y_size+index_y] = (6.0 * u[index_x*y_size+index_y] /
				   (x[index_x+1] - x[index_x-1])
				   - sig * u[(index_x-1)*y_size+index_y]) / p[index_y];
    }

  }

  if (spline_mode == _SPLINE_NATURAL_) {

    for (index_y=0; index_y < y_size; index_y++) {
      qn[index_y]=un[index_y]=0.0;
    }

  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {

      for (index_y=0; index_y < y_size; index_y++) {

	dy_last = 
	  ((x[x_size-3]-x[x_size-1])*(x[x_size-3]-x[x_size-1])*
	   (y_array[index_y*x_size+(x_size-2)]-y_array[index_y*x_size+(x_size-1)])-
	   (x[x_size-2]-x[x_size-1])*(x[x_size-2]-x[x_size-1])*
	   (y_array[index_y*x_size+(x_size-3)]-y_array[index_y*x_size+(x_size-1)]))/
	  ((x[x_size-3]-x[x_size-1])*(x[x_size-2]-x[x_size-1])*(x[x_size-3]-x[x_size-2]));

	qn[index_y]=0.5;

	un[index_y]=
	  (3./(x[x_size-1] - x[x_size-2]))*
	  (dy_last-(y_array[index_y*x_size+(x_size-1)] - y_array[index_y*x_size+(x_size-2)])/
	   (x[x_size-1] - x[x_size-2]));	

      }
    }
    else {
      sprintf(errmsg,"%s(L:%d) Spline mode not identified: %d",__func__,__LINE__,spline_mode);
      return _FAILURE_;
    }
  }
    
  index_x=x_size-1;

  for (index_y=0; index_y < y_size; index_y++) {
    ddy_array[index_y*x_size+index_x] =
      (un[index_y] - qn[index_y] * u[(index_x-1)*y_size+index_y]) /
      (qn[index_y] * ddy_array[index_y*x_size+(index_x-1)] + 1.0);
  }

  for (index_x=x_size-2; index_x >= 0; index_x--) {
    for (index_y=0; index_y < y_size; index_y++) {

      ddy_array[index_y*x_size+index_x] = ddy_array[index_y*x_size+index_x] *
	ddy_array[index_y*x_size+(index_x+1)] + u[index_x*y_size+index_y];

    }
  }

  free(qn);
  free(p);
  free(u);
  free(un);

  return _SUCCESS_;
 }

int array_integrate_all_spline(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_ddy,
		   double * result,
		   char * errmsg) {

  int i;
  double h;

  *result = 0; 

  for (i=0; i < n_lines-1; i++) {

    h = (array[(i+1)*n_columns+index_x]-array[i*n_columns+index_x]);

    *result += 
      (array[i*n_columns+index_y]+array[(i+1)*n_columns+index_y])*h/2.+
      (array[i*n_columns+index_ddy]+array[(i+1)*n_columns+index_ddy])*h*h*h/24.;

  }

  return _SUCCESS_;
}

 /**
 * Not called.
 */
int array_integrate(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   int index_int_y_dx,
		   char * errmsg) {

  int i;
  double sum;

  if ((index_int_y_dx == index_x) || (index_int_y_dx == index_y)) {
    sprintf(errmsg,"%s(L:%d) : Output column &d must differ from input columns %d and %d",__func__,__LINE__,index_int_y_dx,index_x,index_y);
    return _FAILURE_;
  }

  sum=0.;
  *(array+0*n_columns+index_int_y_dx)=sum;

  for (i=1; i<n_lines; i++) {

    sum += 0.5 * (*(array+i*n_columns+index_y) + *(array+(i-1)*n_columns+index_y))
               * (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));

    *(array+i*n_columns+index_int_y_dx)=sum;  
  }
 

  return _SUCCESS_;
}

 /**
 * Called by thermodynamics_init().
 */
int array_integrate_ratio(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y1,
		   int index_y2,
		   int index_int_y1_over_y2_dx,
		   char * errmsg) {

  int i;
  double sum;

  if ((index_int_y1_over_y2_dx == index_x) || (index_int_y1_over_y2_dx == index_y1) || (index_int_y1_over_y2_dx == index_y2)) {
    sprintf(errmsg,"%s(L:%d) : Output column &d must differ from input columns %d, %d and %d",__func__,__LINE__,index_int_y1_over_y2_dx,index_x,index_y1,index_y2);
    return _FAILURE_;
  }

  sum=0.;

  *(array+0*n_columns+index_int_y1_over_y2_dx)=sum;

  for (i=1; i<n_lines; i++) {

    sum += 0.5 * (*(array+i*n_columns+index_y1) / *(array+i*n_columns+index_y2)
		  + *(array+(i-1)*n_columns+index_y1) / *(array+(i-1)*n_columns+index_y2))
      * (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));

    *(array+i*n_columns+index_int_y1_over_y2_dx)=sum;  
  }
 

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   double x,
		   int * last_index,
		   double * result,
		   int result_size, /** from 1 to n_columns */
		   char * errmsg) {

  int inf,sup,mid,i;
  double weight;

  inf=0;
  sup=n_lines-1;

  if (*(array+inf*n_columns+index_x) < *(array+sup*n_columns+index_x)){

    if (x < *(array+inf*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array+inf*n_columns+index_x));
      return _FAILURE_;
    }

    if (x > *(array+sup*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array+sup*n_columns+index_x));
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < *(array+mid*n_columns+index_x)) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < *(array+sup*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array+sup*n_columns+index_x));
      return _FAILURE_;
    }

    if (x > *(array+inf*n_columns+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array+inf*n_columns+index_x));
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > *(array+mid*n_columns+index_x)) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  weight=(x-*(array+inf*n_columns+index_x))/(*(array+sup*n_columns+index_x)-*(array+inf*n_columns+index_x));

  for (i=0; i<result_size; i++)
    *(result+i) = *(array+inf*n_columns+i) * (1.-weight)
      + weight * *(array+sup*n_columns+i);

  *(result+index_x) = x;

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline(
			     double * x_array,
			     int n_lines,
			     double * array,
			     double * array_splined,
			     int n_columns,
			     double x,
			     int * last_index,
			     double * result,
			     int result_size, /** from 1 to n_columns */
			     char * errmsg) {

  int inf,sup,mid,i;
  double h,a,b;
  
  inf=0;
  sup=n_lines-1;
  
  if (x_array[inf] < x_array[sup]){

    if (x < x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    if (x > x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < x_array[sup]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,x_array[sup]);
      return _FAILURE_;
    }

    if (x > x_array[inf]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,x_array[inf]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > x_array[mid]) {sup=mid;}
      else {inf=mid;}

    }

  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) = 
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i) +
      ((a*a*a-a)* *(array_splined+inf*n_columns+i) + 
       (b*b*b-b)* *(array_splined+sup*n_columns+i))*h*h/6.;

  return _SUCCESS_;
}


 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_growing_closeby(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   double x,
		   int * last_index,
		   double * result,
		   int result_size, /** from 1 to n_columns */
		   char * errmsg) {

  int inf,sup,mid,i;
  double weight;

  inf = *last_index;
  sup = *last_index+1;

  while (x < *(array+inf*n_columns+index_x)) {
    inf--;
    if (inf < 0) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,array[index_x]);
      return _FAILURE_;
    }
  }
  sup = inf+1;
  while (x > *(array+sup*n_columns+index_x)) {
    sup++;
    if (sup > (n_lines-1)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,array[(n_lines-1)*n_columns+index_x]);
      return _FAILURE_;
    }
  }
  inf = sup-1;
  
  *last_index = inf;

  weight=(x-*(array+inf*n_columns+index_x))/(*(array+sup*n_columns+index_x)-*(array+inf*n_columns+index_x));

  for (i=0; i<result_size; i++)
    *(result+i) = *(array+inf*n_columns+i) * (1.-weight)
      + weight * *(array+sup*n_columns+i);
  
  *(result+index_x) = x;

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably very close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline_growing_closeby(
					     double * x_array,
					     int n_lines,
					     double * array,
					     double * array_splined,
					     int n_columns,
					     double x,
					     int * last_index,
					     double * result,
					     int result_size, /** from 1 to n_columns */
					     char * errmsg) {

  int inf,sup,mid,i;
  double h,a,b;

  inf = *last_index;
  while (x < x_array[inf]) {
    inf--;
    if (inf < 0) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,x_array[0]);
      return _FAILURE_;
    }
  }
  sup = inf+1;
  while (x > x_array[sup]) {
    sup++;
    if (sup > (n_lines-1)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,x_array[n_lines-1]);
      return _FAILURE_;
    }
  }
  inf = sup-1;
  
  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) = 
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i) +
      ((a*a*a-a)* *(array_splined+inf*n_columns+i) + 
       (b*b*b-b)* *(array_splined+sup*n_columns+i))*h*h/6.;

  return _SUCCESS_;
}

 /**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably close (but maybe not so close) to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
int array_interpolate_spline_growing_hunt(
					     double * x_array,
					     int n_lines,
					     double * array,
					     double * array_splined,
					     int n_columns,
					     double x,
					     int * last_index,
					     double * result,
					     int result_size, /** from 1 to n_columns */
					     char * errmsg) {

  int inf,sup,mid,i,inc;
  double h,a,b;

  inc=1;

  if (x >= x_array[*last_index]) {
    if (x > x_array[n_lines-1]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,
	      x,x_array[n_lines-1]);
      return _FAILURE_;
    }
    /* try closest neighboor upward */
    inf = *last_index;
    sup = inf + inc;
    if (x > x_array[sup]) {
      /* hunt upward */
      while (x > x_array[sup]) {
	inf = sup;
	inc += 1;
	sup += inc;
	if (sup > n_lines-1) {
	  sup = n_lines-1;
	}
      }
      /* bisect */
      while (sup-inf > 1) {
	mid=(int)(0.5*(inf+sup));
	if (x < x_array[mid]) {sup=mid;}
	else {inf=mid;}
      }
    }
   }
  else {
    if (x < x_array[0]) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,
	      x,x_array[0]);
      return _FAILURE_;
    }
    /* try closest neighboor downward */
    sup = *last_index;
    inf = sup - inc;
    if (x < x_array[inf]) {
      /* hunt downward */
      while (x < x_array[inf]) {
	sup = inf;
	inc += 1;
	inf -= inc;
	if (inf < 0) {
	  inf = 0; 
	}
      }
      /* bisect */
      while (sup-inf > 1) {
	mid=(int)(0.5*(inf+sup));
	if (x < x_array[mid]) {sup=mid;}
	else {inf=mid;}
      }
    }
  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x-x_array[inf])/h;
  a = 1-b;

  for (i=0; i<result_size; i++)
    *(result+i) = 
      a * *(array+inf*n_columns+i) +
      b * *(array+sup*n_columns+i) +
      ((a*a*a-a)* *(array_splined+inf*n_columns+i) + 
       (b*b*b-b)* *(array_splined+sup*n_columns+i))*h*h/6.;

  return _SUCCESS_;
}

/** 
 * interpolate linearily to get y_i(x), when x and y_i are in two different arrays
 *
 * Called by transfer_interpolate_sources(); transfer_functions_at_k(); perturb_sources_at_eta().
 */
int array_interpolate_two(
		   double * array_x,
		   int n_columns_x,
		   int index_x,   /** from 0 to (n_columns_x-1) */
		   double * array_y,
		   int n_columns_y,
		   int n_lines,  /** must be the same for array_x and array_y */
		   double x,
		   double * result,
		   int result_size, /** from 1 to n_columns_y */
		   char * errmsg) {

  int inf,sup,mid,i;
  double weight;

  inf=0;
  sup=n_lines-1;

  if (array_x[inf*n_columns_x+index_x] < array_x[sup*n_columns_x+index_x]){

    if (x < array_x[inf*n_columns_x+index_x]) {

      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,array_x[inf*n_columns_x+index_x]);
      return _FAILURE_;
    }

    if (x > array_x[sup*n_columns_x+index_x]) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,array_x[sup*n_columns_x+index_x]);
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x < array_x[mid*n_columns_x+index_x]) {sup=mid;}
      else {inf=mid;}

    }

  }

  else {

    if (x < *(array_x+sup*n_columns_x+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e < x_min=%e",__func__,__LINE__,x,*(array_x+sup*n_columns_x+index_x));
      return _FAILURE_;
    }

    if (x > *(array_x+inf*n_columns_x+index_x)) {
      sprintf(errmsg,"%s(L:%d) : x=%e > x_max=%e",__func__,__LINE__,x,*(array_x+inf*n_columns_x+index_x));
      return _FAILURE_;
    }

    while (sup-inf > 1) {

      mid=(int)(0.5*(inf+sup));
      if (x > *(array_x+mid*n_columns_x+index_x)) {sup=mid;}
      else {inf=mid;}

    }

  }

  weight=(x-*(array_x+inf*n_columns_x+index_x))/(*(array_x+sup*n_columns_x+index_x)-*(array_x+inf*n_columns_x+index_x));

  for (i=0; i<result_size; i++)
    *(result+i) = *(array_y+i*n_lines+inf) * (1.-weight)
      + weight * *(array_y+i*n_lines+sup) ;

/**     *(result+i) = *(array_y+inf*n_columns_y+i) * (1.-weight) */
/**       + weight * *(array_y+sup*n_columns_y+i) ; */

  return _SUCCESS_;
}


/** 
 * Called by transfer_solve().
 */
int array_interpolate_equal(
			    double * array,
			    int n_columns,
			    int n_lines,
			    double x,
			    double x_min,
			    double x_max,
			    double * result,
			    char * errmsg) {
  
  int index_minus,i;
  double x_step,x_minus,weight;

  if (x < x_min) {
    sprintf(errmsg,"%s(L:%d) : x out of bounds: x=%e,x_min=%e",__func__,__LINE__,x,x_min);
    return _FAILURE_;
  }

  if (x > x_max) {
    sprintf(errmsg,"%s(L:%d) : x out of bounds: x=%e,x_max=%e",__func__,__LINE__,x,x_max);
    return _FAILURE_;
  }

  x_step = (x_max-x_min)/(n_lines-1);
  index_minus = (int)((x-x_min)/x_step);
  x_minus = index_minus * x_step;
  weight = (x-x_minus) / x_step;
 
  for (i=0; i<n_columns; i++) 
    result[i] = *(array+n_columns*index_minus+i)*(1.-weight)
      + *(array+n_columns*(index_minus+1)+i)*weight;

  return _SUCCESS_;

}

/** 
 * Called by transfer_solve().
 */
int array_integrate_all(
		   double * array,
		   int n_columns,
		   int n_lines,
		   int index_x,   /** from 0 to (n_columns-1) */
		   int index_y,
		   double *result) {

  int i;
  double sum;

  sum=0.;

  for (i=1; i<n_lines; i++) {

    sum += 0.5 * (*(array+i*n_columns+index_y) + *(array+(i-1)*n_columns+index_y))
               * (*(array+i*n_columns+index_x) - *(array+(i-1)*n_columns+index_x));

  }

  *result = sum;

  return _SUCCESS_;

}

int array_smooth(double * array,
		 int n_columns,
		 int n_lines,
		 int index, /** from 0 to (n_columns-1) */
		 int radius,
		 char * errmsg) {

  double * smooth;
  int i,j,jmin,jmax;
  double weigth;

  smooth=malloc(n_lines*sizeof(double));
  if (smooth == NULL) {
    sprintf(errmsg,"%s(L:%d) Cannot allocate smooth",__func__,__LINE__);
    return _FAILURE_;
  }

  for (i=0; i<n_lines; i++) {
    smooth[i]=0.;
    weigth=0.;
    jmin = max(i-radius,0);
    jmax = min(i+radius,n_lines-1);
    for (j=jmin; j <= jmax; j++) {
      smooth[i] += array[j*n_columns+index];
      weigth += 1.;
    }
    smooth[i] /= weigth;
  }

  for (i=0; i<n_lines; i++)
    array[i*n_columns+index] = smooth[i];

  free(smooth);

  return _SUCCESS_;

}