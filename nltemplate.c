// Note: this template assumes you have 1 solved-for variable. Hence ur[1], dur[3][1]
// If you have N solved-for variables, you need to change to ur[N], dur[3][N]

double cfun(char *label, double *x, double *y, double *z, 
       double ur[1], double dur[3][1], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{
  // your code goes here
}

// ######################################################################

double afun(char *label, double *x, double *y, double *z, 
       double ur[1], double dur[3][1], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{
  // your code goes here
}


// ######################################################################

double ffun(char *label, double *x, double *y, double *z, 
       double ur[1], double dur[3][1], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{
  // your code goes here
}


// ######################################################################

double qfun(char *label, double  *x, double  *y, double  *z, 
                         double *nx, double *ny, double *nz,
       double ur[1], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{
  // your code goes here
}


// ######################################################################

double gfun(char *label, double  *x, double  *y, double  *z, 
                         double *nx, double *ny, double *nz,
       double ur[1], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{
  // your code goes here
}


// ######################################################################

double BCfun(char *label, double *x, double *y, double *z,
       double ur[1], int *nvar, int *istep, 
       int *ireg, int *iregup, double *rnode, double *vec, 
       int *imax, int *jmax, int *kmax, int length)
{
  // your code goes here
}


// ######################################################################

double scanfun(char *label, double  *x, double  *y, double  *z,
                            double *nx, double *ny, double *nz,
       double ur[1], int *nvar, int length)
{
  // your code goes here
}

