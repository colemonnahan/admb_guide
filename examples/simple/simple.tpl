// a modified version of the simple example packaged with ADMB.
DATA_SECTION
 !!CLASS ofstream leapfrog("leapfrog.csv",ios::trunc);
  init_int nobs
  init_vector Y(1,nobs)
  init_vector x(1,nobs)
PARAMETER_SECTION
  init_number a
  init_number b
  sdreport_number aa
  vector pred_Y(1,nobs)
  objective_function_value f
PROCEDURE_SECTION
  pred_Y=a*x+b;
 aa=a;
  f=(norm2(pred_Y-Y));
  f=nobs/2.*log(f);    // make it a likelihood function so that
                       // covariance matrix is correct
  leapfrog << a << "," << b << "," << f << endl;

REPORT_SECTION
 cout << a << ", " << b << endl;

