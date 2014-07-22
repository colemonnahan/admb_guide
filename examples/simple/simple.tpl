
// Modified simple model to include an input for upper bound on b
DATA_SECTION
  init_int nobs
  init_vector Y(1,nobs)
  init_vector x(1,nobs)
  !! ad_comm::change_datafile_name("bounds.txt");
  init_number upb_b
  !! cout << upb_b << endl;
PARAMETER_SECTION
  init_number a
  init_bounded_number b(-10,upb_b)
  sdreport_number aa
  vector pred_Y(1,nobs)
  objective_function_value f
PROCEDURE_SECTION
  pred_Y=a*x+b;
 aa=a;
  f=(norm2(pred_Y-Y));
  f=nobs/2.*log(f);    // make it a likelihood function so that
                       // covariance matrix is correct
