
// Modified simple model to include an input for upper bound on b, also
// added a posfun to see how that affects it
DATA_SECTION
  init_int nobs
  init_vector Y(1,nobs)
  init_vector x(1,nobs)
  !! ad_comm::change_datafile_name("bounds.txt");
  init_number upb_b
  !! cout << upb_b << endl;
INITIALIZATION_SECTION
  a 0
  b 0
PARAMETER_SECTION
  init_number a
  init_bounded_number b(-10,upb_b)
  sdreport_number aa
  vector pred_Y(1,nobs)
  number fpen
  number temp
 objective_function_value f
PROCEDURE_SECTION
  fpen=0;
  b=-posfun(-b, -(upb_b-.0001), fpen);
  pred_Y=a*x+b;
 // cout << fpen<< " " <<  b << " " << upb_b-.01 << endl;
 aa=a;
  f=(norm2(pred_Y-Y));
  f=nobs/2.*log(f);    // make it a likelihood function so that
  f+=100000*fpen;                       // covariance matrix is correct
