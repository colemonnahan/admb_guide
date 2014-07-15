// This is the finance example model modified to have phases read in
DATA_SECTION
  init_int T
  init_vector r(0,T)
  vector sub_r(1,T)
  number h0
  !! ad_comm::change_datafile_name("phases.dat");
  //init_ivector phases(1,4);
  init_int phasea0
  init_int phasea1
  init_int phasea2
  init_int phaseMean
INITIALIZATION_SECTION
  a0 .00016
  a1 .1
  a2 .37
  Mean 0
PARAMETER_SECTION
  init_bounded_number a0(0.0,.01,phasea0)
  init_bounded_number a1(0.0,.5,phasea1)
  init_bounded_number a2(0.0,1.0,phasea2)
  init_bounded_number Mean(-.001,.001,phaseMean)
  vector eps2(1,T)
  vector h(1,T)
  sdreport_number aa
  objective_function_value log_likelihood
PRELIMINARY_CALCS_SECTION
  h0=square(std_dev(r));   // square forms the element-wise square
  sub_r=r(1,T);    // form a subvector so we can use vector operations
  Mean=mean(r);    // calculate the mean of the vector r
PROCEDURE_SECTION
  aa=a0;
  eps2=square(sub_r-Mean);
  h(1)=a0+a2*h0;
  for (int t=2;t<=T;t++)
  {
    h(t)=a0+a1*eps2(t-1)+a2*h(t-1);
  }
  // calculate minus the log-likelihood function and assign it to
  // the object of type objective_function_value
  log_likelihood=.5*sum(log(h)+elem_div(eps2,h));  // elem_div performs
          // element-wise division of vectors
