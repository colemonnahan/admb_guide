
DATA_SECTION
  init_int nobs
  init_vector Y(1,nobs)
  init_vector x(1,nobs)
 !!CLASS ofstream MCMCreport("MCMCreport.csv",ios::app);
PARAMETER_SECTION
 // init_number a(-1)   
  init_number b   
  number a
  sdreport_number aa
  vector pred_Y(1,nobs)
  objective_function_value f
PROCEDURE_SECTION
 a=1.90909098475;
  pred_Y=a*x+b;
 aa=a;
  f=(norm2(pred_Y-Y)); 
  f=nobs/2.*log(f);    // make it a likelihood function so that
  MCMCreport << a <<"," << b << endl;
