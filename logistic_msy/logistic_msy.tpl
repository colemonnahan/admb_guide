// Me playing with logistic model as practice. CCM 6/19/2012
DATA_SECTION
  init_int num_years
  init_vector catches(1,num_years)
  init_int num_obs
  init_vector Y_years(1,num_obs)
  init_vector Y_obs(1,num_obs)
  init_vector Y_SD(1,num_obs)
  init_int temp
  init_int check
  !! cout << check << endl;
  int j
  
  // mceval report; i.e. this tells ADMB where to write the MCMC
  // output in this phase
  !!CLASS ofstream MCMCreport("ADMB.csv",ios::app);
  // !! ad_comm::change_datafile_name("prior_number.dat"); 
  //  init_int prior_number
  //  !! cout << "prior number is " << prior_number << endl;

INITIALIZATION_SECTION
 r .05 
  MSY 480

PARAMETER_SECTION
  init_bounded_number r(0,.5)
  init_bounded_number MSY(.1,1000)
  number K
  vector Y_pred(1,num_years)
  vector Y_final(1,num_years)
  number fpen
  number temp 			// used in intermediate calc
  objective_function_value NLL
  sdreport_number sd_K

PROCEDURE_SECTION
  j=1;
  NLL=0;
  K=4*MSY/r;
  Y_pred[1]=K;
  sd_K=K;
 fpen=0.0;
  for (int i=2; i<=num_years; i++)
  {
  Y_pred(i)=Y_pred(i-1)+Y_pred(i-1)*r*(1-Y_pred(i-1)/K)-catches(i-1);
  Y_pred(i)=posfun(Y_pred(i), 1, fpen);
  if(Y_years(j)==i)
    {
    temp=square(log(Y_obs(j))- log(Y_pred(i)))/(2*square(Y_SD(j)));
    NLL+=temp;
    j++;
    }
  }
 NLL+=100*fpen; 

  // If desired, add the contribution of the priors to get a scaled
  // posterior. Note that not adding a prior implies uniform priors on
  // all parameters.
 // if(prior_number==1)  f+= pow(r-.062,2)/(2*pow(.029,2)); // prior on r is normal if this is one, otherwise uniform so ignore it

  // If in MCMC phase, print out the variable values (i.e. this is one iteration)
 if(mceval_phase()) MCMCreport << r << "," << K << endl;

REPORT_SECTION

 cout << "NLL is " << NLL << endl;
  cout << "r is " << r << endl;
  cout << "K is " << K << endl;
  cout << "fpen is " << fpen << endl;