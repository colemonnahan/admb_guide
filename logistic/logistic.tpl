// Me playing with logistic model as practice. CCM 6/19/2012


DATA_SECTION
  init_int num_years
  init_vector catches(1,num_years)
 !! catches=catches;
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
  !!CLASS ofstream MCMCreport("MCMC.csv",ios::app);
 // !! ad_comm::change_datafile_name("prior_number.dat"); 
  // init_int prior_number
   //!! cout << "prior number is " << prior_number << endl;

INITIALIZATION_SECTION
  logr -3.5
  logK 10.4631

PARAMETER_SECTION
  init_number logr
  init_number logK 
  number K 
  number r
  vector Y_pred(1,num_years)
  vector Y_temp(1,num_years)
  vector Y_final(1,num_years)
  number fpen
  number temp 			// used in intermediate calc
  objective_function_value NLL
  sdreport_number sd_K
  sdreport_number sd_r

PROCEDURE_SECTION
  K=mfexp(logK);
  r=mfexp(logr);
  j=1;
  NLL=0;
  Y_pred[1]=K;
  sd_K=K;
  sd_r=r;

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
 //NLL+=0.5*log(2*3.141593*square(.039))+square(r-.042)/(2*square(.039)); // prior on r is normal (for now)
  // If in MCMC phase, print out the variable values (i.e. this is one iteration)
 if(mceval_phase()) MCMCreport << r << "," << K << endl;

REPORT_SECTION
  cout << "NLL is " << NLL << endl;
  cout << "r is " << r << endl;
  cout << "K is " << K << endl;
  cout << "fpen is " << fpen << endl;
