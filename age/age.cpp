#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <age.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_MCMCreport = new ofstream("MCMC.csv",ios::app);;
  num_years.allocate("num_years");
  catches.allocate(1,num_years,"catches");
  num_obs.allocate("num_obs");
  Y_years.allocate(1,num_obs,"Y_years");
  Y_obs.allocate(1,num_obs,"Y_obs");
  Y_SD.allocate(1,num_obs,"Y_SD");
  max_age.allocate("max_age");
  z.allocate("z");
 ad_comm::change_datafile_name("age.ctl");
  phase_K.allocate("phase_K");
  phase_r.allocate("phase_r");
  phase_S0.allocate("phase_S0");
  phase_Splus.allocate("phase_Splus");
}

void model_parameters::initializationfunction(void)
{
  K.set_initial_value(5000);
  r.set_initial_value(.05);
  S0.set_initial_value(.8);
  Splus.set_initial_value(.95);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  K.allocate(1000,10000,phase_K,"K");
  r.allocate(0,1,phase_r,"r");
  S0.allocate(0,1,phase_S0,"S0");
  Splus.allocate(0,1,phase_Splus,"Splus");
  age_pred.allocate(1,num_years,1,max_age+1,"age_pred");
  #ifndef NO_AD_INITIALIZE
    age_pred.initialize();
  #endif
  Nplus.allocate(1,num_years,"Nplus");
  #ifndef NO_AD_INITIALIZE
    Nplus.initialize();
  #endif
  Y_final.allocate(1,num_obs,"Y_final");
  #ifndef NO_AD_INITIALIZE
    Y_final.initialize();
  #endif
  fmax.allocate("fmax");
  #ifndef NO_AD_INITIALIZE
  fmax.initialize();
  #endif
  f0.allocate("f0");
  #ifndef NO_AD_INITIALIZE
  f0.initialize();
  #endif
  N0.allocate("N0");
  #ifndef NO_AD_INITIALIZE
  N0.initialize();
  #endif
  NLL.allocate("NLL");
  Kreport.allocate("Kreport");
}

void model_parameters::userfunction(void)
{
  ofstream& MCMCreport= *pad_MCMCreport;
  // Initialize basic params, note that j is just a counter for use in
  // subseting years in which there are abundance estimates. There is
  // probably a much better way to do this. In R it would be
  // Y_final[Y_years] which would then match up to Y_obs.
   j=1;
  NLL=0;
  // Initialize the age structured matrix. Remember that age 0
  // (calves) are in column 1, and likewise, so that all ages are
  // shifted by 1.
   f0=(1-Splus)/(pow(Splus,max_age-1)*S0);
   fmax=f0+r;
   N0=K*pow(Splus,max_age-1)*f0;
   age_pred(1,1)=N0;
   age_pred(1,2)=N0*S0;
   for(int i=3;i<=max_age; i++)  age_pred(1,i)=age_pred(1,i-1)*Splus;
   age_pred(1,max_age+1)=N0/f0;
   // Nplus is the age 1+ animals from the previous year, *not* the
   // plus group!
   // To calculate 1+ group sum the row and subtract off calves from
   // the row sum since they aren't counted in the 1+ group.
   Nplus(1)=sum(row(age_pred,1))-age_pred(1,1);
  // ---------------
   // The first row of the matrix is now ready to go and we can simply
   // loop through each year and make the necessary calculations.
   for(int y=2;y<=num_years;y++)
   {
   dvariable fpen=0.0;
   // the first two age classes still need special calcs, and need to
   // calculate the calves last, by assumption of the model
   age_pred(y,2)=age_pred(y-1,1)*S0; // no fishing pressure on these so no catches
   for(int age=3;age<=max_age+1;age++)
       {
          age_pred(y,age)=(age_pred(y-1,age-1)-
	             catches(y-1)*age_pred(y-1,age-1)/Nplus(y-1))*Splus;
       }
   // Account for the plus group in the last age class by adding the
   // previous plus group
   age_pred(y,max_age+1)+=(age_pred(y-1,max_age+1)-
          catches(y-1)*age_pred(y-1,max_age+1)/Nplus(y-1))*Splus;
   // I assume that calves are born at the end of the year so I need
   // to calculate them last
   Nplus(y)=posfun(sum(row(age_pred,y))-age_pred(y,1), 100, fpen);
   NLL+=10000*fpen;
   if(Nplus(y)<=0) cout << "negative biomass";
   age_pred(y,1)=age_pred(y,max_age+1)*(f0+(fmax-f0)*(1-pow(Nplus(y)/K,z)));
   // within the loop check if this is a year with an abundance
   // estimate and if so add to the NLL
       if(Y_years(j)==y)
         {
	 Y_final(j)=Nplus(y);
	  NLL+=pow(log(Y_obs(j))- log(Y_final(j)),2)/(2*Y_SD(j)*Y_SD(j));
	  j++;
        }
   } // end of loop through years
  // The matrix is now fully complete and the negative loglikelihood
  // is calculated. MLE part is complete.
  // ---------------
  // ---------------
  // Bayesian calculations
  // If desired, add the contribution of the priors to get a scaled
  // posterior. Note that not adding a prior implies uniform priors on
  // all parameters.
  NLL+= pow(r-.042,2)/(2*pow(.019,2)); // prior on r is normal (for now)
  NLL+= pow(S0-.8,2)/(2*pow(.1,2)); // prior on S0 is normal (for now)
  NLL+= pow(Splus-.9,2)/(2*pow(.1,2)); // prior on S0 is normal (for now)
  //cout << Nplus(num_years) << "," << NLL << endl;
  // If in MCMC phase, print out the variable values (i.e. this is one iteration)
 if(mceval_phase())
   {
     MCMCreport << K << "," << r << "," << S0 << "," << Splus << endl;
   }
 // end of Bayesian part
 // ---------------
 // ------------------------------------------------------------
 // ------------------------------------------------------------
 // REPORT_SECTION
 //  report << "Y" << endl << Y_final << endl;
 // ------------------------------------------------------------
 // End of file
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_MCMCreport;
  pad_MCMCreport = NULL;
}

void model_parameters::report(void){}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
  #if defined(__GNUDOS__) || defined(DOS386) || defined(__DPMI32__)  || \
     defined(__MSVC32__)
      if (!arrmblsize) arrmblsize=150000;
  #else
      if (!arrmblsize) arrmblsize=25000;
  #endif
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
