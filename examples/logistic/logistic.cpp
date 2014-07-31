#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <logistic.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  num_years.allocate("num_years");
  catches.allocate(1,num_years,"catches");
  num_obs.allocate("num_obs");
  Y_years.allocate(1,num_obs,"Y_years");
  Y_obs.allocate(1,num_obs,"Y_obs");
  Y_SD.allocate(1,num_obs,"Y_SD");
  temp.allocate("temp");
  check.allocate("check");
 cout << check << endl;
  pad_MCMCreport = new ofstream("ADMB.csv",ios::app);;
}

void model_parameters::initializationfunction(void)
{
  r.set_initial_value(0.025);
  K.set_initial_value(5000);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  r.allocate(0,.15,"r");
  K.allocate("K");
  Y_pred.allocate(1,num_years,"Y_pred");
  #ifndef NO_AD_INITIALIZE
    Y_pred.initialize();
  #endif
  Y_temp.allocate(1,num_years,"Y_temp");
  #ifndef NO_AD_INITIALIZE
    Y_temp.initialize();
  #endif
  Y_final.allocate(1,num_years,"Y_final");
  #ifndef NO_AD_INITIALIZE
    Y_final.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  temp.allocate("temp");
  #ifndef NO_AD_INITIALIZE
  temp.initialize();
  #endif
  NLL.allocate("NLL");
  sd_K.allocate("sd_K");
}

void model_parameters::userfunction(void)
{
  ofstream& MCMCreport= *pad_MCMCreport;
  j=1;
  NLL=0;
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
 //NLL+=0.5*log(2*3.141593*square(.039))+square(r-.042)/(2*square(.039)); // prior on r is normal (for now)
  // If in MCMC phase, print out the variable values (i.e. this is one iteration)
 if(mceval_phase()) MCMCreport << r << "," << K << endl;
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  cout << "NLL is " << NLL << endl;
  cout << "r is " << r << endl;
  cout << "K is " << K << endl;
  cout << "fpen is " << fpen << endl;
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
