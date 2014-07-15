#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <forest.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nsteps.allocate("nsteps");
  k.allocate("k");
  a.allocate(1,k+1,"a");
  freq.allocate(1,k,"freq");
 sum_freq=sum(freq);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_tau.allocate(-14,15,2,"log_tau");
  log_nu.allocate(-15,4,"log_nu");
  beta.allocate(.1,1.0,-1,"beta");
  log_sigma.allocate(-5,3,"log_sigma");
  tau.allocate("tau");
 tau.set_stepnumber(25);
  nu.allocate("nu");
  sigma.allocate("sigma");
  S.allocate(1,k+1,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  f.allocate("f");
}

void model_parameters::initializationfunction(void)
{
  log_tau.set_initial_value(0);
  beta.set_initial_value(0.6666667);
  log_nu.set_initial_value(0);
  log_sigma.set_initial_value(-2);
}

void model_parameters::userfunction(void)
{
  tau=exp(log_tau);
  nu=exp(log_nu);
  sigma=exp(log_sigma);
   funnel_dvariable Integral;
   int i;
   for (i=1;i<=k+1;i++)
   {
     a_index=i;
     ad_begin_funnel();
     Integral=adromb(&model_parameters::h,-3.0,3.0,nsteps);
     S(i)=Integral;
   }
   f=0.0;
   for (i=1;i<=k;i++)
   {
     dvariable ff=0.0;
     dvariable diff=posfun((S(i)-S(i+1))/S(i),.000001,ff);
     f-=freq(i)*log(1.e-50+S(i)*diff);
     f+=ff;
   }
   f+=sum_freq*log(1.e-50+S(1));
}

dvariable model_parameters::h(const dvariable& z)
{
  dvariable tmp;
  tmp=exp(-.5*z*z + tau*(-1.+exp(-nu*pow(a(a_index),beta)*exp(sigma*z))) );  
  return tmp;
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
  report << "nsteps = " << std::setprecision(10) <<  nsteps << endl;
  report << "f = " << std::setprecision(10) <<  f << endl;
  report << "a" << endl << a << endl;
  report << "freq" << endl << freq << endl;
  report << "S" << endl << S << endl;
  report << "S/S(1)" << endl << std::ios::fixed << std::setprecision(6) << S/S(1) << endl;
  report << "tau "  << tau << endl; 
  report << "nu "  << nu << endl; 
  report << "beta "  << beta << endl; 
  report << "sigma "  << sigma << endl; 
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

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
