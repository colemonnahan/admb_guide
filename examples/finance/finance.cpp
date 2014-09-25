#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <finance.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  T.allocate("T");
  r.allocate(0,T,"r");
  sub_r.allocate(1,T);
 ad_comm::change_datafile_name("phases.dat");
  phasea0.allocate("phasea0");
  phasea1.allocate("phasea1");
  phasea2.allocate("phasea2");
  phaseMean.allocate("phaseMean");
}

void model_parameters::initializationfunction(void)
{
  a0.set_initial_value(.00016);
  a1.set_initial_value(.1);
  a2.set_initial_value(.37);
  Mean.set_initial_value(0);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  a0.allocate(0.0,.01,phasea0,"a0");
  a1.allocate(0.0,.5,phasea1,"a1");
  a2.allocate(0.0,1.0,phasea2,"a2");
  Mean.allocate(-.001,.001,phaseMean,"Mean");
  eps2.allocate(1,T,"eps2");
  #ifndef NO_AD_INITIALIZE
    eps2.initialize();
  #endif
  h.allocate(1,T,"h");
  #ifndef NO_AD_INITIALIZE
    h.initialize();
  #endif
  aa.allocate("aa");
  log_likelihood.allocate("log_likelihood");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  h0=square(std_dev(r));   // square forms the element-wise square
  sub_r=r(1,T);    // form a subvector so we can use vector operations
  Mean=mean(r);    // calculate the mean of the vector r
}

void model_parameters::userfunction(void)
{
  log_likelihood =0.0;
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::report(const dvector& gradients){}

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
    if (!arrmblsize) arrmblsize=15000000;
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
