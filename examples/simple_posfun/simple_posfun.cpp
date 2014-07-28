#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <simple_posfun.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nobs.allocate("nobs");
  Y.allocate(1,nobs,"Y");
  x.allocate(1,nobs,"x");
 ad_comm::change_datafile_name("bounds.txt");
  upb_b.allocate("upb_b");
 cout << upb_b << endl;
}

void model_parameters::initializationfunction(void)
{
  a.set_initial_value(0);
  b.set_initial_value(0);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  a.allocate("a");
  b.allocate(-10,upb_b,"b");
  aa.allocate("aa");
  pred_Y.allocate(1,nobs,"pred_Y");
  #ifndef NO_AD_INITIALIZE
    pred_Y.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  temp.allocate("temp");
  #ifndef NO_AD_INITIALIZE
  temp.initialize();
  #endif
  f.allocate("f");
}

void model_parameters::userfunction(void)
{
  fpen=0;
  b=-posfun(-b, -(upb_b-.0001), fpen);
  pred_Y=a*x+b;
 // cout << fpen<< " " <<  b << " " << upb_b-.01 << endl;
 aa=a;
  f=(norm2(pred_Y-Y));
  f=nobs/2.*log(f);    // make it a likelihood function so that
  f+=100000*fpen;                       // covariance matrix is correct
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

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
