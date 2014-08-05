#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <simple.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_leapfrog = new ofstream("leapfrog.csv",ios::trunc);;
  nobs.allocate("nobs");
  Y.allocate(1,nobs,"Y");
  x.allocate(1,nobs,"x");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  a.allocate("a");
  b.allocate("b");
  aa.allocate("aa");
  pred_Y.allocate(1,nobs,"pred_Y");
  #ifndef NO_AD_INITIALIZE
    pred_Y.initialize();
  #endif
  f.allocate("f");
}

void model_parameters::userfunction(void)
{
  ofstream& leapfrog= *pad_leapfrog;
  pred_Y=a*x+b;
 aa=a;
  f=(norm2(pred_Y-Y));
  f=nobs/2.*log(f);    // make it a likelihood function so that
                       // covariance matrix is correct
  leapfrog << a << "," << b << "," << f << endl;
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
 cout << a << ", " << b << endl;
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_leapfrog;
  pad_leapfrog = NULL;
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
