#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <truncreg.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nobs.allocate("nobs");
  m.allocate("m");
  trunc_flag.allocate("trunc_flag");
  data.allocate(1,nobs,1,m+1,"data");
  Y.allocate(1,nobs);
  X.allocate(1,nobs,1,m);
  Y=column(data,1);
  for (int i=1;i<=nobs;i++)
  {
    X(i)=data(i)(2,m+1).shift(1);
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  sigma.allocate("sigma");
  vhat.allocate("vhat");
  #ifndef NO_AD_INITIALIZE
  vhat.initialize();
  #endif
  log_a.allocate(-5.0,5.0,"log_a");
  a.allocate("a");
  u.allocate(1,m,"u");
  f.allocate("f");
}

void model_parameters::userfunction(void)
{
  a=exp(log_a);
  dvar_vector pred=X*u;
  dvar_vector res=Y-pred;
  dvariable r2=norm2(res); 
  vhat=r2/nobs; 
  dvariable v=a*vhat;
  sigma=sqrt(v);
  dvar_vector spred=pred/sigma;
  f=0.0;
  switch (trunc_flag)
  {
  case -1:  // left_truncated
    {
      for (int i=1;i<=nobs;i++)
      {
        f+=log(1.00001-cumd_norm(-spred(i)));
      }
    }
    break;
  case 1:   // right truncated
    {
      for (int i=1;i<=nobs;i++)
      {
        f+=log(0.99999*cumd_norm(-spred(i)));
      }
    }
    break;
  case 0:   // no truncation
    break;
  default:
    cerr << "Illegal value for truncation flag" << endl;
    ad_exit(1);
  }
  f+=0.5*nobs*log(v)+0.5*r2/v;
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
  report << "#u " << endl << u << endl;
  report << "#sigma " << endl << sigma << endl;
  report << "#a " << endl << a << endl;
  report << "#vhat " << endl << vhat << endl;
  report << "#shat " << endl << sqrt(vhat) << endl;
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
