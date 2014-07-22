// TOok the original cpp file and chopped out bits I don't care about so it
// is more readable  
void function_minimizer::computations1(int argc,char * argv[])
{
  tracing_message(traceflag,"B1");
 
  int on=-1;
  int nopt=-1;
 
  set_runtime();
 
  if ( (on=option_match(argc,argv,"-hbf",nopt))>-1)
    {
      gradient_structure::Hybrid_bounded_flag=1;
    }
 
  // Sets the maximum number of function evaluation as determined from the
  // command line
  if ( (on=option_match(argc,argv,"-maxfn",nopt))>-1)
    {
      if (nopt ==1)
	{
	  set_runtime_maxfn(argv[on+1]);
	}
      else
	{
	  cerr << "Wrong number of options to -mafxn -- must be 1"
	    " you have " << nopt << endl;
	}
    }
 

  stddev_params::get_stddev_number_offset();
 
  tracing_message(traceflag,"C1");
 
 
  repeatminflag=0;
  do
    {
      if (option_match(argc,argv,"-noest") == -1)
	{
	  if (!function_minimizer::have_constraints)
	    {
	      minimize();
	    }
	  else
	    {
	      constraints_minimize();
	    }
	}
      else
	{
	  initial_params::current_phase=initial_params::max_number_phases;
	}
      tracing_message(traceflag,"D1");
 
      //double ratio=100.*gradient_structure::max_last_offset/12000.0;
      tracing_message(traceflag,"E1");

      if (option_match(argc,argv,"-est") == -1)
	{
	  if (!quit_flag)
	    {
	      int on=-1;
	      int on1=-1;
	      on=option_match(argc,argv,"-nohess");
	      on1=option_match(argc,argv,"-noest");
	      if (on==-1 && on1==-1)
		{
		    {
		      depvars_routine();
		      hess_inv();
		      if (spminflag==0)
			{
			  sd_routine();
			}
		    }
		}
	      else
		{
		  initial_params::sd_phase=1;
		}
	      if (spminflag==0)
		{
		  if ( (on=option_match(argc,argv,"-lprof"))>-1)
		    {
		      if (likeprof_params::num_likeprof_params)
			{
			    {
			      likeprof_routine(ffbest);
			    }
			}
		    }
		  int nopt=0;
		  int on2=-1;
		  int nopt2=-1;
 
		  // stuff for mcmc
		  //cout << "checking for mcmc" << endl;
		  if ( (on=option_match(argc,argv,"-mcmc",nopt))>-1 ||
		       (on=option_match(argc,argv,"-mcmc2",nopt))>-1)
		    {
		      if ( (on2=option_match(argc,argv,"-mcmc2",nopt2))>-1)
			mcmc2_flag=1;
		      else
			mcmc2_flag=0;
 
			{
			  mcmc_computations();
			}
		    }
		  initial_params::sd_phase=0;
		}
	    }
	}
      // end of if -est not specified
    } // end of do loop
  while(spminflag || repeatminflag);
}
