
// Modified simple model to include an input for upper bound on b
DATA_SECTION
  init_int nobs
  init_vector Y(1,nobs)
  init_vector x(1,nobs)
  !! ad_comm::change_datafile_name("bounds.txt");
  init_number upb_b
  !! cout << upb_b << endl;
PARAMETER_SECTION
  init_number a
  init_bounded_number b(-10,upb_b)
  sdreport_number aa
  vector pred_Y(1,nobs)
  objective_function_value f
PROCEDURE_SECTION
  pred_Y=a*x+b;
 aa=a;
  f=(norm2(pred_Y-Y));
  f=nobs/2.*log(f);    // make it a likelihood function so that
                       // covariance matrix is correct
C:\Program Files\R\R-3.0.2\bin\x64;c:/ADMB/admb101-gcc452-win32/bin;c:/gnu/GCC452-WIN32/bin;c:/gnu/GDB/bin;C:\Program Files (x86)\GNU Emacs 24.3\bin

c:\Rtools\bin;c:\Rtools\gcc-4.6.3\bin;C:\WINDOWS\system32;C:\WINDOWS;C:\WINDOWS\System32\Wbem;C:\WINDOWS\System32\WindowsPowerShell\v1.0\;C:\Program Files (x86)\UWICK\SSH Tectia\SSH Tectia AUX;C:\Program Files (x86)\UWICK\SSH Tectia\SSH Tectia AUX/Support binaries;C:\Program Files (x86)\UWICK\SSH Tectia\SSH Tectia Broker;C:\Program Files (x86)\UWICK\SSH Tectia\SSH Tectia Client;C:\Program Files\MiKTeX 2.9\miktex\bin\x64\