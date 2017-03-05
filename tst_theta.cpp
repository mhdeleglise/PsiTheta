#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <mpfr.h>
#include "psitheta.h"


using namespace std;


int main(int argc,char *argv[])
{
  long x;
  long decimal_prec= (argc > 2)? atol(argv[2]) : 25;
  istringstream x_stream(argv[1],istringstream::in);
  x_stream >> x;
  
  Theta_computer theta_comp;
  mpfr_t theta_res;
  mpfr_init_set_si(theta_res,0,MPFR_RNDN);
  theta_comp.theta(theta_res, x, decimal_prec);

  mpfr_out_str(stdout,10,decimal_prec,theta_res, MPFR_RNDN);printf("\n");

  return 0;

}
