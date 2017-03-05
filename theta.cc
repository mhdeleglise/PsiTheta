#include<iostream>
#include<sstream>
#include<iomanip>
#include"psitheta.h"
#include"sieve_slices.h"

void Theta_computer::psi_theta_diff(mpfr_t res, long x, long decimal_relative_error) {
  mpfr_t sum;
  mpfr_init2(res,ceil(3.3*(1+decimal_relative_error)));
  mpfr_set_default_prec(ceil(3.3*decimal_relative_error + 0.5*log2(x)));
  mpfr_init_set_si(_logp, 0, MPFR_RNDN);
  mpfr_init_set_si(_pfloat, 0, MPFR_RNDN);
  mpfr_init_set_si(sum, 0, MPFR_RNDN);
  prime_generator pg(100000,1000000,1);
  long p=pg.next_prime();
  long xp=x/p;
  while (p <= xp) {
    mpfr_set_si(_pfloat, p, MPFR_RNDN);
    mpfr_log(_logp, _pfloat, MPFR_RNDN);
    while (p <= xp) {
      while (p <= xp) {
	mpfr_add(sum, sum, _logp, MPFR_RNDN);
	xp /= p;
      }
    }
    p = pg.next_prime();
    xp = x/p;
  }
  mpfr_set(res, sum, MPFR_RNDN);
  mpfr_clear(_logp);
  mpfr_clear(_pfloat);
  mpfr_clear(sum);
}


void Theta_computer::theta(mpfr_t res, long x, long decimal_relative_error) {
  mpfr_init2(res,ceil(3.3*(1+decimal_relative_error)));
  mpfr_t res1,res2;
  mpfr_init(res1);
  mpfr_init(res2);
  mpfr_init2(res,ceil(3.3*(1+decimal_relative_error)));
  mpfr_init2(res,ceil(3.3*(1+decimal_relative_error)));

  Psi_computer psi_comp;
  psi_comp.psi(res1, x, decimal_relative_error);
  psi_theta_diff(res2, x, decimal_relative_error);
  mpfr_set_prec(res,ceil(3.3*(1+decimal_relative_error)));
  mpfr_sub(res,res1,res2,MPFR_RNDN);
}

