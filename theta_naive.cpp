#include<cmath>
#include"sieve_slices.h"
#include<mpfr.h>


int main(int argc, char* argv[]) {

  long N= atol(argv[1]);
  long nbits= (argc > 2) ? atol(argv[2]) : 128;
  mpfr_set_default_prec(nbits);

  cout << "N = " << N << endl;

  mpfr_t tmp;
  mpfr_t res;
  mpfr_init_set_si(res,0,MPFR_RNDN);
  mpfr_init_set_si(tmp,0,MPFR_RNDN);
  
  long wsize=1000000;
  cout << "Creation d'un prime generator de taille " << wsize << endl;
  prime_generator pg(100000,wsize);
  
  long cnte=0;
  
  long p = pg.next_prime();
  mpfr_set_si(res,1,MPFR_RNDN);
  mpfr_log(res,res,MPFR_RNDN);
  while (p <= N) {
    cnte++;
    mpfr_set_si(tmp,p,MPFR_RNDN);
    mpfr_log(tmp,tmp,MPFR_RNDN);
    mpfr_add(res,res,tmp,MPFR_RNDN);
    p=pg.next_prime();
  }
  mpfr_printf("theta= %.35Re\n",res);
  return 0;

}
