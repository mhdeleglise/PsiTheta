#ifndef psi_theta_h
#define psi_theta_h
#include<math.h>
#include<mpfr.h>

class Psi_computer{

public:
  Psi_computer() {
    _sqrt_xm=NULL;
    _xm=NULL;
    _sumlog_table=NULL;
    _psi_tab=NULL;
    _mu=NULL;
    _pi=NULL;
    _p=NULL;
  };
  
  void psi(mpfr_t psi_value, long x, long decimal_relative_error, int v=0);
  void set_verbose(int v) {verbose=v;}
  ~Psi_computer();
private:
  // number of bits to use while computing to ensure a correct answer with decimal_prec decimal_places
  // 0.66*log(x) is an approximation of the total number of artithmetic operations done
  int nnbits(long x, int decimal_relative_error) {return ceil(3.5*decimal_relative_error + 0.65*log(x) + 10);} 
  void set_binary_precision(int nbits) {
    mpfr_set_default_prec(nbits);
    init_mpfr_vars(nbits);
  }

  // This function returns the range of the last bernoulli number needed to
  // get the precison requested 
  //long imax_bernoulli(long x, long decimal_relative_error);
  
  
  static const int MAX_SMALL_X=300;
  static int const SUMLOG_MAX_INDEX=50000;
  static int const MSTART_MACLAURIN=SUMLOG_MAX_INDEX+1;
  static int const BERNOULLI_SIZE=21;
  static int _small_primes[63];


  int imaxb;
  int verbose;

  int *_pi;
  int *_p;
  long* _mu;
  mpfr_t *_psi_tab;
  mpfr_t *_sumlog_table;
  long *_xm;
  long *_sqrt_xm;

  
  long pimax;
  long pmax;
  long maxp;
  long a;
  long b;
  
  mpfr_t _tmp1;
  mpfr_t _big_sumlog;
  mpfr_t _small_sumlog;
  mpfr_t _sumlog;
  mpfr_t _mcarre;
  mpfr_t _ncarre;
  mpfr_t _n;
  mpfr_t _m;
  mpfr_t _n1;
  mpfr_t _m1;
  mpfr_t _logn;
  mpfr_t _logm;
  mpfr_t _den;
  mpfr_t _ti;
  
  mpfr_t _bernoulli[BERNOULLI_SIZE];
  
  mpfr_t _last_psi;
  mpfr_t _psi_u;
  mpfr_t _S_1;
  mpfr_t _S_2;
  mpfr_t _S_3;
  mpfr_t _S_4;
  mpfr_t _lambda;
  mpfr_t _zzz;
  mpfr_t _longres;

  void init_mpfr_vars(int nbits);
  void clear_mpfr_vars();
  void init_p(long x);
  void init_pi(long x);
  void init_mu(long x);
  void build_mu(long a,long b);
  void init_psi_tab(long l);
  void clear_psi_tab(long l);
  void build_psi(long a,long b);
  void init_bernoulli();
  void init_bernoulli(long x);
  void clear_bernoulli();
  void init_sumlog();
  void clear_sumlog();
  void euler_mac_laurin_sum_log(long m, long n);
  void small_sumlog(long n);
  void sumlog(long n);
  void init_xm(long x,long u);
  void init(long x,long l,long u);

  inline long     mu(long n)     {return _mu[long(n)];}

  void compute_lambda(long n);
  void compute_small_psi(mpfr_t psi_value, long x, long nbits);
};


template <class T_INT> inline T_INT MAX(T_INT a,T_INT b) {return (a>b)? a:b;}
template <class T_INT> inline T_INT MIN(T_INT a,T_INT b) {return (a<b)? a:b;}
template <class T_INT> inline T_INT ABS(T_INT a) {return (a>=0)? a:(-a);}
template <class T_INT> inline T_INT ODD(T_INT a) {return long(a)&1L;}


template <class T_INT> T_INT power(T_INT x,int n)
{
  T_INT r;
  if (ODD(n)) r=x;
  else        r=1;
  while ((n>>=1))
    {
      x *= x;
      if (ODD(n)) r*=x;
    }
  return r;
}

template <class T_INT> inline T_INT root(T_INT x,int n)
{
  T_INT u,v;
  u = 1;
  v = x/n+1;
  while (ABS(v-u)>1)
    {
      u=v;
      v=x;
      for (int i=0;i<n-1;i++) v/=u;
      v = (u*(n-1)+v)/n;
    }
  return MIN(u,v);
}

inline long root(long x,int n)
{
  long u,v;
  u = 1;
  v = (x-1)/n+1;
  while (ABS(v-u)>1)
    {
      u=v;
      v=x;
      for (int i=0;i<n-1;i++) v/=u;
      v = (u*(n-1)+v)/n;
    }
  return MIN(u,v);
}


class Theta_computer{

public:
  void psi_theta_diff(mpfr_t res, long x, long decimal_relative_error);
  void theta(mpfr_t res, long x, long decimal_relative_error);

private:
  mpfr_t _logp;
  mpfr_t _pfloat;
};
  
#endif
