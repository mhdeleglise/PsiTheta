#include"psitheta.h"
#include<iostream>
#include<cstdlib>

using namespace std;


int Psi_computer::_small_primes[63]= {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,\
113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307};

void Psi_computer::init_mpfr_vars(int nnbits) {
  mpfr_set_default_prec(nnbits);
  mpfr_init_set_si(_tmp1,0, MPFR_RNDN);
  mpfr_init_set_si(_big_sumlog,0, MPFR_RNDN);
  mpfr_init_set_si(_small_sumlog,0, MPFR_RNDN);
  mpfr_init_set_si(_sumlog,0, MPFR_RNDN);
  mpfr_init_set_si(_mcarre,0, MPFR_RNDN);
  mpfr_init_set_si(_ncarre,0, MPFR_RNDN);
  mpfr_init_set_si(_n,0, MPFR_RNDN);
  mpfr_init_set_si(_m,0, MPFR_RNDN);
  mpfr_init_set_si(_n1,0, MPFR_RNDN);
  mpfr_init_set_si(_m1,0, MPFR_RNDN);
  mpfr_init_set_si(_logn,0, MPFR_RNDN);
  mpfr_init_set_si(_logm,0, MPFR_RNDN);
  mpfr_init_set_si(_den,0, MPFR_RNDN);
  mpfr_init_set_si(_ti,0, MPFR_RNDN);
  mpfr_init_set_si(_last_psi,0, MPFR_RNDN);
  mpfr_init_set_si(_psi_u,0, MPFR_RNDN);
  mpfr_init_set_si(_S_1,0, MPFR_RNDN);
  mpfr_init_set_si(_S_2,0, MPFR_RNDN);
  mpfr_init_set_si(_S_3,0, MPFR_RNDN);
  mpfr_init_set_si(_S_4,0, MPFR_RNDN);
  mpfr_init_set_si(_lambda,0, MPFR_RNDN);
  mpfr_init_set_si(_zzz,0, MPFR_RNDN);
  mpfr_init_set_si(_longres,0, MPFR_RNDN);
  init_bernoulli();
}

void Psi_computer::clear_mpfr_vars() {
  mpfr_clears(_tmp1, _big_sumlog, _small_sumlog, _sumlog, _mcarre, _ncarre, _n, _m, _n1, _m1, _logn, _logm, _den, _ti, _last_psi,\
	      _psi_u, _S_1, _S_2, _S_3, _S_4, _lambda, _zzz, _longres, (mpfr_ptr) 0);
}

void Psi_computer::init_pi(long x)
{
  pimax=0;
  if (x<=0) {cerr << "bad x in init_pi, x = " << x << endl;exit(1);}
  _pi = new int[x+1];
  if (_pi==NULL) {cerr << "no memory for _pi !" << endl;exit(1);}
  long i;
  _pi[0] = 0; _pi[1] = 0; 
  for (i=2;i<=x;i++) _pi[i] = 1;
  long maxp = root(x,2);
  for (long p=2;p<=maxp;p++)
    {
      if (_pi[p])
	for (long m=p*p;m<=x;m+=p) _pi[m] = 0;
    }
  for (i=1;i<=x;i++) _pi[i] += _pi[i-1];
  pimax = _pi[x];
}

void Psi_computer::init_p(long x)
{
  pmax=0;
  maxp=0;
  int i,pi;
  _p = new int[pimax+1];
  if (_p==NULL) {cerr << "no memory for _p !" << endl; exit(1);}
  _p[0] = 0;
  pi = 0;
  for (i=2;i<=x;i++) 
    if (_pi[i]-_pi[i-1]) _p[++pi] = i;
  pmax = _p[pimax];
  maxp = x;
}

void Psi_computer::init_mu(long x)
{
  _mu = new long[x+1];
  if (_mu==NULL) {cerr << "no memory for _mu !" << endl; exit(1);}
}


void Psi_computer::build_mu(long a,long b)
{
  long i;
  long l = b-a;
  long sqrt_b = root(b,2);
  if (maxp<sqrt_b)
    {
      cerr << "Insufficient size of primetable" << endl;
      cerr << "maxp = " << maxp << ", sqrt(b) = " << sqrt_b << endl;
      exit(1);
    }
  for (i=0;i<l;i++) _mu[i] = 1;
  for (i=1;(i<=pimax)&&(_p[i]<=sqrt_b);i++)
    {
      long p = _p[i];
      long p2 = p; p2 *= p;
      // eliminate square factors
      long n  = -(a%p2);
      for (n=(n<0)?n+p2:n;n<l;n+=p2) _mu[long(n)] = 0;
      // compute _mu[]
      n = -(a%p);
      for (n=(n<0)?n+p:n;n<l;n+=p) _mu[long(n)] *= -p;
    }
  for (i=0;i<l;i++)
    {
      if (_mu[i]!=0)
	{
	  if (ABS(_mu[i])<(a+i)) _mu[i] = -_mu[i];
	  if (_mu[i]>0) _mu[i] =  1;
	  else        _mu[i] = -1;
	}
    }
}

void Psi_computer::init_psi_tab(long l)
{
  _psi_tab = new mpfr_t[l+1];
  if (_psi_tab==NULL) {cerr << "no memory for _psi_tab !" << endl; exit(1);}
  for (int i=0; i <= l; i++) {
    mpfr_init_set_si(_psi_tab[i],i, MPFR_RNDN);
  }
}

void Psi_computer::clear_psi_tab(long l) {
  for (int i=0; i <= l; i++)
    mpfr_clear(_psi_tab[i]);
}

void Psi_computer::build_psi(long a,long b)
{
  long n;
  long l = b-a;
  long sqrt_b = root(b,2);
  if (maxp<sqrt_b)
    {
      cerr << "Insufficient size of primetable" << endl;
      cerr << "maxp = " << maxp << ", sqrt(b) = " << sqrt_b << endl;
      exit(1);
    }
  for (n=0;n<l;n++) mpfr_set_si(_psi_tab[n],1, MPFR_RNDN);
  if (a==0) {
    mpfr_set_si(_psi_tab[0],0,MPFR_RNDN);
    mpfr_set_si(_psi_tab[1],0,MPFR_RNDN);
      }
  if (a==1)
    mpfr_set_si(_psi_tab[0],0,MPFR_RNDN);

  for (long k=1;(k<=pimax)&&(_p[k]<=sqrt_b);k++)
    {
      long p = _p[k];
      long m = MAX(p*p,a-(a%p));
      for (m = (m<a)?m+p:m; m<b; m+=p)	mpfr_set_si(_psi_tab[long(m-a)],0,MPFR_RNDN);
      long i=2;
      for (m=p*p;m<a;m*=p) i++;
      while (m<b)
	{
	  mpfr_set_si(_psi_tab[long(m-a)],i, MPFR_RNDN);
	  m*=p;i++;
	}
    }
  for (n=0;n<l;n++)
    {
      if (mpfr_zero_p(_psi_tab[n])) {
	mpfr_set(_psi_tab[n],_last_psi,MPFR_RNDN);
      }
      else  {
	// _psi_tab[n] = LOG(a+n)/_psi_tab[n]+lastpsi;
	mpfr_set_si(_tmp1, a+n, MPFR_RNDN);
	mpfr_log(_tmp1, _tmp1, MPFR_RNDN);
	mpfr_div(_tmp1,_tmp1,_psi_tab[n], MPFR_RNDN);
	mpfr_add(_psi_tab[n], _tmp1, _last_psi, MPFR_RNDN);
      }
      mpfr_set(_last_psi, _psi_tab[n], MPFR_RNDN);
    }
}

void Psi_computer::init_bernoulli(long x) {
  //cout << "In init_bernoulli(x) x= " << x << endl;
  long m=MSTART_MACLAURIN;
  long n=x/root(x,3);
  mpfr_set_si(_m,m, MPFR_RNDN);
  mpfr_set_si(_n,n,MPFR_RNDN);
  mpfr_set_si(_ncarre, n, MPFR_RNDN); 
  mpfr_set_si(_mcarre, m, MPFR_RNDN);
  mpfr_mul_si(_ncarre, _ncarre, n, MPFR_RNDN);
  mpfr_mul_si(_mcarre, _mcarre, m, MPFR_RNDN);
  mpfr_set_si(_tmp1,1,MPFR_RNDN);
  mpfr_div(_mcarre,_tmp1, _mcarre, MPFR_RNDN);
  mpfr_div(_ncarre,_tmp1, _ncarre, MPFR_RNDN);
  
  long signe = -1 ;
  mpfr_set(_n1, _n, MPFR_RNDN);
  mpfr_set(_m1, _m, MPFR_RNDN);

  long negexp;
  for (int i=1;i<BERNOULLI_SIZE;i++)
    { 
      mpfr_mul(_n1, _n1, _ncarre, MPFR_RNDN);
      mpfr_mul(_m1, _m1, _mcarre, MPFR_RNDN);
      signe = -signe;
      mpfr_set_d(_den, signe * (2.0*i)*(i+i-1), MPFR_RNDN);
      mpfr_div(_ti,_bernoulli[i],_den, MPFR_RNDN);
      mpfr_sub(_tmp1, _n1, _m1, MPFR_RNDN);
      mpfr_mul(_ti, _ti, _tmp1, MPFR_RNDN);
      mpfr_get_d_2exp(&negexp,_ti, MPFR_RNDN);
      if (abs(negexp) > mpfr_get_default_prec()) {
	//cout << "Here negexp= " << negexp << "  >  " << mpfr_get_default_prec() << "   we stop with imax = " << i << endl;
	imaxb=i;
	return;
      }
    }
  cout << "It is not absolutly sure than we get the request precision\n";
  exit(0);
}



void Psi_computer::init_bernoulli() {
  mpfr_init_set_si(_bernoulli[0],0,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[1],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[1], _bernoulli[1], 6, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[2],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[2], _bernoulli[2], 30, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[3],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[3], _bernoulli[3], 42, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[4],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[4], _bernoulli[4], 30, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[5],5,MPFR_RNDN);
  mpfr_div_si(_bernoulli[5], _bernoulli[5], 66, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[6],691,MPFR_RNDN);
  mpfr_div_si(_bernoulli[6], _bernoulli[6], 2730, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[7],7,MPFR_RNDN);
  mpfr_div_si(_bernoulli[7], _bernoulli[7], 6, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[8],3617,MPFR_RNDN);
  mpfr_div_si(_bernoulli[8],_bernoulli[8], 510,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[9],43867,MPFR_RNDN);
  mpfr_div_si(_bernoulli[9],_bernoulli[9], 798,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[10],174611,MPFR_RNDN);
  mpfr_div_si(_bernoulli[10],_bernoulli[10], 330,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[11],854513,MPFR_RNDN);
  mpfr_div_si(_bernoulli[11],_bernoulli[11], 138,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[12],236364091,MPFR_RNDN);
  mpfr_div_si(_bernoulli[12],_bernoulli[12], 2730,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[13],8553103,MPFR_RNDN);
  mpfr_div_si(_bernoulli[13],_bernoulli[13], 6,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[14],23749461029,MPFR_RNDN);
  mpfr_div_si(_bernoulli[14],_bernoulli[14], 870,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[15],8615841276005,MPFR_RNDN);
  mpfr_div_si(_bernoulli[15],_bernoulli[15], 14322,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[16],7709321041217,MPFR_RNDN);
  mpfr_div_si(_bernoulli[16],_bernoulli[16], 510,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[17],2577687858367,MPFR_RNDN);
  mpfr_div_si(_bernoulli[17],_bernoulli[17], 6,MPFR_RNDN);
  mpfr_init_set_str(_bernoulli[18],"26315271553053477373",10,MPFR_RNDN);
  mpfr_div_si(_bernoulli[18],_bernoulli[18], 1919190,MPFR_RNDN);
  mpfr_init_set_str(_bernoulli[19],"2929993913841559",10,MPFR_RNDN);
  mpfr_div_si(_bernoulli[19],_bernoulli[19], 6,MPFR_RNDN);
  mpfr_init_set_str(_bernoulli[20],"261082718496449122051",10,MPFR_RNDN);
  mpfr_div_si(_bernoulli[20],_bernoulli[20], 13530,MPFR_RNDN);
}


void Psi_computer::clear_bernoulli() {
  for (int i=0; i < BERNOULLI_SIZE; i++)
    mpfr_clear(_bernoulli[i]);
}

void Psi_computer::init_sumlog()
{
  _sumlog_table = new mpfr_t[SUMLOG_MAX_INDEX+1];
  mpfr_init_set_si(_sumlog_table[0],0,MPFR_RNDN);
  for (long i=1;i<=SUMLOG_MAX_INDEX;i++) {
    mpfr_init_set_si(_sumlog_table[i], i, MPFR_RNDN);
    mpfr_log(_sumlog_table[i], _sumlog_table[i], MPFR_RNDN);
    mpfr_add(_sumlog_table[i], _sumlog_table[i], _sumlog_table[i-1], MPFR_RNDN);
  }
}

void Psi_computer::clear_sumlog() {
  for (int i=0; i <= SUMLOG_MAX_INDEX; i++)
    mpfr_clear(_sumlog_table[i]);
}


void Psi_computer::euler_mac_laurin_sum_log(long m, long n)
{
  //printf("In euler_mac m =   %ld      n= %ld   \n",m,n);
  mpfr_set_si(_big_sumlog,0, MPFR_RNDN);
  mpfr_set_si(_n, n, MPFR_RNDN);
  mpfr_set_si(_m, m, MPFR_RNDN);
  mpfr_set_si(_ncarre, n, MPFR_RNDN); 
  mpfr_set_si(_mcarre, m, MPFR_RNDN);
  mpfr_mul_si(_ncarre, _ncarre, n, MPFR_RNDN);
  mpfr_mul_si(_mcarre, _mcarre, m, MPFR_RNDN);
  mpfr_set_si(_tmp1,1,MPFR_RNDN);
  mpfr_div(_mcarre,_tmp1, _mcarre, MPFR_RNDN);
  mpfr_div(_ncarre,_tmp1, _ncarre, MPFR_RNDN);

  mpfr_log(_logn, _n, MPFR_RNDN);
  mpfr_mul_d(_logn, _logn, n+0.5, MPFR_RNDN); 
  mpfr_log(_logm, _m, MPFR_RNDN);
  mpfr_mul_d(_logm, _logm, m-0.5, MPFR_RNDN);
  mpfr_add(_big_sumlog, _big_sumlog, _logn, MPFR_RNDN);
  mpfr_sub(_big_sumlog, _big_sumlog, _logm, MPFR_RNDN);
  mpfr_sub_si(_big_sumlog,_big_sumlog, n-m, MPFR_RNDN);

  long signe = -1 ;
  mpfr_set(_n1, _n, MPFR_RNDN);
  mpfr_set(_m1, _m, MPFR_RNDN);

  for (int i=1;i<imaxb;i++)
    { 
      mpfr_mul(_n1, _n1, _ncarre, MPFR_RNDN);
      mpfr_mul(_m1, _m1, _mcarre, MPFR_RNDN);
      signe = -signe;
      mpfr_set_d(_den, signe * (2.0*i)*(i+i-1), MPFR_RNDN);
      mpfr_div(_ti,_bernoulli[i],_den, MPFR_RNDN);
      mpfr_sub(_tmp1, _n1, _m1, MPFR_RNDN);
      mpfr_mul(_ti, _ti, _tmp1, MPFR_RNDN);
      mpfr_add(_big_sumlog, _big_sumlog, _ti, MPFR_RNDN);
    }
  //mpfr_get_d_2exp(&negexp,_ti, MPFR_RNDN);
  //mpfr_printf("BERNOULLI_SIZE=%d  m=%ld   n=%ld  negexp= %ld  ti= %.30Re\n",BERNOULLI_SIZE,m,n,negexp,_ti);
}


void Psi_computer::small_sumlog(long n)  
{
  mpfr_set(_small_sumlog, _sumlog_table[n], MPFR_RNDN);
}


void Psi_computer::sumlog(long n)
{
  if (n<=SUMLOG_MAX_INDEX) {
    small_sumlog(n);
    mpfr_set(_sumlog, _small_sumlog, MPFR_RNDN);
    return;
  }
  euler_mac_laurin_sum_log(SUMLOG_MAX_INDEX+1,n);
  small_sumlog(SUMLOG_MAX_INDEX);
  mpfr_add(_sumlog,_big_sumlog,_small_sumlog, MPFR_RNDN);
}

void Psi_computer::init_xm(long x,long u)
{
  _xm = new long[u+1];
  if (_xm==NULL) {cerr << "no memory for _xm !" << endl; exit(1);}
  _sqrt_xm = new long[u+1];
  if (_sqrt_xm==NULL) {cerr << "no memory for _sqrt_xm !" << endl; exit(1);}
  for (long m=1;m<=u;m++) {_xm[m] = x/m; _sqrt_xm[m]=root(_xm[m],2);}
}


void Psi_computer::init(long x,long l,long u)
{
  a=0;
  b=0;
  init_pi(root(x/u,2));
  init_p(root(x/u,2));
  init_mu(u+1);
  build_mu(0,u+1);
  init_xm(x,u);
  init_psi_tab(l);
  //time_msg((char*)"init_psi");
  init_sumlog();
  //time_msg((char*)"init_sumlog");
}


void Psi_computer::compute_lambda(long n) {
  mpfr_set(_lambda, _psi_tab[n-a], MPFR_RNDN);
  mpfr_sub(_lambda, _lambda, _psi_tab[n-a-1], MPFR_RNDN);
}

void Psi_computer::compute_small_psi(mpfr_t psi_value, long x, long decimal_prec) {
  int nbits=ceil(3.3*decimal_prec);
  mpfr_inits2(nbits+8, _longres, _tmp1, (mpfr_ptr) 0);
  mpfr_set_si(_longres,0,MPFR_RNDN);
  int ip=0;
  long p=_small_primes[ip]; 
  while (p <= x) {
    mpfr_set_si(_tmp1, p, MPFR_RNDN);
    mpfr_log(_tmp1, _tmp1, MPFR_RNDN);
    long q=p;
    while (q <= x) {
      mpfr_add(_longres, _longres, _tmp1, MPFR_RNDN);
      q *= p;
    }
    p=_small_primes[++ip];
  }
  mpfr_set(psi_value, _longres, MPFR_RNDN);
}

void Psi_computer::psi(mpfr_t psi_value, long x, long decimal_prec, int v)
{
  verbose = v;
  mpfr_set_prec(psi_value, nnbits(x,decimal_prec));
  if (verbose)
    cout << "precision " << nnbits(x,decimal_prec) << endl;
  if (x <= MAX_SMALL_X) {
    compute_small_psi(psi_value,x,decimal_prec);
    return;
  }
  long  x3 = root(x,3);
  long   u = x3;
  long   l = u*8;
  //cout << "Nombre de bits necessaires pour les clacul intermediaires " << nnbits(x,decimal_prec) << endl;
  //mpfr_set_default_prec(nnbits(decimal_prec));
  //cout << "binary precision set to " << nnbits(x,decimal_prec) << endl;
  set_binary_precision(nnbits(x,decimal_prec));
  init_bernoulli(x);
  init(x,l,u);
  a = 1;
  b = a+l;
  mpfr_set_si(_last_psi,0,MPFR_RNDN);
  build_psi(a,b);
  mpfr_set(_last_psi, _psi_tab[long(b-1-a)], MPFR_RNDN);
  mpfr_set(_psi_u, _psi_tab[long(u-a)], MPFR_RNDN);
  mpfr_set(_S_1, _psi_u, MPFR_RNDN);
  for (long n=1;n<=u;n++)
    {
      if (mu(n)==1)
	{
	  long xn = x/n;
	  mpfr_set_si(_sumlog, 0, MPFR_RNDN);
	  sumlog(xn);
	  mpfr_add(_S_2, _S_2, _sumlog, MPFR_RNDN);
	  mpfr_set_si(_zzz, xn/u,MPFR_RNDN);
	  mpfr_sub_si(_zzz,_zzz, u/n,MPFR_RNDN);
	  mpfr_mul(_tmp1, _psi_u, _zzz, MPFR_RNDN);
	  mpfr_sub(_S_4, _S_4, _tmp1, MPFR_RNDN);
	}
      if (mu(n)==-1)
	{
	  long xn = x/n;
	  sumlog(xn);
	  mpfr_sub(_S_2, _S_2, _sumlog, MPFR_RNDN);
	  mpfr_set_si(_zzz, xn/u,MPFR_RNDN);
	  mpfr_sub_si(_zzz,_zzz, u/n,MPFR_RNDN);
	  mpfr_mul(_tmp1, _psi_u, _zzz, MPFR_RNDN);
	  mpfr_add(_S_4, _S_4,_tmp1, MPFR_RNDN);
	}
    }
  for (long m=2;m<=u;m++)
    {
      compute_lambda(m);
      if (! mpfr_zero_p(_lambda))
	{
	  long xm = _xm[m];
	  long S  = 0;
	  for (long n=1;n<=u;n++)
	    {
	      if (mu(n)==1)  S += xm/n;
	      if (mu(n)==-1) S -= xm/n;
	    }
	  mpfr_set_si(_tmp1, S, MPFR_RNDN);
	  mpfr_mul_si(_lambda, _lambda, S, MPFR_RNDN);
	  mpfr_add(_S_3, _S_3, _lambda, MPFR_RNDN);
	}
    }
  long xu = x/u;
  while (a<xu)
    {
      for (long m = 1; m <= u; m++)
	{
	  if (mu(m)==0) continue;
	  long xm = _xm[m];
	  long n1 = MAX(long(u/m),xm/b)+1;
	  long n = n1;
	  long n2;
	  n2 = MIN(xm/a,long(xm/u));
	  if (mu(m)==1) 
	    {
	      while (n <= n2)
		{
		  long q = xm/n;
		  long nextn = xm/q + 1;
		  if (nextn <= n2) {
		    mpfr_set_si(_tmp1, nextn-n, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_add(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  else {
		    mpfr_set_si(_tmp1, n2-n+1, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_add(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  n = nextn;
		}
	    }
	  else
	    {
	      while (n <= n2)
		{
		  long q = xm/n;
		  long nextn = xm/q + 1;
		  if (nextn <= n2) {
		    // _S_4 -= psi(q)*LONG_DOUBLE(nextn-n);
		    mpfr_set_si(_tmp1, nextn-n, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_sub(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  else {
		    mpfr_set_si(_tmp1, n2-n+1, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_sub(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  n = nextn;
		}
	    }
	}
      a+=l;
      if (a>=xu) break;
      b = MIN(long(a+l),xu);
      build_psi(a,b);
      mpfr_set(_last_psi, _psi_tab[b-1-a], MPFR_RNDN);
    }
  mpfr_add(_longres, _longres, _S_1, MPFR_RNDN);
  mpfr_add(_longres, _longres, _S_2, MPFR_RNDN);
  mpfr_sub(_longres, _longres, _S_3, MPFR_RNDN);
  mpfr_sub(_longres, _longres, _S_4, MPFR_RNDN);
  mpfr_set(psi_value,_longres, MPFR_RNDN);
  clear_mpfr_vars();
  clear_psi_tab(l);
  clear_sumlog();
  clear_bernoulli();
}


Psi_computer::~Psi_computer() {
  if (_sqrt_xm)
    delete [] _sqrt_xm;
  if (_xm)
    delete [] _xm;
  if (_sumlog_table)
    delete [] _sumlog_table;
  if (_psi_tab)
    delete [] _psi_tab;
  if (_mu)
    delete [] _mu;
  if (_p)
    delete [] _p;
  if (_pi)
    delete [] _pi;
}

