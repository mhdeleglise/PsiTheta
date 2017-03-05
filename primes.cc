#include<cstdlib>
#include<iostream>
#include"./basic_templates.h"
#include"./primes.h"
using namespace std;



void Prime_table::init_pi(long x)
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


void Prime_table::init_p(long x)
{
  maxp=0;
  int i,pi;
  _p = new int[pimax+2];
  if (_p==NULL) {cerr << "no memory for _p !" << endl; exit(1);}
  _p[0] = 0;
  pi = 0;
  for (i=2;i<=x;i++) 
    if (_pi[i]-_pi[i-1]) _p[++pi] = i;
  maxp = _p[pimax];
  _p[pimax+1]=2147483647;
}


void Prime_table::init(long x)
{
    init_pi(x);
    init_p(x);
}

void Prime_table::display()
{
  cout << "Prime_table number_of_primes() = " << number_of_primes() << "   maxprime= " << maxprime() << endl;
}

Prime_table::~Prime_table() {
  delete [] _p;
  delete [] _pi;
}
