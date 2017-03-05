#ifndef primes_h
#define primes_h

class Prime_table{
public:
  Prime_table() {_pi=NULL; _p=NULL;};
  Prime_table(long x) {init_pi(x); init_p(x);}
  long pi(long x) {return _pi[x];}
  long prime(int i) {return _p[i];}
  long maxprime() {return maxp;}
  long number_of_primes() {return pimax;}
  void init(long x);
  void display();
  ~Prime_table();
private:
  void init_pi(long x);
  void init_p(long x);
  long pimax;
  int *_pi;
  int *_p;
  long maxp;
};

#endif
