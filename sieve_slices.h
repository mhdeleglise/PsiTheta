#ifndef sieve_by_slice_h
#define sieve_by_slice_h
#include<cstdlib>
#include<iostream>
#include "bit_table.h"
#include "arith_basic_functions.h"
#include "primes.h"

const int NUMBER_OF_PRESIEVED_PRIMES=2;

template<class btable, class longint> class sieve_by_slice:
  public Prime_table,
  public btable
  {
    public:
    longint window_start;
    longint window_end;
    long window_size;
    
    long left_index;
    long right_index;
    void create(long window_size, longint startx, sieve_type t);
    sieve_by_slice(long maxp): Prime_table(maxp) {};
    sieve_by_slice(long maxp, long window_size, longint startx, sieve_type t):
      Prime_table(maxp) { create(window_size, startx, t);}
    void set_around(longint x);
    void display(int how_many = 0);
    void display_counts();
    ~sieve_by_slice() {};

    long index(long128 u) {return lower_index64(u-window_start);}
    long index(long u) {return lower_index64(u-window_start);}
    void sieve_by(long p);
    virtual void eratosthenes();
    int get_window_size() const{return window_size; }
    longint get_window_start() const {return window_start;}
    longint get_window_end() const { return window_start + window_size - 1; } 
    longint get_next_window_start() const { return window_start + window_size;}
    longint get_integer(long i) const{
      return window_start+value(i); }
    int belong_to_window(longint x);
    void shift_window_forward();
    void shift_window_backward();
    longint count(longint x);
    longint get_first_prime();
    longint get_next_prime();
    longint get_next_prime_without_shifting();
    longint get_previous_prime();
    longint get_next_prime(longint x);
    longint get_previous_prime(longint x);
    void set_indexes(longint x);
    // set index first prime at x position so that next_prime will be the first prime
    // bigger than x and previous_prime the last prime smaller than x
    protected:
    sieve_type sieve_t;
    longint last_total;  
};


class prime_generator: public sieve_by_slice<bit_table, long> {

public:
  prime_generator(long maxp, long ws) : sieve_by_slice(maxp, ws,0,AUTO_SIEVE) {set_indexes(1);}

  prime_generator(long maxp, long ws, long starting_from): sieve_by_slice(maxp,ws,starting_from,AUTO_SIEVE) {
    set_indexes(starting_from); // next_prime() will be the first prime bigger than starting_from
  }

  long next_prime() {return sieve_by_slice<bit_table, long>::get_next_prime();}
  long prev_prime() {return sieve_by_slice<bit_table, long>::get_previous_prime();}
  long next_prime(long x) {return sieve_by_slice<bit_table, long>::get_next_prime(x);}
  long prev_prime(long x) {return sieve_by_slice<bit_table, long>::get_previous_prime(x);}
  
  void display() {
    cout << "Prime generator of size " << sieve_by_slice<bit_table, long>::window_size\
	 << ".   Primes will start from = " << sieve_by_slice::get_integer(left_index) << endl;
    cout << "window_start= " << sieve_by_slice::window_start << endl;
    sieve_by_slice::display();
  }

};

#include"sieve_slices.hpp"

#endif
