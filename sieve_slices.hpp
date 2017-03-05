#ifndef sieve_by_slice_hpp
#define sieve_by_slice_hpp
//#include "sieve_slices.h"
#include"primes.h"
extern int tiny_primes[];

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::create(long slice_size, 
					longint startx, sieve_type t)
{
  sieve_t= t;
  window_size=next_mult<longint>(slice_size,_period );
  window_start=(startx/window_size)*window_size;
  window_end=window_start+window_size-1;
  btable::create(1 + lower_index64(window_size));
  btable::unset_bit(0);
  last_total = 0;
  if (sieve_t == AUTO_SIEVE)
    {
#ifdef DEBUG_SV
      cout << "sieve_by_slice<btable>::create: AUTO_SIEVE is set: Eratosthenes_sieve will be done\n";
#endif
      eratosthenes();
    }
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::set_around(longint x)
{
  //cout << "set_around x= " << x << "  window_start= " << window_start <<  endl;
  if ((x < window_start) || (x > window_end)) {
    long q = x / window_size;
    window_start = q * window_size;
    window_end   = window_start + window_size-1;
    btable::fill();
    btable::unset_bit(0);
    if (sieve_t == AUTO_SIEVE)
      {
#ifdef DEBUG_SV
	cout << "sieve_by_slice<btable>::create: AUTO_SIEVE is set: Eratosthenes_sieve will be done\n";
#endif
	eratosthenes();
      }
  }
}


template<class btable, class longint> void 
sieve_by_slice<btable, longint>::sieve_by(long p)
{
  //cout << "sieve by " << p << endl;
  long inc = p<<1;
  long true_inc = (p<<2) + (p<<1);
  
  int r  = p % 6;
  longint start;
  long index;
  
  if (r==1)
    {
      // j=0
      start = p;
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);

      // j=1
      start= p+(p<<2);   // 5p
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);
    }
  else // r=5
    {
      //j=0;
      start= p+(p<<2);
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);

      //j=1
      start= p;
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);
    }
}


template<class btable, class longint> void 
sieve_by_slice<btable, longint>::shift_window_forward()
{

#ifdef DEBUG_SV
  cout << "In sieve_slice::shift_window_forward forward start = " << window_start << endl;
  cout << "                                           end  = " << window_start+window_size-1 << endl;
#endif
  
  if (btable::init_counters()) {
    if (sieve_t==AUTO_SIEVE)
      last_total = count(window_start+window_size-1);
  }
  window_start += window_size;
  window_end += window_size;
  btable::fill();
  btable::unset_bit(0);
  if (sieve_t == AUTO_SIEVE)
    {
#ifdef DEBUG_SV
      cout << "sieve_slice::forward, AUTO_SIEVE is set: Eratosthenes sieve will be done" << endl;
#endif
      eratosthenes();
    }
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::shift_window_backward()
{
  window_start -= window_size;
  window_end   -= window_size;
  btable::fill();
  btable::unset_bit(0);
  if (sieve_t == AUTO_SIEVE)
    {
      eratosthenes();
      last_total -= btable::count((int)lower_index64(window_size-1));  
    }
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::display(int how_many)
{
  cout << "Sieve_linear::display size= " << window_size << "   [" << window_start << ", " << window_start + window_size - 1 << "]"  \
       << "    sieve_type= " << sieve_t;
  cout << "    last_total= " << last_total << endl;
  cout << "    btable::get_bit_size= " << btable::get_bit_size() << endl;
  int cnte=0;
  for (int i = 1; i < btable::get_bit_size(); i++)
    if (btable::get_bit(i))
      {
	cnte++;
	cout.width(8);
	cout << get_integer(i);
	if (cnte==how_many)
	  break;
      }
  cout << "\n\n";
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::display_counts() {
    cout << "Sieve_linear::display(): Les comptes\n";
    for ( int i = 0 ; i < btable::bit_size ; i++)
      {
	cout.width(8);
	cout << count(get_integer(i));
      }
    cout << "\n\n";
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::eratosthenes()
{
#ifdef DEBUG_SV
  cout << "In ERATOS\n";
#endif
  if (sieve_t == NO_SIEVE)
    {
      cout << "Eratosthenes is not defined for this bit table!\n";
      exit(0);
    }
  long p;
  //long ip =1;
  long ip =3;

  if (!window_start) btable::unset_bit(1);
  double bound_p = sqrt((double)get_window_end()+2);
  double bound1_p = min(double(maxprime()),bound_p);
  //cout << "bound1_p= " << bound1_p << endl;
  while (p = prime(ip++), p <= bound1_p)
   {
     sieve_by(p);
     if (p >= window_start)
       btable::set_bit(lower_index64(p-window_start));
   }

  if (bound_p > bound1_p) {
    cerr << "maxprime() = " << maxprime()\
	 << "  is too small to sieve until window_end = " << get_window_end() << endl;
    cerr << "We finish eratosthenes_sieve using a prime_generator " << endl;
    prime_generator pg(maxprime(),long(min(1.1*bound_p, 4000000000.0)), maxprime());
    p=pg.next_prime();
    while (p < bound_p) {
      sieve_by(p);
      if (p >= window_start) 
	btable::set_bit(lower_index64(p-window_start)); 
      p=pg.next_prime();
    }
  }

  btable::init_counters();
  left_index = 0;
  right_index = btable::get_bit_size();
#ifdef DEBUG_SV
  cout << "OUT ERATOS\n";
#endif
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_first_prime()
{
  left_index=0;
  for (;;)
    {
      while(++left_index < btable::get_bit_size())
	{
	  if (btable::get_bit(left_index))
	    return get_integer(left_index);
	}
      shift_window_forward();
    }
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_next_prime()
{
  //cout << "In get_next_prime left_index= " << left_index << endl;
    if (left_index < 0)
    switch (left_index) {
    case -1:
      left_index++;
      return 3;
    case -2:
      left_index++;
      return 2;
    }

    //cout << "sieve_t = " << sieve_t << endl;
    //cout << "wstart= " << get_window_start() << "    wend= " << get_window_end() << endl;

    if (sieve_t == NO_SIEVE)
      {
	cout << "sieve_by_slice de type NO_SIEVE\n";
	cout << "sieve_by_slice::get_next_prime n'est pas défini\n";
	exit(0);
      }
    for (;;)
      {
	while(++left_index < btable::get_bit_size())
	  {
	    if (btable::get_bit(left_index))
	      return get_integer(left_index);
	  }
	shift_window_forward();
      }
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_next_prime_without_shifting()
{
  for (;;)
    {
      while(++left_index < btable::get_bit_size())
	{
	  if (btable::get_bit(left_index))
	    return get_integer(left_index);
	}
      return 0;
    }
}


template<class btable, class longint> longint  
sieve_by_slice<btable, longint>::get_previous_prime()
{
  while (right_index >= 0) {
    while ((--right_index) >= 0)
      {
	if (btable::get_bit(right_index)) {
	  return get_integer(right_index);
	}
      }
    if (window_start) 
      {
	shift_window_backward();
      }
  }

  if (right_index>=-NUMBER_OF_PRESIEVED_PRIMES-1) {
    right_index -= 1;
    return tiny_primes[NUMBER_OF_PRESIEVED_PRIMES+right_index+2];
  }
  return 0;
}

template<class btable, class longint> void sieve_by_slice<btable, longint>::set_indexes(longint x)
{
  //cout << "In set_indexes x= " << x << endl;
  set_around(x);
  
  //  if ((x < 5) && (window_start==0)){
  if (x < 5) {
    switch (x) {
    case 4:
      right_index = -1;
      left_index  =  0;
      break;

    case 3:
      right_index =  -2; 
      left_index  =  0;
      break;

    case 2:
      right_index =  -3; 
      left_index  = -1;
      break;

    case 1:
    case 0:
      right_index = -3;
      left_index  = -2;
      break;
    }
    //cout << "left_index set to " << left_index << endl;
    //cout << "right_index set to " << right_index << endl;
    return;
  }
 
  left_index = lower_index64(x - window_start);

  right_index = (int)lower_index64(x - window_start);
  if (get_integer(right_index) < x)
    right_index++;
}


template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_previous_prime(longint x)
{
  set_around(x);
  right_index = (int)lower_index64(x - window_start);
  if (get_integer(right_index) < x)
    right_index++;
  return get_previous_prime();
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_next_prime(longint x)
{
  left_index = (int)lower_index64(x - window_start);
  return get_next_prime();
}


template<class btable, class longint> longint
sieve_by_slice<btable, longint>::count(longint x)
{
  if (sieve_t == AUTO_SIEVE)
    while (x >= window_start + window_size) 
      {
	shift_window_forward();
      }

  return(last_total + 
	 btable::count((int)lower_index64(x-window_start))); 
}

template<class btable, class longint> int
sieve_by_slice<btable, longint>::belong_to_window(longint x)
{
  return (x >= window_start) && (x < window_start + window_size);
}
#endif
