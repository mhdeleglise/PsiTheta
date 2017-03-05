// -*- C++ -*-
#ifndef bit_table_h
#define bit_table_h
#include<vector>
#include"arith_basic_functions.h"

using namespace std;

#define WORDSIZE1 63
#define WORDSIZE2  6


//enum count_type{NO_COUNT, COUNT, DYNAMIC};
enum sieve_type{NO_SIEVE, AUTO_SIEVE};


class bit_table {
public:
  void create(long);
  bit_table() {};
  bit_table(long n) { create(n);}
  void display(long n=0) const;
  long bit_size; 
  inline long get_bit_size() const {return bit_size; }
  inline unsigned long get_bit(long) const;
  virtual int init_counters() {return 0;};
  long count(long x) const {return 0;}
  inline void unset_bit(long);
  inline void set_bit(long);
  //virtual ~bit_table();

private:
  //  unsigned long* words_array;
  vector<unsigned long> words_array;
  inline unsigned long get_bit_mask(long n) const { 
    return ((unsigned long)1) << (n & WORDSIZE1); }

protected:
  long last_cnte;
  long word_size;
  unsigned long get_word(long x) const { 
    return words_array[x >> WORDSIZE2];} 
  inline unsigned long get_bit_mask2(long n) const { 
    return ((unsigned long)2) << (n & WORDSIZE1); }
   inline int number_bit_at_right(long x) const {
    return number_bit(get_word(x) & ~(get_bit_mask2(x)-1));}
  inline int get_word_index(long x) const { 
    return x >> WORDSIZE2;}
  inline long get_word_size() const {
    return word_size;}
  inline unsigned long get_nth_word(long n)const { 
    return words_array[n];}
  void delete_words_array();
  void fill();
};

inline void 
bit_table::unset_bit(long x)
{
#ifdef DEBUG_SAFE
  if ((x < 0) || (x >= bit_size))
    {
      cout << "In bit_table::unset_bit: incorrect x = " << x << "\n";
      error();
    }
#endif
  words_array[get_word_index(x)] &= ~get_bit_mask(x);
}

inline void 
bit_table::set_bit(long x)
{
#ifdef DEBUG_SAFE
  if ((x < 0) || (x >= bit_size))
    {
      cout << "In bit_table::set_bit: incorrect x = " << x << "\n";
      error();
    }
#endif
  words_array[get_word_index(x)] |= get_bit_mask(x);
}

inline unsigned long 
bit_table::get_bit(long x) const
{
#ifdef DEBUG_SAFE
  if (x >= bit_size)
    {
      cout << "In bit_table::get_bit: incorrect x = " << x << "\n";
      error();
    }
#endif
  return words_array[get_word_index(x)] & get_bit_mask(x);
}
#endif

