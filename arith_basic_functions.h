#ifndef arith_basic_functions_h
#define arith_basic_functions_h
#include<math.h>
#include<cmath>
#include<cstdlib>
#include<limits.h>
#include<stdint.h>

typedef __int128_t  long128;
std::ostream& operator<<(std::ostream& stream, const long128& x);
long128 atolong128(char* s, int base=10);
double pow128(long128 x, double e);

const int _number_of_classes=2;
const int _period=6;
extern int _class[2];
extern int _index[6];

inline long value(long i)
{
  if (!i) return 0;

  long q=(i-1) / 2;
  long r=(i-1) %2;

  return  + q*_period + _class[r];
}


struct ldiv64 { 
  long q;
  long r; 
  inline ldiv64(long a);// Quotient et reste de la division par 6
};

class ldiv128 { 
public:
  long128 q;
  long r; 
  inline ldiv128(long128 a); // Quotient et reste de la division par 6
};

inline long lower_index64(long x)
{
#ifdef DEBUG_SAFE
  if (x < 0)
    {
      cout << "lower_index: incorrect x = " << x << "\n";
      error();
    }
#endif 
  ldiv64 qr(x);
  return 2*(qr.q) + _index[(int)qr.r];
}

inline long128 lower_index128(long128 x)
{
#ifdef DEBUG_SAFE
  if (x < 0)
    {
      cout << "lower_index: incorrect x = " << x << "\n";
      error();
    }
#endif 
  ldiv128 qr(x);
  return 2*(qr.q) + _index[(int)qr.r];
}


inline int number_bit(unsigned long a)
{
  const unsigned long mask1 = 0x5555555555555555;
  const unsigned long mask2 = 0x3333333333333333;

  a = (a&mask1) + ((a>>1)&mask1);
  a = (a&mask2) + ((a>>2)&mask2);

  unsigned long b = (a&0xFFFFFFFF)+(a>>32);
  b = (b&0x0F0F0F0F)+((b>>4)&0x0F0F0F0F);
  b = (b&0xFFFF)+(b>>16);
  b = (b&0xFF)+(b>>8);
  return int(b);
}

long min(long a, long b);

long squareroot128(long128 x);
long squareroot64(long x);
int cuberoot128(long128 x);
double log128(long128 x);
double exp128(long128 x);
double pow128(long128 x, double e);

template<class longint> longint prev_mult(longint a, int periode) {
  return a - a%periode;
}

template<class longint> longint next_mult(longint b, int periode) {
  int r= b % periode;
  return (r==0) ? b : b - b%periode + periode;
}

// Quotient et reste de la division par 6
inline ldiv64::ldiv64(long a)
{
  long y=a;
  q=y;
  while (y > 0) {
    y >>= 2;
    q+=y;
  }
  q >>= 3;
  long q6 = (q<<1) + (q<<2);
  r= a-q6;
  while (r >= 6) {
    q+=1;
    r-=6;
  }
}

inline ldiv128::ldiv128(long128 a)
{
  long128 y=a;
  q=y;
  while (y > 0) {
    y >>= 2;
    q+=y;
  }
  q >>= 3;
  long q6 = (q<<1) + (q<<2);
  r= a-q6;
  while (r >= 6) {
    q+=1;
    r-=6;
  }
}

inline int compare_int(int a, int b)
{
  if (a == b) return 0;
  return (a < b) ? -1 : 1;
}

inline int compare_abs_int(int a, int b)
{
  int abs_a = abs(a);
  int abs_b = abs(b);
  if (abs_a == abs_b) return 0;
  return (abs_a < abs_b) ? -1 : 1;
}

inline int max32(int a, int b)
{
  if (a < b) return b;
  return a;
}

inline double max_double(double a, double b)
{
  if (a < b) return b;
  return a;
}

inline long max64(long a, long b)
{
  if (a < b) return b;
  return a;
}

inline int min32(int a, int b)
{
  if (b < a) return b;
  return a;
}

inline long min64(long a, long b)
{
  if (b < a) return b;
  return a;
}

inline long atolong(char* s)
{
  return atol(s);
}

#endif

