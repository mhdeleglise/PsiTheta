#ifndef basic_templates_h
#define basic_templates_h

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
#endif
