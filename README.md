# Computing psi(x) and theta(x) Tchebychev functions in time  O(x^(2/3+epsilon))

theta function is defined by theta(x) = sum_{p <= x, p prime } log p

psi is defined by  psi(x) = sum_{p^n <=x, p prime} log p

Let us define delta(x) by
delta(x) = sum_{p^n <= x, p prime, n >= 2} log p
wich is easily computed using O(x^(1/2+epsilon)) time.

theta(x) is computed using formula theta(x) = psi(x) - delta(x)

# Howto use theses files
After downloading the sources in a directory, open a command window
and execute

> make

We will get the two executable programms tst_psi  and   tst_theta

The command

> tst_theta 1000000000000 100

computes theta(1000000000000) with 100 decimal places, rounded to the
nearest value. If you omit the number of decimal places, the default
value is 25.

Examples:

> tst_theta 1234567890

gives

> 1.234518946373189641192934e9

While

> tst_theta 1234567890 100

gives

> 1.234518946373189641192934340523247504502720741909842124778080972936152059737170653995781075599578447e9

It works on a PC using linux Ubunto or a macbook. 
GNU MP and GNU MPFR must be installed on this computer.

# Reférences

M. Deléglise and J. Rivat, 
Computing Psi(x), 
Mathematics  of  Computation.
Volume 67, Number 224, October 1998, Pages 1691-1696.
