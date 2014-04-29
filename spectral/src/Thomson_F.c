
#include <R.h>
#include <Rmath.h>


void R_Thomson_F (double *stat, long *index, 
		   double *x, long *N,
		   double *tapers, long *K,
		   double *H0, double *mt_real, double *mt_imag)
{
  double C1_real = 0.0, C1_imag = 0.0;
  double mtr = 0.0, mti = 0.0;
  double SS = 0.0, RSS = 0.0, n1 = (double)(*K-1);
  double a1, a2, mult;
  long   k, r, t;

  mult = 2.0 * M_PI * (double)(*index) / (double)(*N);

  for (k=0, r=0; k<(*K); k++)
  {
    SS += H0[k]*H0[k];

    mtr = mti = 0.0;
    for (t=0; t<(*N); t++, r++)
    {
      a1   = tapers[r] * x[t];
      a2   = (double)t * mult;
      mtr += a1 * cos(a2);
      mti -= a1 * sin(a2);
    }
    mt_real[k] = mtr;
    mt_imag[k] = mti;

    C1_real += mtr * H0[k];
    C1_imag += mti * H0[k];
  }
  
  C1_real /= SS;
  C1_imag /= SS;

  for (k=0; k<(*K); k++)
  {
    a1 = mt_real[k] - C1_real * H0[k];
    a2 = mt_imag[k] - C1_imag * H0[k];
    RSS += (a1*a1 + a2*a2);
  }
  
  *stat = n1 * (C1_real*C1_real + C1_imag*C1_imag) * SS / RSS;
}

