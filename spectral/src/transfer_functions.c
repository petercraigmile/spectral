
#include <Rmath.h>


void squared_gain (double *sq_gain,
		   double *fs, int *N,
		   double *mult,		   
		   double *filter, int *L) {

  int i, l;
  double w, cos_term, sin_term, cw, sw;
  double sin_lw, cos_lw, old_sin_lw, old_cos_lw;
  double K;

  K = 2.0 * M_PI * (*mult);

  for (i=0; i<(*N); i++)  {

    w = K * fs[i];

    cw = cos(w);
    sw = sin(w);

    cos_lw = 1.0;
    sin_lw = 0.0;

    cos_term = filter[0];
    sin_term = 0.0;

    for (l=1; l<(*L); l++) {

      old_cos_lw = cos_lw;
      old_sin_lw = sin_lw;

      cos_lw = old_cos_lw * cw - old_sin_lw * sw;
      sin_lw = old_sin_lw * cw + old_cos_lw * sw;

      cos_term += filter[l] * cos_lw;
      sin_term += filter[l] * sin_lw;
    }

    sq_gain[i] = cos_term * cos_term + sin_term * sin_term;
  }
}
  


