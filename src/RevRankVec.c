#include <R.h> /* for NA_REAL, includes math.h */

/* x, y: two decreasingly ordered vectors of length nx and ny; rv: a
vector of length nx, rv[i] is the reverse rank of x[i] in y (number
of y[j] >= x[i]). This function is useful for computing resampling
based p-values efficiently. */
void vec_rev_rank(double *x, double *y, int *nx, int *ny,
              int *rv)
{
  int i, j, k, count=0;
  for (i=0; i<*nx; i++){
    for (j=count; j<*ny; j++){
      if (x[i] <= y[j]) {
        count += 1;
      } else {
        rv[i] = count; break;
      }
    }
    if (count>=*ny) {
      for (k=i; k<*nx; k++) {
        rv[k] = *ny;
      }
      break;
    }
  }
}
