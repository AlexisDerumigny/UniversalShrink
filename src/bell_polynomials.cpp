#include <Rcpp.h>
using namespace Rcpp;


const int MAXN = 100;
long long C[MAXN + 1][MAXN + 1];

void precompute_combinations() {
  for (int n = 0; n <= MAXN; n++) {
    C[n][0] = C[n][n] = 1;
    for (int r = 1; r < n; r++) {
      C[n][r] = C[n-1][r-1] + C[n-1][r];
    }
  }
}

// long long fact[MAXN + 1];
// 
// void precompute_factorials() {
//   fact[0] = 1;
//   for (int i = 1; i <= MAXN; i++)
//     fact[i] = fact[i - 1] * i;
// }


// [[Rcpp::export]]
NumericMatrix bellPolynomials(NumericVector x, int verbose) {
  
  precompute_combinations();
  
  int m = x.length();
  
  // Elements of this matrix are B(n,k)
  NumericMatrix B(m + 1, m + 1);
  
  B(0,0) = 1;
  
  for (int n = 1; n <= m; n++) {
    B(n, 0) = 0;
    B(0, n) = 0;
    for (int k = 1; k <= n; k++) {
      double temp = 0;
      for (int j = 1; j <= n - k + 1; j++) {
        double temp_add = C[n-1][j-1] * x(j - 1) * B(n-j, k-1);
        if (verbose > 0) {
          Rcout << n << " " << k << " " << j << " " << temp_add << "\n";
          Rcout << C[n-1][j-1] << " " << x(j - 1) << " " << B(n-j, k-1) << "\n";
        }
        temp = temp + temp_add;
      }
      B(n, k) = temp;
    }
  }
  
  return B;
}



/*** R
bellPolynomials(c(1:20))
*/


