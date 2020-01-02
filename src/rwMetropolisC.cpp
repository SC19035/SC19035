#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector rwMetropolisC(double sigma, double x0, double N) {
	NumericVector x(N+1);
	x[0] = x0;
	NumericVector u = runif(N);
	double k = 0;
	double y = 0;
	for (int i = 2; i < N+1; i++) {
		y = rnorm(1, x[i-2], sigma)[0];
		if (u[i-2] <= exp(-((abs(y))-(abs(x[i-2])))) ) {
			x[i-1] = y;
		} 
		else {
			x[i-1] = x[i-2];
			k++;
		}
	}
	x[N] = k;
	return(x);
}
















