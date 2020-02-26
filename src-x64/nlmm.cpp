#define ARMA_DONT_PRINT_ERRORS
/*#include <Rcpp.h>*/
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
double dmvnrm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd = false){
	double out = 0;
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    for (int i = 0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
        out += constants - 0.5 * arma::sum(z%z) + rootisum;     
    }  
      
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List C_wdensgl(NumericMatrix w, NumericVector beta, NumericMatrix Sigma1, List Sigma2, NumericVector alpha, List Y, NumericMatrix x, NumericMatrix z, int M, int N, IntegerVector ni, int P, int Q){

	double val = 0;
	int start = 0;
	double k1 = 1/alpha[0];
	double k2 = 1/alpha[1];
	double dg = 0;
	double dn = 0;
	List Omega(M);
	
	for (int i = 0; i < M; ++i){// loop over subjects
		int n = ni[i];
		
		arma::mat y = Y[i];
		arma::mat sigma2 = Sigma2[i];
		arma::rowvec mu(n,arma::fill::zeros);
		arma::mat zS(n,Q,arma::fill::zeros);
		arma::mat zSz(n,n,arma::fill::zeros);
		arma::mat omega(n,n,arma::fill::zeros);
		
		for(int j = 0; j < n; ++j){			
			// vector x %*% beta
			for(int k = 0; k < P; ++k){
				mu[j] += x(start + j,k)*beta[k];
			}
			// matrix Z %*% Sigma1 %*% t(Z)
			for(int kk = 0; kk < Q; ++kk){
				for(int k = 0; k < Q; ++k){
					zS(j,kk) += z(start + j,k)*Sigma1(k,kk);
				}
			}
			for(int jj = 0; jj < n; ++jj){
				for(int k = 0; k < Q; ++k){
					zSz(j,jj) += zS(j,k)*z(start + jj,k);
				}
			}
			// matrix Omega
			for(int jj = 0; jj < n; ++jj){
				omega(j,jj) = w(i,0)*zSz(j,jj) + w(i,1)*sigma2(j,jj);
			}
		}
		Omega[i] = omega;
		// Multivariate normal density
		dn = dmvnrm_arma(y, mu, omega, true);
		//Rcout << "dn = " << dn << std::endl;
		// Gamma density
		dg = -log(tgamma(k1)) + (k1-1) * log(w(i,0)) - w(i,0) - log(tgamma(k2)) + (k2-1) * log(w(i,1)) - w(i,1);
		//Rcout << "dg = " << dg << std::endl;
		// Log-likelihood
		val += dg + dn;
		// Increment index
		start += n;
	}
	
	List ans;
	ans["val"] = val;
	ans["Omega"] = Omega;
	
	return ans;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double C_ll(NumericMatrix knots, NumericMatrix weights, NumericVector beta, NumericMatrix Sigma1, List Sigma2, List Y, NumericMatrix x, NumericMatrix z, int M, int N, IntegerVector ni, int P, int Q, int K){

	NumericVector val(M);
	double ans = 0;
	
	for(int gq2 = 0; gq2 < K; ++gq2){
		for(int gq1 = 0; gq1 < K; ++gq1){
			int start = 0;
			for (int i = 0; i < M; ++i){// loop over subjects
				int n = ni[i];
				arma::mat y = Y[i];
				arma::mat sigma2 = Sigma2[i];
				arma::rowvec mu(n,arma::fill::zeros);
				arma::mat zS(n,Q,arma::fill::zeros);
				arma::mat zSz(n,n,arma::fill::zeros);
				arma::mat omega(n,n,arma::fill::zeros);

				for(int j = 0; j < n; ++j){			
					// vector x %*% beta
					for(int k = 0; k < P; ++k){
						mu[j] += x(start + j,k)*beta[k];
					}
					// matrix Z %*% Sigma1 %*% t(Z)
					for(int kk = 0; kk < Q; ++kk){
						for(int k = 0; k < Q; ++k){
							zS(j,kk) += z(start + j,k)*Sigma1(k,kk);
						}
					}
					for(int jj = 0; jj < n; ++jj){
						for(int k = 0; k < Q; ++k){
							zSz(j,jj) += zS(j,k)*z(start + jj,k);
						}
					}
					// matrix Omega
					for(int jj = 0; jj < n; ++jj){
						omega(j,jj) = knots(gq1,0)*zSz(j,jj) + knots(gq2,1)*sigma2(j,jj);
					}
				}
				// Integrated log-likelihood
				val[i] += dmvnrm_arma(y, mu, omega, false) * weights(gq1,0) * weights(gq2,1);
				//Rcout << omega << std::endl;
				// Increment index
				start += n;
			} // end loop over subjects
		} // end loop over W_{1}
	} // end loop over W_{2}
	
	for (int i = 0; i < M; ++i){
		ans += log(val[i]);
	}
	
	return ans;
}

