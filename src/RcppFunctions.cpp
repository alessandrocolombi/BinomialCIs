// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppGSL.h>

// Include file with basic libraries to include
#include "headers.h"
#include "recurrent_traits.h"

#include "mysample.h"


using namespace Rcpp;

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	log_stable_sum
//------------------------------------------------------------------------------------------------------------------------------------------------------
double log_stable_sum(const std::vector<double>& a, const bool is_log, const double& val_max, const unsigned int& idx_max)
{

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;

	if(is_log==TRUE){ // a contains values in log scale

		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (val_max +
				std::log(1 +
					    std::accumulate(   a.cbegin(), a.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )  +
					    std::accumulate(   a.cbegin()+idx_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )
				        )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;
		// Do not checks if values of a are strictly positive
		return ( std::log(val_max) +
				 std::log(1 +
					      std::accumulate(   a.cbegin(), a.cbegin()+idx_max, 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )  +
					      std::accumulate(   a.cbegin()+idx_max+1, a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   )
			             )
			   );
	}
}

double log_stable_sum(const Rcpp::NumericVector& a, const bool is_log, const double& val_max){

	double inf = std::numeric_limits<double>::infinity();

	if(a.size() == 0)
		return 0.0;

	// Do not checks if it is really the max value
	if(is_log){ // a contains values in log scale

		if(val_max == -inf) // need to handle the case where all values are -inf
			return -inf;

		return (val_max +
					std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(x - val_max );}   )   )
			   );
	}
	else{

		if(val_max < 0)
			throw std::runtime_error("log_stable_sum, if is_log is FALSE, the maximum value can not be negative ");
		if(val_max == 0)
			return 0.0;

		return ( std::log(val_max) +
				 std::log( std::accumulate(   a.cbegin(), a.cend(), 0.0, [&val_max](double& acc, const double& x){return acc + exp(std::log(x) - std::log(val_max) );}   ) )
			   );
	}
}

// In this version of the formula, the maximum value is computed
// [[Rcpp::export]]
double log_stable_sum(const Rcpp::NumericVector& a, const bool is_log){
	if(a.size() == 0)
		return 0.0;

	// Computes maximum value
	auto it_max{std::max_element(a.cbegin(), a.cend())};
	double val_max{*it_max};
	// Calls the specialized version
	return log_stable_sum(a,is_log,val_max);
}



//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Factorials and Pochammer
//------------------------------------------------------------------------------------------------------------------------------------------------------

//' log Raising Factorial (old)
//'
//' This function computes the logarithm of the rising factorial \code{(a)_n} implementing it from scratch.
//' Notation is log( Gamma(a+n)/Gamma(a) )
// [[Rcpp::export]]
double log_raising_factorial_old(const unsigned int& n, const double& a)
{
	if(n==0)
		return 0.0;
	if(a<0)
		throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale");
	else if(a==0.0){
		return -std::numeric_limits<double>::infinity();
	}
	else{

		double val_max{std::log(a+n-1)};
		double res{1.0};
		if (n==1)
			return val_max;
		for(std::size_t i = 0; i <= n-2; ++i){
			res += std::log(a + (double)i) / val_max;
		}
		return val_max*res;

	}
}

//' log Raising Factorial
//'
//' This function computes the logarithm of the rising factorial \code{(a)_n} implemented
//' as the difference of lgamma functions
//' Notation is log( Gamma(a+n)/Gamma(a) )
// [[Rcpp::export]]
double log_raising_factorial(const unsigned int& n, const double& a)
{
	if(n < 0)
		throw std::runtime_error("Error in log_raising_factorial: n can not be negative ");
	if(n==0)
		return 0.0;
	if(a<0)
		throw std::runtime_error("Error in log_raising_factorial, can not compute the raising factorial of a negative number in log scale");
	else if(a==0.0){
		return -std::numeric_limits<double>::infinity();
	}
	else{
		return std::lgamma((double)n + a) - std::lgamma(a);
	}
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Riemann and Hurwitz zeta functions
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double RiemannZeta(const double& a)
{
	double inf = std::numeric_limits<double>::infinity();
	if(a <= 0)
		return inf;
	else 
		return gsl_sf_zeta(a);
}

// [[Rcpp::export]]
double HurwitzZeta(const double& a, const unsigned int& m)
{
	double inf = std::numeric_limits<double>::infinity();
	if(a <= 0)
		return inf;
	else 
		return gsl_sf_hzeta(a, m);
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Features - Frequentist
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double compute_UBFreq_BeBe(const int& n, const double& alpha_lev, const double& beta, const double& Shat)
{
	if(n <= 0)
		throw std::runtime_error("Error in compute_log_UBFreq_BeBe: n must be strictly positive ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_log_UBFreq_BeBe: alpha_lev must be in (0,1) ");
	if(alpha_lev-beta <= 0 )
		throw std::runtime_error("Error in compute_log_UBFreq_BeBe: alpha_lev-beta must be strictly positive");
	if(1.0/beta <= 0 )
		throw std::runtime_error("Error in compute_log_UBFreq_BeBe: 1/beta must be strictly positive");
	if(Shat < 0)
		throw std::runtime_error("Error in compute_log_UBFreq_BeBe: Shat must be positive ");

	double res{1.0/(double)n};
	// Mod 17/09/25. Ci eravamo persi un quadrato?
	double Warg = (double)n/(alpha_lev - beta) * 
					(  std::sqrt( (std::log(1.0/beta))/(2.0*(double)n) ) + std::sqrt( (std::log(1.0/beta))/(2.0*(double)n) + Shat ) ) *
					(  std::sqrt( (std::log(1.0/beta))/(2.0*(double)n) ) + std::sqrt( (std::log(1.0/beta))/(2.0*(double)n) + Shat ) ) ;
	double temp{gsl_sf_lambert_W0(Warg)};
	
	if(temp < 0)
		throw std::runtime_error("Error in compute_log_UBFreq_BeBe: the argument of the Lambert W function can not be negative");

	res *= temp;
	return res;				
}


