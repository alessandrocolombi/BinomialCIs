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
	if(n==0)
		return 0.0;
	if(a<0)
		throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale");
	else if(a==0.0){
		return -std::numeric_limits<double>::infinity();
	}
	else{
		return std::lgamma((double)n + a) - std::lgamma(a);
	}
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Stick-Breaking
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix r_SB( const int& Nrep, const int& Natoms, 
					const double& alpha, const double& theta, 
					const int& seed)
{
	sample::rbeta beta; 
	sample::GSL_RNG engine(seed);

	if(Nrep <= 0)
		throw std::runtime_error("Error in r_SB: Nrep must be positive");
	if(Natoms < 1)
		throw std::runtime_error("Error in r_SB: Natoms must be at least one");
	if(alpha < 0)
		throw std::runtime_error("Error in r_SB: alpha must be positive");
	if(theta < -alpha)
		throw std::runtime_error("Error in r_SB: theta must be >= alpha");

	NumericMatrix res(Nrep,Natoms);
	for(int b = 0; b < Nrep; b++){
		double cumprod{1.0};
		NumericVector pvec(Natoms);
		for(int i = 0; i < Natoms; i++){
			//Rcpp::Rcout<<"("<<b<<","<<i<<"); cumprod = "<<cumprod<<std::endl;
			double pi{ beta(engine, 1.0-alpha, theta + double(i+1)*alpha) };
			pvec[i] =  pi * cumprod ; // save i-th value
			res(b,i) = pi * cumprod ;
			cumprod *= (1.0-pi); // update cumulative sum 
		}
	}

	return res;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Upper Bound - Oracle
//------------------------------------------------------------------------------------------------------------------------------------------------------


NumericMatrix r_SB(const int& Nrep, const int& Natoms, 
									 const double& alpha, const double& theta, 
									 sample::GSL_RNG& engine)
{
	sample::rbeta beta; 

	if(Nrep <= 0)
		throw std::runtime_error("Error in r_SB: Nrep must be positive");
	if(Natoms < 1)
		throw std::runtime_error("Error in r_SB: Natoms must be at least one");
	if(alpha < 0)
		throw std::runtime_error("Error in r_SB: alpha must be positive");
	if(theta < -alpha)
		throw std::runtime_error("Error in r_SB: theta must be >= alpha");

	NumericMatrix res(Nrep,Natoms);
	for(int b = 0; b < Nrep; b++){
		double cumprod{1.0};
		NumericVector pvec(Natoms);
		for(int i = 0; i < Natoms; i++){
			//Rcpp::Rcout<<"("<<b<<","<<i<<"); cumprod = "<<cumprod<<std::endl;
			double pi{ beta(engine, 1.0-alpha, theta + double(i+1)*alpha) };
			pvec[i] =  pi * cumprod ; // save i-th value
			res(b,i) = pi * cumprod ;
			cumprod *= (1.0-pi); // update cumulative sum 
		}
	}

	return res;
}

// PDparams is a (Nexp x 4) matrix whose cols are: alpha, theta, n, Kn
// [[Rcpp::export]]
NumericVector logUB_Oracle( const int& Nexp, const int& Nrep, 
														const IntegerVector& Natoms,
								  					const NumericMatrix& PDparams, 
								  					const double& alpha_lev, const int& seed)
{
	sample::rbeta beta; 
	sample::GSL_RNG engine(seed);
	
	if(Nexp <= 0)
		throw std::runtime_error("Error in logUB_Oracle: Nexp must be positive");
	if(Nrep < 1)
		throw std::runtime_error("Error in logUB_Oracle: Nrep must be at least one");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in logUB_Oracle: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");
	if(Nexp !=  PDparams.nrow())
		throw std::runtime_error("Error in logUB_Oracle: the number of rows of PDparams must be equal to Nexp");
	if(4 != PDparams.ncol())
		throw std::runtime_error("Error in logUB_Oracle: the number of cols of PDparams must be equal to 4");
	if(Natoms.length() != Nexp)
		throw std::runtime_error("Error in logUB_Oracle: Natoms must be a vector of length Nexp");


	NumericVector log_res(Nexp);

	for(int ii = 0; ii < Nexp; ii++)
	{
		double alpha_post{PDparams(ii,0)}; //alpha
		double theta_post{PDparams(ii,1) + PDparams(ii,3)*PDparams(ii,0)}; // theta + K*alpha
		NumericMatrix p_pyp = r_SB(Nrep, Natoms[ii], alpha_post, theta_post, engine);

		double a_beta_post{PDparams(ii,1) + PDparams(ii,3)*PDparams(ii,0)}; // theta + K*alpha
		double b_beta_post{PDparams(ii,2) - PDparams(ii,3)*PDparams(ii,0)}; // n - K*alpha

		NumericVector log_res_vec(Nrep); // log(beta*max(p_pyp))
		for(int bb = 0; bb < Nrep; bb++)
		{
			double log_beta{ std::log( beta(engine,a_beta_post,b_beta_post) ) };
			NumericVector p_pyp_bb( p_pyp(bb, _ ) ); 
			log_res_vec[bb] = log_beta + std::log( Rcpp::max(p_pyp_bb) );

			//Check for User Interruption
			try{
			    Rcpp::checkUserInterrupt();
			}
			catch(Rcpp::internal::InterruptedException e){
			    //Print error and return
			    throw std::runtime_error("Execution stopped by the user");
			}
		}

		// Quantile calculation:
		std::sort(log_res_vec.begin(), log_res_vec.end());
		int idx = std::floor( (1.0 - alpha_lev) * (double)(Nrep - 1) );  // quantile index
		log_res[ii] = log_res_vec[idx];
	}

	return log_res;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Upper Bound - Proposed - Markov
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double compute_log_UBMarkov( const int& Rmax, const double& alpha, const double& theta, const int& Kn, const int& n,
						     const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();

	if(Rmax < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov: Rmax must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in compute_log_UBMarkov: Kn must be smaller or equal to n");
	if(alpha < 0)
		throw std::runtime_error("Error in compute_log_UBMarkov: alpha must be positive");
	if(theta < -alpha)
		throw std::runtime_error("Error in compute_log_UBMarkov: theta must be >= alpha");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in compute_log_UBMarkov: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");

	double log_res{inf};
	for (int r = 1; r <= Rmax; r++)
	{
		// 1/r * ( log(ExpVal) - log(alpha_lev) )
		double log_res_r{ std::log(theta + (double)Kn*alpha) };
		log_res_r += log_raising_factorial(r-1,1.0-alpha) - log_raising_factorial(r,(double)n+theta);
		log_res_r *= 1.0/double(r);
		log_res_r -= std::log(alpha_lev)/double(r);



		// Find minumum
		if(log_res_r < log_res)
			log_res = log_res_r;
	}

	return log_res;
}

// PDparams is a (Nexp x 4) matrix whose cols are: alpha, theta, n, Kn
// [[Rcpp::export]]
NumericVector logUB_Markov( const int& Nexp, const int& Rmax, 
								  					const NumericMatrix& PDparams, 
								  					const double& alpha_lev )
{

	if(Nexp <= 0)
		throw std::runtime_error("Error in logUB_Markov: Nexp must be positive");
	if(Rmax < 1)
		throw std::runtime_error("Error in logUB_Markov: Rmax must be at least one");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in logUB_Markov: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");
	if(Nexp !=  PDparams.nrow())
		throw std::runtime_error("Error in logUB_Markov: the number of rows of PDparams must be equal to Nexp");
	if(4 != PDparams.ncol())
		throw std::runtime_error("Error in logUB_Markov: the number of cols of PDparams must be equal to 4");

	NumericVector log_res(Nexp);

	for(int ii = 0; ii < Nexp; ii++)
	{
		log_res[ii] = compute_log_UBMarkov( Rmax, PDparams(ii,0), PDparams(ii,1), 
																				(int)PDparams(ii,3), (int)PDparams(ii,2), 
																				alpha_lev );
	}

	return log_res;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Upper Bound - Proposed - Cantelli
//------------------------------------------------------------------------------------------------------------------------------------------------------


// [[Rcpp::export]]
double log_ExpMr(const int& r, const double& alpha, const double& theta, const int& Kn, const int& n)
{
	if(r < 1)
		throw std::runtime_error("Error in log_ExpMr: r must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in log_ExpMr: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_ExpMr: Kn must be smaller or equal to n");
	if(alpha < 0)
		throw std::runtime_error("Error in log_ExpMr: alpha must be positive");
	if(theta < -alpha)
		throw std::runtime_error("Error in log_ExpMr: theta must be >= alpha");
	
	// 1/r * ( log(ExpVal) - log(alpha_lev) )
	double log_res_r{ std::log(theta + (double)Kn*alpha) };
	log_res_r += log_raising_factorial(r-1,1.0-alpha) - log_raising_factorial(r,(double)n+theta);
	return log_res_r;
}

// [[Rcpp::export]]
double log_2ndMomMr(const int& r, const double& alpha, const double& theta, const int& Kn, const int& n)
{
	if(r < 1)
		throw std::runtime_error("Error in log_2ndMomMr: r must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in log_2ndMomMr: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_2ndMomMr: Kn must be smaller or equal to n");
	if(alpha < 0)
		throw std::runtime_error("Error in log_2ndMomMr: alpha must be positive");
	if(theta < -alpha)
		throw std::runtime_error("Error in log_2ndMomMr: theta must be >= alpha");

		double log_res1_r{ 0.0 };
		log_res1_r += log_raising_factorial(2*r, theta+(double)Kn*alpha)  - log_raising_factorial(2*r,(double)n+theta);
		
		Rcpp::NumericVector log_vec(2);
		log_vec[0] = log_raising_factorial(2*r-1,1.0-alpha) - log_raising_factorial(2*r-1,theta+(double)Kn*alpha+1.0);
		log_vec[1] = log_raising_factorial(r-1  ,1.0-alpha) - log_raising_factorial(r-1  ,theta+(double)Kn*alpha+alpha+1.0) ;
		log_vec[1] += gsl_sf_lnbeta((double)r-alpha, (double)r + theta + alpha + (double)Kn*alpha + 1.0) - 
									gsl_sf_lnbeta(1.0-alpha, theta + alpha + (double)Kn*alpha + 1.0);
		double log_res2_r = log_stable_sum(log_vec, TRUE);

		return log_res1_r + log_res2_r;
}

// [[Rcpp::export]]
double log_VarMr(const int& r, const double& alpha, const double& theta, const int& Kn, const int& n)
{
	if(r < 1)
		throw std::runtime_error("Error in log_VarMr: r must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in log_VarMr: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_VarMr: Kn must be smaller or equal to n");
	if(alpha < 0)
		throw std::runtime_error("Error in log_VarMr: alpha must be positive");
	if(theta < -alpha)
		throw std::runtime_error("Error in log_VarMr: theta must be >= alpha");

	double lExp  = log_ExpMr(r, alpha, theta, Kn, n);
	double lMom2 = log_2ndMomMr(r, alpha, theta, Kn, n);
	double ldiff = lMom2 - 2.0*lExp; 
	if( ldiff < 0 )
		throw std::runtime_error("Error in log_VarMr: negative variance");
	if( ldiff < 1e-16 )
		throw std::runtime_error("Error in log_VarMr: variance is zero or close to zero");

	// Compute log Variance in log scale (stable operation)
	double log_var{lMom2};
	log_var += std::log( 1.0 - std::exp(-ldiff) );
	return log_var;
}

// [[Rcpp::export]]
double compute_log_UBCantelli( const int& Rmax, 
							   const double& alpha, const double& theta, const int& Kn, const int& n,
							   const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();

	if(Rmax < 1)
		throw std::runtime_error("Error in compute_log_UBCantelli: Rmax must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in compute_log_UBCantelli: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in compute_log_UBCantelli: Kn must be smaller or equal to n");
	if(alpha < 0)
		throw std::runtime_error("Error in compute_log_UBCantelli: alpha must be positive");
	if(theta < -alpha)
		throw std::runtime_error("Error in compute_log_UBCantelli: theta must be >= alpha");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in compute_log_UBCantelli: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");

	double log_res{inf};
	double log_res_r{0.0};
	for (int r = 1; r <= Rmax; r++)
	{

		Rcpp::NumericVector log_vec(2);
		log_vec[0] = log_ExpMr(r, alpha, theta, Kn, n);
		log_vec[1] = 0.5*( std::log(1.0-alpha_lev) - std::log(alpha_lev) + log_VarMr(r, alpha, theta, Kn, n) );
		log_res_r = 1.0/(double)r * log_stable_sum(log_vec, TRUE);


		// Find minumum
		if(log_res_r < log_res)
			log_res = log_res_r;
	}

	return log_res;
}

// PDparams is a (Nexp x 4) matrix whose cols are: alpha, theta, n, Kn
// [[Rcpp::export]]
NumericVector logUB_Cantelli( const int& Nexp, const int& Rmax, 
							  const NumericMatrix& PDparams, 
							  const double& alpha_lev )
{

	if(Nexp <= 0)
		throw std::runtime_error("Error in logUB_Cantelli: Nexp must be positive");
	if(Rmax < 1)
		throw std::runtime_error("Error in logUB_Cantelli: Rmax must be at least one");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in logUB_Cantelli: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");
	if(Nexp !=  PDparams.nrow())
		throw std::runtime_error("Error in logUB_Cantelli: the number of rows of PDparams must be equal to Nexp");
	if(4 != PDparams.ncol())
		throw std::runtime_error("Error in logUB_Cantelli: the number of cols of PDparams must be equal to 4");

	NumericVector log_res(Nexp);

	for(int ii = 0; ii < Nexp; ii++)
	{
		log_res[ii] = compute_log_UBCantelli( Rmax, PDparams(ii,0), PDparams(ii,1), 
																				  (int)PDparams(ii,3), (int)PDparams(ii,2), 
																				  alpha_lev );
	}

	return log_res;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Upper Bound - Painsky
//------------------------------------------------------------------------------------------------------------------------------------------------------

double compute_log_UBPainsky( const int& Rmax, const int& n, const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();

	if(Rmax < 1)
		throw std::runtime_error("Error in compute_log_UBPainsky: Rmax must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in compute_log_UBPainsky: n must be at least one");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in compute_log_UBPainsky: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");

	double log_res{inf};
	for (int r = 1; r <= Rmax; r++)
	{
		// qstar = (r-1)/(r-1+n)
		// 1/r * ( ( qstar^(r-1) * (1-qstar)^(n)  ) / alfa )
		double qstar{ ((double)r - 1.0)/( (double)r - 1.0 + (double)n) };
		double log_res_r{ ( ((double)r - 1.0)*std::log(qstar) + (double)n*std::log(1.0 - qstar)  ) };
		log_res_r *= 1.0/double(r);
		log_res_r -= std::log(alpha_lev)/double(r);

		// Find minumum
		if(log_res_r < log_res)
			log_res = log_res_r;
	}

	return log_res;
}
// [[Rcpp::export]]
NumericVector logUB_Painsky( const int& Nexp, const int& Rmax, const IntegerVector& n_vec, const double& alpha_lev )
{

	if(Nexp <= 0)
		throw std::runtime_error("Error in logUB_Painsky: Nexp must be positive");
	if(Rmax < 1)
		throw std::runtime_error("Error in logUB_Painsky: Rmax must be at least one");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in logUB_Painsky: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");
	if(Nexp != n_vec.length())
		throw std::runtime_error("Error in logUB_Painsky: the length of n_vec must be equal to Nexp");

	NumericVector log_res(Nexp);

	for(int ii = 0; ii < Nexp; ii++)
	{
		log_res[ii] = compute_log_UBPainsky( Rmax, n_vec[ii], alpha_lev );
	}

	return log_res;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Param. Est. - PYP
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double log_eppfPYP( const int& n, const int& Kn, const std::vector<int>& n_j, 
					const double& alpha, const double& theta )
{
	double inf = std::numeric_limits<double>::infinity();

	if(n < 1)
		throw std::runtime_error("Error in log_eppfPYP: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_eppfPYP: Kn must be smaller or equal to n");
	if( std::accumulate(n_j.cbegin(), n_j.cend(), 0) != n )
		throw std::runtime_error("Error in log_eppfPYP: the sum of elements in n_j must be equal to n");
	if(n_j.size() != Kn)
		throw std::runtime_error("Error in log_eppfPYP: the length of n_j must match the number of species Kn");
	std::for_each(n_j.cbegin(), n_j.cend(), [](const int& x){ 
		if(x<=0)
			throw std::runtime_error("Error in log_eppfPYP: invalid number of species");  
	});

	// Check domain
	if(theta < -alpha){
		Rcpp::Rcout<<"Caso proibito (theta < -alpha)"<<std::endl;
		return -std::exp(20);
	}
	if(alpha < 0 || alpha > 1 - 1e-16){
		Rcpp::Rcout<<"Caso proibito (alpha < 0 or alpha > 1)"<<std::endl;
		return -std::exp(20);
	}
	if(alpha < 1e-16 ){
		// Dirichlet process case:
		double res{ (double)Kn*std::log(theta) }; // alpha^Kn
		res -= log_raising_factorial(n, theta); // 1 / (theta)_n
		for(std::size_t j = 0; j < n_j.size(); j++){
			res += lgamma(n_j[j]); // (n_j[j]-1)!
		}
		if( res == inf || std::isnan(res) || res == -inf){
			Rcpp::Rcout<<"Error in log_eppfPYP: NaN, Inf or -Inf returned "<<std::endl;
			Rcpp::Rcout<<"res = "<<res<<std::endl;
			Rcpp::Rcout<<"alpha = "<<alpha<<std::endl;
			Rcpp::Rcout<<"theta = "<<theta<<std::endl;
			Rcpp::Rcout<<"n = "<<n<<std::endl;
			Rcpp::Rcout<<"Kn = "<<Kn<<std::endl;
			Rcpp::Rcout<<"Stampo n_j: ";		
			for(auto __v : n_j)
				Rcpp::Rcout<<__v<<", ";
			Rcpp::Rcout<<std::endl;
			throw std::runtime_error("Error ");
		}
		return res;
	}
	// Poisson-Dirichlet eppf
	double res{ (double)Kn*std::log(alpha) }; // alpha^Kn
	res += log_raising_factorial(Kn, theta/alpha) - log_raising_factorial(n, theta); // (theta/alpha)_Kn / (theta)_n
	for(std::size_t j = 0; j < n_j.size(); j++){
		res += log_raising_factorial(n_j[j] - 1, 1.0 - alpha); // (1-alpha)_(n_j[j]-1)
	}
	if( res == inf || std::isnan(res) || res == -inf){
		Rcpp::Rcout<<"Error in log_eppfPYP: NaN, Inf or -Inf returned "<<std::endl;
		Rcpp::Rcout<<"res = "<<res<<std::endl;
		Rcpp::Rcout<<"alpha = "<<alpha<<std::endl;
		Rcpp::Rcout<<"theta = "<<theta<<std::endl;
		Rcpp::Rcout<<"n = "<<n<<std::endl;
		Rcpp::Rcout<<"Kn = "<<Kn<<std::endl;
		Rcpp::Rcout<<"Stampo n_j: ";		
		for(auto __v : n_j)
			Rcpp::Rcout<<__v<<", ";
		Rcpp::Rcout<<std::endl;
		throw std::runtime_error("Error ");
	}
	return res;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Finite species model
//------------------------------------------------------------------------------------------------------------------------------------------------------

double log_dPois1(const int& x, const double& Lambda)
{
	if(x <= 0)
		return -std::numeric_limits<double>::infinity(); 
	else
		return (-Lambda + (double)(x-1)*std::log(Lambda) - gsl_sf_lnfact(x-1) );
}

// [[Rcpp::export]]
double compute_logV( const int& Kn, const int& n, 
				     const double& gamma, const double& Lambda, 
				     unsigned int M_max )
{
	if(n < 1)
		throw std::runtime_error("Error in compute_log_Vprior: n must be at least one");
	if(Kn < 1)
		throw std::runtime_error("Error in compute_log_Vprior: Kn must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in compute_log_Vprior: Kn must be smaller or equal to n");
	// if here, 1 <= Kn <= n 

	if(gamma < 1e-16)
		throw std::runtime_error("Error in compute_log_Vprior: gamma can not be smaller than 0 (if here, gamma < 1e-16)");
	if(Lambda < 0)
		throw std::runtime_error("Error in compute_log_Vprior: Lambda can not be smaller than 0 ");
	if(M_max < 1 || M_max > 1e5)
		throw std::runtime_error("Error in compute_log_Vprior: invalid number of M_max ");

	// Initialize vector of results
	std::vector<double> log_vect_res(M_max+1, -std::numeric_limits<double>::infinity() );
	// Initialize quantities to find the maximum
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	// Start the loop, let us compute all the elements
	for(int Mstar=0; Mstar <= M_max; ++Mstar){
		log_vect_res[Mstar] = log_raising_factorial( Kn, (double)(Mstar+1) ) +
							  log_dPois1(Mstar+Kn, Lambda) - 
							  log_raising_factorial( n, gamma*(double)(Mstar + Kn) ) ;
		// Check if it is the new maximum
        if(log_vect_res[Mstar]>val_max){
        	idx_max = Mstar;
        	val_max = log_vect_res[Mstar];
        }

	}
	// Formula to compute the log of all the sums in a stable way
	return log_stable_sum(log_vect_res, TRUE, val_max, idx_max);
}

// [[Rcpp::export]]
double log_eppfFD( const int& n, const int& Kn, const std::vector<int>& n_j, 
				   const double& gamma, const double& Lambda, unsigned int M_max )
{
	double inf = std::numeric_limits<double>::infinity();

	// Check n and Kn
	if(n < 1)
		throw std::runtime_error("Error in log_eppfFD: n must be at least one");
	if(Kn < 1)
		throw std::runtime_error("Error in log_eppfFD: Kn must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_eppfFD: Kn must be smaller or equal to n");
	// if here, 1 <= Kn <= n 

	// Check n_j
	if(n_j.size() != Kn)
		throw std::runtime_error("Error in log_eppfFD: the length of n_j must match the number of species Kn");
	std::for_each(n_j.cbegin(), n_j.cend(), [](const int& x){ 
		if(x<=0)
			throw std::runtime_error("Error in log_eppfFD: invalid number of species");  
	});

	// Check parameters
	if(M_max < 1 || M_max > 1e4)
		throw std::runtime_error("Error in log_eppfFD: invalid number of M_max ");
	if(gamma < 1e-16){
		Rcpp::Rcout<<"Caso proibito"<<std::endl;
		return -std::exp(20);
	}
	if(Lambda < 0){
		Rcpp::Rcout<<"Caso proibito"<<std::endl;
		return -std::exp(20);
	}

	
	double logV{ compute_logV( Kn, n, gamma, Lambda, M_max ) }; 
	double res{logV};
	for(std::size_t j = 0; j < n_j.size(); j++){
		res += log_raising_factorial(n_j[j], gamma); 
	}
	if(res > 0){
		// increase M_max and repeat
		log_eppfFD( n, Kn, n_j, gamma, Lambda, 2*M_max );

		//Check for User Interruption
		try{
		    Rcpp::checkUserInterrupt();
		}
		catch(Rcpp::internal::InterruptedException e){
		    //Print error and return
		    throw std::runtime_error("Execution stopped by the user");
		}

	}
	if( res == inf || std::isnan(res) || res == -inf){
		Rcpp::Rcout<<"Error in log_eppfFD: NaN, Inf or -Inf returned "<<std::endl;
		Rcpp::Rcout<<"res = "<<res<<std::endl;
		Rcpp::Rcout<<"logV = "<<logV<<std::endl;
		Rcpp::Rcout<<"gamma = "<<gamma<<std::endl;
		Rcpp::Rcout<<"Lambda = "<<Lambda<<std::endl;
		Rcpp::Rcout<<"n = "<<n<<std::endl;
		Rcpp::Rcout<<"Kn = "<<Kn<<std::endl;
		Rcpp::Rcout<<"Stampo n_j: ";		
		for(auto __v : n_j)
			Rcpp::Rcout<<__v<<", ";
		Rcpp::Rcout<<std::endl;
		throw std::runtime_error("Error ");
	}
	return res;
}


double log_qMpost( const int& m, const int& n, const int& Kn, 
				   const double& gamma, const double& Lambda, 
				   const double& logV,
				   unsigned int M_max )
{
	if(m < 0)
		throw std::runtime_error("Error in log_qMpost: m must be positive or zero");
	if(n < 1)
		throw std::runtime_error("Error in log_qMpost: n must be at least one");
	if(Kn < 1)
		throw std::runtime_error("Error in log_qMpost: Kn must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_qMpost: Kn must be smaller or equal to n");
	// if here, 1 <= Kn <= n 

	if(gamma < 1e-16)
		throw std::runtime_error("Error in log_qMpost: gamma can not be smaller than 0 (if here, gamma < 1e-16)");
	if(Lambda < 0)
		throw std::runtime_error("Error in log_qMpost: Lambda can not be smaller than 0 ");
	if(M_max < 1 || M_max > 1e5)
		throw std::runtime_error("Error in log_qMpost: invalid number of M_max ");

	//double logV{ compute_logV( Kn, n, gamma, Lambda, M_max ) }; 
	double res{ 0.0 };
	res -= logV;
	res += log_raising_factorial( Kn, (double)(m+1) ) +
		   log_dPois1(m + Kn, Lambda) - 
		   log_raising_factorial( n, gamma*(double)(m + Kn) ) ;

    return res;
}



// [[Rcpp::export]]
double log_qMpost( const int& m, const int& n, const int& Kn, 
				   const double& gamma, const double& Lambda, 
				   unsigned int M_max )
{
	if(m < 0)
		throw std::runtime_error("Error in log_qMpost: m must be positive or zero");
	if(n < 1)
		throw std::runtime_error("Error in log_qMpost: n must be at least one");
	if(Kn < 1)
		throw std::runtime_error("Error in log_qMpost: Kn must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_qMpost: Kn must be smaller or equal to n");
	// if here, 1 <= Kn <= n 

	if(gamma < 1e-16)
		throw std::runtime_error("Error in log_qMpost: gamma can not be smaller than 0 (if here, gamma < 1e-16)");
	if(Lambda < 0)
		throw std::runtime_error("Error in log_qMpost: Lambda can not be smaller than 0 ");
	if(M_max < 1 || M_max > 1e5)
		throw std::runtime_error("Error in log_qMpost: invalid number of M_max ");

	double logV{ compute_logV( Kn, n, gamma, Lambda, M_max ) }; 
    return  log_qMpost( m, n, Kn, gamma, Lambda, logV, M_max );
}



// [[Rcpp::export]]
double log_ExpMr_FD(const int& r, const double& gamma, const double& Lambda, const int& Kn, const int& n, unsigned int M_max)
{
	double inf = std::numeric_limits<double>::infinity();

	if(r < 1)
		throw std::runtime_error("Error in log_ExpMr_FD: r must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in log_ExpMr_FD: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in log_ExpMr_FD: Kn must be smaller or equal to n");
	if(gamma < 0)
		throw std::runtime_error("Error in log_ExpMr_FD: gamma must be positive");
	if(Lambda < 0)
		throw std::runtime_error("Error in log_ExpMr_FD: Lambda must be positive");
	
	// Compute qMstar
	/*
	double logV{ compute_logV( Kn, n, gamma, Lambda, M_max ) }; 
	int Mstar_max{500};
	bool sum_not_one{TRUE};
	std::vector<double> lqMpost;
	while( sum_not_one && Mstar_max < 10000){
		//Rcpp::Rcout<<"Mstar_max = "<<Mstar_max<<std::endl;
		double cumsum{0.0};
		lqMpost.resize(Mstar_max);
		for(int m = 0; m < Mstar_max; m++){
			lqMpost[m] = log_qMpost( m, n, Kn, gamma, Lambda, logV, M_max );
			cumsum += std::exp(lqMpost[m]);
		}
		//Rcpp::Rcout<<"qMpost somma a "<<cumsum<<std::endl;
		if( (1.0 - cumsum) < 1e-4  )
			sum_not_one = FALSE;
		else
			Mstar_max += 500;

		//Check for User Interruption
		try{
		    Rcpp::checkUserInterrupt();
		}
		catch(Rcpp::internal::InterruptedException e){
		    //Print error and return
		    throw std::runtime_error("Execution stopped by the user");
		}

	}
	*/
	// Compute qMstar
	double logV{ compute_logV( Kn, n, gamma, Lambda, M_max ) }; 
	std::vector<double> lqMpost;
	lqMpost.reserve(500);
	double cumsum{0.0};
	int m{0};
	bool sum_not_one{TRUE};
	while( sum_not_one && m < 10000){
		lqMpost.push_back( log_qMpost( m, n, Kn, gamma, Lambda, logV, M_max ) );
		cumsum += std::exp(lqMpost[m]);

		if( (1.0 - cumsum) < 1e-4 )
			sum_not_one = FALSE;

		m++;
		//Check for User Interruption
		try{
		    Rcpp::checkUserInterrupt();
		}
		catch(Rcpp::internal::InterruptedException e){
		    //Print error and return
		    throw std::runtime_error("Execution stopped by the user");
		}
	}
	int Mstar_max = lqMpost.size();
	// Compute log expected value
	double lPoch_gamma_r = log_raising_factorial(r, gamma); // (gamma)_r
	
	// initialization
	std::vector<double> log_vect_res(Mstar_max, -inf );
	unsigned int idx_max{0};
	double val_max(log_vect_res[idx_max]);

	for(int Mstar=0; Mstar < Mstar_max; Mstar++){
		log_vect_res[Mstar] = std::log(Mstar) + lqMpost[Mstar] + 
							  lPoch_gamma_r - 
							  log_raising_factorial( r, (double)n + gamma*(double)(Mstar + Kn) );
		// Check if it is the new maximum
	    if(log_vect_res[Mstar]>val_max){
	    	idx_max = Mstar;
	       	val_max = log_vect_res[Mstar];
	       	//Rcpp::Rcout<<"idx_max: "<<idx_max<<"; val_max = "<<val_max<<std::endl;
	    }
	    //Rcpp::Rcout<<"log_vect_res["<<Mstar<<"]: "<<log_vect_res[Mstar]<<std::endl;
	}

	// Formula to compute the log of all the sums in a stable way
	return log_stable_sum(log_vect_res, TRUE, val_max, idx_max);	
}

// [[Rcpp::export]]
double compute_log_UBMarkov_FD( const int& Rmax, const double& gamma, const double& Lambda, 
								const int& Kn, const int& n,
						        const double& alpha_lev, unsigned int M_max)
{
	double inf = std::numeric_limits<double>::infinity();

	if(Rmax < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_FD: Rmax must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_FD: n must be at least one");
	if(Kn > n)
		throw std::runtime_error("Error in compute_log_UBMarkov_FD: Kn must be smaller or equal to n");
	if(gamma < 0)
		throw std::runtime_error("Error in compute_log_UBMarkov_FD: alpha must be positive");
	if(Lambda < 0)
		throw std::runtime_error("Error in compute_log_UBMarkov_FD: theta must be >= alpha");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in compute_log_UBMarkov_FD: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");

	double log_res{inf};
	for (int r = 1; r <= Rmax; r++)
	{
		// 1/r * ( log(ExpVal) - log(alpha_lev) )
		double lExp = log_ExpMr_FD(r, gamma, Lambda, Kn, n, M_max);
		//Rcpp::Rcout<<"lExp:"<<std::endl<<lExp<<std::endl;
		double log_res_r{ lExp - std::log(alpha_lev) };
		//Rcpp::Rcout<<"log_res_r:"<<std::endl<<log_res_r<<std::endl;
		log_res_r *= 1.0/double(r);
		//Rcpp::Rcout<<"log_res_r:"<<std::endl<<log_res_r<<std::endl;

		// Find minumum
		if(log_res_r < log_res){
			//Rcpp::Rcout<<"Min r: "<<r<<std::endl;
			log_res = log_res_r;
		}
	}

	return log_res;
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
	double Warg = (double)n/(alpha_lev - beta) * 
					(  std::sqrt( (std::log(1.0/beta))/(2.0*(double)n) ) + std::sqrt( (std::log(1.0/beta))/(2.0*(double)n) + Shat ) );
	double temp{gsl_sf_lambert_W0(Warg)};
	
	if(temp < 0)
		throw std::runtime_error("Error in compute_log_UBFreq_BeBe: the argument of the Lambert W function can not be negative");

	res *= temp;
	return res;				
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Features - Poisson process
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double log_ExpMr_BeBePois(const int& r, const double& alpha, const double& c, const double& gamma, const int& n)
{
	if(r < 1)
		throw std::runtime_error("Error in log_ExpMr_BeBePois: r must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in log_ExpMr_BeBePois: n must be at least one");
	if(alpha < 0.0)
		throw std::runtime_error("Error in log_ExpMr_BeBePois: alpha must be positive");
	if(alpha >= 1.0)
		throw std::runtime_error("Error in log_ExpMr_BeBePois: alpha must be strictly smaller than 1");
	if(c <= -alpha)
		throw std::runtime_error("Error in log_ExpMr_BeBePois: c must be > alpha");
	return std::log(gamma) + gsl_sf_lnbeta( (double)r - alpha, c + alpha + (double)n );
}

// [[Rcpp::export]]
double compute_log_UBMarkov_BeBePois( const int& Rmax, 
									  const double& alpha, const double& c, const double& gamma, 
									  const int& n, const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();

	if(Rmax < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBePois: Rmax must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBePois: n must be at least one");
	if(alpha < 0.0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBePois: alpha must be positive");
	if(alpha >= 1.0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBePois: alpha must be strictly smaller than 1");
	if(c <= -alpha)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBePois: c must be > alpha");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBePois: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");

	double log_res{inf};
	for (int r = 1; r <= Rmax; r++){
		// 1/r * ( log(ExpVal) - log(alpha_lev) )
		double log_res_r{ 0.0 };
		log_res_r = 1.0/double(r) * ( log_ExpMr_BeBePois(r, alpha, c, gamma, n) - std::log(alpha_lev) );
		// Find minumum
		if(log_res_r < log_res)
			log_res = log_res_r;
	}

	return log_res;
}

// [[Rcpp::export]]
double log_phi_n(const double& alpha, const double& c, const double& gamma, const int& n)
{
	if(n < 1)
		throw std::runtime_error("Error in log_phi_n: n must be at least one");
	if(gamma < 0.0)
		throw std::runtime_error("Error in log_phi_n: gamma must be positive");
	if(alpha < 0.0)
		throw std::runtime_error("Error in log_phi_n: alpha must be positive");
	if(alpha >= 1.0)
		throw std::runtime_error("Error in log_phi_n: alpha must be strictly smaller than 1");
	if(c <= -alpha)
		throw std::runtime_error("Error in log_phi_n: c must be > alpha");

	if(c > 1e-5 & alpha > 1e-5){
		// Use the direct formula
		double res{ std::log(gamma) + gsl_sf_lngamma(1.0 - alpha) - std::log(alpha) };
		res += -gsl_sf_lngamma(c+1.0) - gsl_sf_lngamma((double)n + c + 1.0);
		Rcpp::NumericVector lres_vec(2);
		lres_vec[0] = std::log((double)n + c) + gsl_sf_lngamma(c+1.0) + gsl_sf_lngamma((double)n + c + alpha);
		lres_vec[1] = std::log(c) + gsl_sf_lngamma(c+alpha) + gsl_sf_lngamma((double)n + c + 1.0);
		res += lres_vec[0] + std::log( 1.0 - std::exp( lres_vec[1] - lres_vec[0] ) );
		return res;
	}
	else{
		// Compute the whole sum
		Rcpp::NumericVector lres_vec(n);
		for(int i = 0; i < n; i++){
			lres_vec[i] = std::log(gamma) + gsl_sf_lnbeta(1.0 - alpha, c + alpha + i);
		}
		return log_stable_sum(lres_vec, TRUE);
	}

}

// [[Rcpp::export]]
double log_efpfBeBePois( const int& n, const int& Kn, const std::vector<int>& n_j, 
					     const double& alpha, const double& c, const double& gamma )
{
	double inf = std::numeric_limits<double>::infinity();

	if(n < 1)
		throw std::runtime_error("Error in log_efpfBeBePois: n must be at least one");
	if(Kn <= 0)
		throw std::runtime_error("Error in log_efpfBeBePois: Kn must be at least one");
	if(n_j.size() != Kn)
		throw std::runtime_error("Error in log_efpfBeBePois: the length of n_j must match the number of features Kn");
	std::for_each(n_j.cbegin(), n_j.cend(), [n](const int& x){ 
		if(x<=0)
			throw std::runtime_error("Error in log_efpfBeBePois: invalid number of features (n_j[j] <= 0)");
		if(x > n )  
			throw std::runtime_error("Error in log_efpfBeBePois: invalid number of features (n_j[j] > n)");
	});

	// Check domain
	if(c <= -alpha){
		Rcpp::Rcout<<"Caso proibito (c < -alpha)"<<std::endl;
		return -std::exp(20);
	}
	if(alpha < 0 || alpha > 1 - 1e-16){
		Rcpp::Rcout<<"Caso proibito (alpha < 0 or alpha > 1)"<<std::endl;
		return -std::exp(20);
	}
	if(gamma <= 0){
		Rcpp::Rcout<<"Caso proibito (gamma <= 0)"<<std::endl;
		return -std::exp(20);
	}

	double res{ 0.0 }; 
	res = -std::exp( log_phi_n(alpha,c,gamma,n) ) + Kn*std::log(gamma);
	for(std::size_t j = 0; j < n_j.size(); j++){
		res += gsl_sf_lnbeta( (double)n_j[j] - alpha, c + alpha + (double)(n - n_j[j]) );
	}
	if( res == inf || std::isnan(res) || res == -inf){
		Rcpp::Rcout<<"Error in log_efpfBeBePois: NaN, Inf or -Inf returned "<<std::endl;
		Rcpp::Rcout<<"res = "<<res<<std::endl;
		Rcpp::Rcout<<"alpha = "<<alpha<<std::endl;
		Rcpp::Rcout<<"c = "<<c<<std::endl;
		Rcpp::Rcout<<"n = "<<n<<std::endl;
		Rcpp::Rcout<<"Kn = "<<Kn<<std::endl;
		Rcpp::Rcout<<"Stampo n_j: ";		
		for(auto __v : n_j)
			Rcpp::Rcout<<__v<<", ";
		Rcpp::Rcout<<std::endl;
		throw std::runtime_error("Error: no good!");
	}
	return res;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Features - Mixed Poisson process
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double log_ExpMr_BeBeMixPois( const int& r, const double& alpha, const double& c, 
							  const int& n, const int& Kn,
							  const double& u, const double& v )
{
	if(r < 1)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: r must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: n must be at least one");
	if(Kn <= 0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: Kn must be at least one");
	if(alpha < 0.0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: alpha must be positive");
	if(alpha >= 1.0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: alpha must be strictly smaller than 1");
	if(c <= -alpha)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: c must be > alpha");
	if(u < 0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: u must be positive");
	if(v < 0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixPois: v must be positive");

	double res{ gsl_sf_lnbeta( (double)r - alpha, c + alpha + (double)n ) };
	res += std::log(u + (double)Kn) - std::log( v + std::exp(log_phi_n(alpha,c,1.0,n)) );
	return res;
}

// [[Rcpp::export]]
double compute_log_UBMarkov_BeBeMixPois( const int& Rmax, const double& alpha, const double& c,
									     const int& n, const int& Kn,
							             const double& u, const double& v, const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();

	if(Rmax < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: Rmax must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: n must be at least one");
	if(Kn <= 0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: Kn must be at least one");
	if(alpha < 0.0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: alpha must be positive");
	if(alpha >= 1.0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: alpha must be strictly smaller than 1");
	if(c <= -alpha)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: c must be > alpha");
	if(u < 0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: u must be positive");
	if(v < 0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: v must be positive");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");

	double log_res{inf};
	for (int r = 1; r <= Rmax; r++)
	{
		// 1/r * ( log(ExpVal) - log(alpha_lev) )
		double log_res_r{ 0.0 };
		log_res_r = 1.0/double(r) * ( log_ExpMr_BeBeMixPois( r, alpha, c, n, Kn, u, v ) - std::log(alpha_lev) );
		// Find minumum
		if(log_res_r < log_res)
			log_res = log_res_r;
	}

	return log_res;
}

// [[Rcpp::export]]
double log_efpfBeBeMixPois( const int& n, const int& Kn, const std::vector<int>& n_j, 
					        const double& alpha, const double& c,
					        const double& mu_gamma, const double& var_gamma )//const double& u, const double& v )
{
	double inf = std::numeric_limits<double>::infinity();

	if(n < 1)
		throw std::runtime_error("Error in log_efpfBeBeMixPois: n must be at least one");
	if(Kn <= 0)
		throw std::runtime_error("Error in log_efpfBeBeMixPois: Kn must be at least one");
	if(n_j.size() != Kn)
		throw std::runtime_error("Error in log_efpfBeBeMixPois: the length of n_j must match the number of features Kn");
	std::for_each(n_j.cbegin(), n_j.cend(), [n](const int& x){ 
		if(x<=0)
			throw std::runtime_error("Error in log_efpfBeBeMixPois: invalid number of features (n_j[j] <= 0)");
		if(x > n )  
			throw std::runtime_error("Error in log_efpfBeBeMixPois: invalid number of features (n_j[j] > n)");
	});

	// Check domain
	if(c <= -alpha){
		//Rcpp::Rcout<<"Caso proibito (c < -alpha)"<<std::endl;
		return -std::exp(20);
	}
	if(alpha < 0 || alpha > 1 - 1e-16){
		//Rcpp::Rcout<<"Caso proibito (alpha < 0 or alpha > 1)"<<std::endl;
		return -std::exp(20);
	}
	if(mu_gamma <= 0 || var_gamma <= 0){
		//Rcpp::Rcout<<"Caso proibito (mu_gamma or var_gamma < 1e-16)"<<std::endl;
		return -std::exp(20);
	}

	// Gamma reparametrization
	double u{mu_gamma*mu_gamma/var_gamma}; // shape hyperparameter
	double v{mu_gamma/var_gamma};          // rate hyperparameter

	double res{ 0.0 }; 
	double phi_n = std::exp( log_phi_n(alpha,c,1.0,n) );
	//Rcpp::Rcout<<"---------------"<<std::endl;
	//Rcpp::Rcout<<"phi_n = "<<phi_n<<std::endl;
	res = log_raising_factorial(Kn,u) + u * ( std::log(v) - std::log(v+phi_n) ) - (double)Kn * std::log(v+phi_n);
	for(std::size_t j = 0; j < n_j.size(); j++){
		res += gsl_sf_lnbeta( (double)n_j[j] - alpha, c + alpha + (double)(n - n_j[j]) );
	}
	if( res == inf || std::isnan(res) || res == -inf){
		Rcpp::Rcout<<"Error in log_efpfBeBeMixPois: NaN, Inf or -Inf returned "<<std::endl;
		Rcpp::Rcout<<"res = "<<res<<std::endl;
		Rcpp::Rcout<<"alpha = "<<alpha<<std::endl;
		Rcpp::Rcout<<"c = "<<c<<std::endl;
		Rcpp::Rcout<<"mu_gamma = "<<mu_gamma<<std::endl;
		Rcpp::Rcout<<"var_gamma = "<<var_gamma<<std::endl;
		Rcpp::Rcout<<"u = "<<u<<std::endl;
		Rcpp::Rcout<<"v = "<<v<<std::endl;
		Rcpp::Rcout<<"n = "<<n<<std::endl;
		Rcpp::Rcout<<"Kn = "<<Kn<<std::endl;
		Rcpp::Rcout<<"Stampo n_j: ";		
		for(auto __v : n_j)
			Rcpp::Rcout<<__v<<", ";
		Rcpp::Rcout<<std::endl;
		throw std::runtime_error("Error in log_efpfBeBeMixPois");
	}
	return res;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Features - Mixed Binomial process (Negative Binomial prior)
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double log_ExpMr_BeBeMixNBin( const int& r, const double& a, const double& b, 
							  const int& n, const int& Kn,
							  const double& r_nb, const double& p_nb )
{
	if(r < 1)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixNBin: r must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixNBin: n must be at least one");
	if(Kn <= 0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixNBin: Kn must be at least one");
	if(a < 0.0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixNBin: a must be positive");
	if(b < 0.0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixNBin: c must be positive");
	if(r_nb < 0)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixNBin: r_nb must be positive");
	if(p_nb < 0 || p_nb > 1)
		throw std::runtime_error("Error in log_ExpMr_BeBeMixNBin: p_nb must be in (0,1)");

	double kappa_n{std::exp( log_raising_factorial(n,b) - log_raising_factorial(n, a+b) ) };
	double r_nb_prime = r_nb + Kn;
	double p_nb_prime = 1.0 - kappa_n*(1.0 - p_nb);
	double l_postmean_M{ std::log(r_nb_prime) + std::log(1.0 - p_nb_prime) - std::log(p_nb_prime)  };
	double res = l_postmean_M + log_raising_factorial(r,a) + log_raising_factorial(n,b) - log_raising_factorial(r+n, a+b) ;
	return res;
}

// [[Rcpp::export]]
double compute_log_UBMarkov_BeBeMixNBin( const int& Rmax, const double& a, const double& b,
									     const int& n, const int& Kn,
							             const double& r_nb, const double& p_nb, 
							             const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();

	if(Rmax < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: Rmax must be at least one");
	if(n < 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: n must be at least one");
	if(Kn <= 0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: Kn must be at least one");
	if(a < 0.0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: a must be positive");
	if(b < 0.0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: b must be positive");
	if(r_nb < 0)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: r_nb must be positive");
	if(p_nb < 0 || p_nb > 1)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: p_nb must be in (0,1)");
	if(alpha_lev < 1e-16 || alpha_lev > 1 - 1e-16)
		throw std::runtime_error("Error in compute_log_UBMarkov_BeBeMixPois: alpha_lev must be at scalar in the range (1e-16,1-1e-16)");

	double log_res{inf};
	for (int r = 1; r <= Rmax; r++)
	{
		// 1/r * ( log(ExpVal) - log(alpha_lev) )
		double log_res_r{ 0.0 };
		log_res_r = 1.0/double(r) * ( log_ExpMr_BeBeMixNBin( r, a, b, n, Kn, r_nb,  p_nb ) - std::log(alpha_lev) );
		// Find minumum
		if(log_res_r < log_res)
			log_res = log_res_r;
	}

	return log_res;
}

// [[Rcpp::export]]
double log_efpfBeBeMixNBin( const int& n, const int& Kn, const std::vector<int>& n_j, 
					        const double& a, const double& b,
					        const double& mu_nb, const double& var_nb )
{
	double inf = std::numeric_limits<double>::infinity();

	if(n < 1)
		throw std::runtime_error("Error in log_efpfBeBeMixPois: n must be at least one");
	if(Kn <= 0)
		throw std::runtime_error("Error in log_efpfBeBeMixPois: Kn must be at least one");
	if(n_j.size() != Kn)
		throw std::runtime_error("Error in log_efpfBeBeMixPois: the length of n_j must match the number of features Kn");
	std::for_each(n_j.cbegin(), n_j.cend(), [n](const int& x){ 
		if(x<=0)
			throw std::runtime_error("Error in log_efpfBeBeMixPois: invalid number of features (n_j[j] <= 0)");
		if(x > n )  
			throw std::runtime_error("Error in log_efpfBeBeMixPois: invalid number of features (n_j[j] > n)");
	});

	// Check domain
	if(a <= 0 || b <= 0){
		return -std::exp(20);
	}
	if(mu_nb <= 0 || var_nb <= 0){
		//Rcpp::Rcout<<"Caso proibito (mu_nb or var_nb < 1e-16)"<<std::endl;
		return -std::exp(20);
	}
	if(mu_nb >= var_nb){
		Rcpp::Rcout<<"Caso proibito (mu_nb >= var_nb)"<<std::endl;
		return -std::exp(20);
	}
	// Negative binomial reparametrization
	double p_nb{mu_nb/var_nb};               // prob hyperparameter
	double r_nb{mu_nb*mu_nb/(var_nb-mu_nb)}; // size hyperparameter

	double res{ gsl_sf_lngamma((double)Kn + r_nb) - gsl_sf_lngamma((double)Kn - 1.0) - gsl_sf_lngamma(r_nb) };
	double temp{ std::exp( log_raising_factorial(n,b) - log_raising_factorial(n,a+b) ) };
	res += r_nb*std::log(p_nb) + (double)Kn*std::log(1.0 - p_nb) - (r_nb+(double)Kn) * std::log( 1.0 - (1.0-p_nb)*temp );
	res -= (double)Kn*log_raising_factorial(n,a+b);
	for(std::size_t j = 0; j < n_j.size(); j++){
		res += log_raising_factorial(n_j[j], a) + log_raising_factorial((double)n - n_j[j], b);
	}
	if( res == inf || std::isnan(res) || res == -inf){
		Rcpp::Rcout<<"Error in log_efpfBeBeMixPois: NaN, Inf or -Inf returned "<<std::endl;
		Rcpp::Rcout<<"res = "<<res<<std::endl;
		Rcpp::Rcout<<"a = "<<a<<std::endl;
		Rcpp::Rcout<<"b = "<<b<<std::endl;
		Rcpp::Rcout<<"mu_nb = "<<mu_nb<<std::endl;
		Rcpp::Rcout<<"var_nb = "<<var_nb<<std::endl;
		Rcpp::Rcout<<"r_nb = "<<r_nb<<std::endl;
		Rcpp::Rcout<<"p_nb = "<<p_nb<<std::endl;
		Rcpp::Rcout<<"n = "<<n<<std::endl;
		Rcpp::Rcout<<"Kn = "<<Kn<<std::endl;
		Rcpp::Rcout<<"Stampo n_j: ";		
		for(auto __v : n_j)
			Rcpp::Rcout<<__v<<", ";
		Rcpp::Rcout<<std::endl;
		throw std::runtime_error("Error in log_efpfBeBeMixPois");
	}
	return res;
}
// --------------------------------------------------------------------------------------------
// Test functions
// --------------------------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector prova(Rcpp::NumericVector x)
{
  return x+x;
}

