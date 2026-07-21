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

double log_stable_sum(const Rcpp::NumericVector& a, const bool is_log, const double& val_max)
{

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
double log_stable_sum(const Rcpp::NumericVector& a, const bool is_log)
{
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

// [[Rcpp::export]]
double LambertW(const double& x){
	return gsl_sf_lambert_W0(x);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Sample while controlling S
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector sim_CntS_Zipfs(double S, double s) {
  if(S <= 0.0)
    throw std::runtime_error("Error in sim_CntS_Constant. S must be positive");
  if( s > 1.0){
  	double Smax = RiemannZeta(s);
  	if(S > Smax)
  		throw std::runtime_error("Error in sim_CntS_Constant. It is impossible to achieve such S given the requested s");
  }

  const int Mmax = 1000000;
  double cumS = 0.0;
  int j = 1;
  std::vector<double> w;
  w.reserve(1000); // preallocate modestly to avoid frequent reallocations
  
  bool flag = true;
  
  while(flag) {
    double newval = std::pow(j + 1.0, -s);
    w.push_back(newval);
    cumS += newval;
    
    if(cumS >= S)
      flag = false;
    
    if((j + 1) > Mmax)
      flag = false;
    
    j++;
  }
  
  double total = std::accumulate(w.begin(), w.end(), 0.0);
  if (total < S){
  	Rcpp::Rcout<<"total = "<<total<<std::endl;
  	throw std::runtime_error("Error in sim_CntS_Constant. The sum is smaller than S after 1 million points");
  }
    
  
  return wrap(w);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Sample from truncated Beta process
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix r_BetaPr( const int& Nrep, const int& Natoms, const double& gamma, const int& seed)
{
	sample::rbeta beta; 
	sample::GSL_RNG engine(seed);

	if(Nrep <= 0)
		throw std::runtime_error("Error in r_BetaPr: Nrep must be positive");
	if(Natoms < 1)
		throw std::runtime_error("Error in r_BetaPr: Natoms must be at least one");
	if(gamma <= 0)
		throw std::runtime_error("Error in r_BetaPr: gamma must be > 0");

	NumericMatrix res(Nrep,Natoms);
	for(int b = 0; b < Nrep; b++){
		double cumprod{1.0};
		for(int i = 0; i < Natoms; i++){
			//Rcpp::Rcout<<"("<<b<<","<<i<<"); cumprod = "<<cumprod<<std::endl;
			double beta_i =  beta(engine, gamma, 1.0);
			res(b,i) = beta_i * cumprod ; // save i-th value
			cumprod *= beta_i; // update cumulative product 
		}
	}

	return res;
}

// Note: if sigma < 1e-5, then it is consider as 0
// [[Rcpp::export]]
NumericMatrix r_3ParamBetaPr( const int& Nrep, const int& Natoms, 
							  const double& gamma, const double& sigma, const double& c,
							  const int& seed)
{
	sample::rgamma rgamma; 
	sample::GSL_RNG engine(seed);

	if(Nrep <= 0)
		throw std::runtime_error("Error in r_3ParamBetaPr: Nrep must be positive");
	if(Natoms < 1)
		throw std::runtime_error("Error in r_3ParamBetaPr: Natoms must be at least one");
	if(gamma <= 0)
		throw std::runtime_error("Error in r_3ParamBetaPr: gamma must be > 0");
	if(sigma < 0 || sigma >= 1)
		throw std::runtime_error("Error in r_3ParamBetaPr: sigma must be in [0,1)");
	if(c <= -sigma)
		throw std::runtime_error("Error in r_3ParamBetaPr: c must be > -sigma");

	auto invPsi = [gamma, sigma, c](const double& xi){
		if(xi < 0)
			throw std::runtime_error("Error in invPsi: xi can not be negative");

		double arg = (xi)/(gamma*c);
		double res;
		if( sigma < 1e-5 ){
			res = gsl_expm1(arg);
		}
		else{
			arg *= sigma;
			res = gsl_expm1( gsl_log1p( arg ) ); 
		}
		if(res < 0)
			throw std::runtime_error("Error in invPsi: res can not be negative");
		return res;
	};

	NumericMatrix res(Nrep,Natoms);
	for(int b = 0; b < Nrep; b++){
		for(int i = 0; i < Natoms; i++){
			// double xi = rgamma(engine, 1.0, 1.0);
			double xi = rgamma(engine, i+1.0, 1.0);
			double ti = invPsi(xi);
			double z = rgamma(engine, 1.0 - sigma, 1.0/(1.0 + ti));
			double y = rgamma(engine, c+sigma, 1.0);
			double z_over_y = z/y;
			res(b,i) = (z_over_y)/(z_over_y + 1.0);
		}
	}

	return res;
}


// [[Rcpp::export]]
NumericMatrix r_iid3ParamBetaPr( const int& Nrep, const int& Natoms, 
								 const double& gamma, const double& sigma, const double& c,
							     const int& seed)
{
	sample::rgamma rgamma; 
	sample::runif  runif; 
	sample::GSL_RNG engine(seed);

	if(Nrep <= 0)
		throw std::runtime_error("Error in r_3ParamBetaPr: Nrep must be positive");
	if(Natoms < 1)
		throw std::runtime_error("Error in r_3ParamBetaPr: Natoms must be at least one");
	if(gamma <= 0)
		throw std::runtime_error("Error in r_3ParamBetaPr: gamma must be > 0");
	if(sigma < 0 || sigma >= 1)
		throw std::runtime_error("Error in r_3ParamBetaPr: sigma must be in [0,1)");
	if(c <= -sigma)
		throw std::runtime_error("Error in r_3ParamBetaPr: c must be > -sigma");

	auto invPsi = [gamma, sigma, c](const double& xi){
		if(xi < 0)
			throw std::runtime_error("Error in invPsi: xi can not be negative");

		double arg = (xi)/(gamma*c);
		double res;
		if( sigma < 1e-5 ){
			res = gsl_expm1(arg);
		}
		else{
			arg *= sigma;
			res = gsl_expm1( gsl_log1p( arg ) ); 
		}
		if(res < 0)
			throw std::runtime_error("Error in invPsi: res can not be negative");
		return res;
	};

	NumericMatrix res(Nrep,Natoms);
	double t = invPsi((double)Natoms);
	for(int b = 0; b < Nrep; b++){
		for(int i = 0; i < Natoms; i++){
			double G = rgamma(engine, 1.0 - sigma, 1.0);
			double u = runif(engine);
			double z;
			if(sigma < 1e-5){
				z = G * std::pow(t+1.0, -(1.0 - u));
			}
			else{
				double arg = gsl_expm1( sigma * std::log(t+1.0) );
				double loga = std::log1p((1.0 - u) * arg);
				z = G * std::exp(- (loga/sigma) );
			}
			double y = rgamma(engine, c+sigma, 1.0);
			double z_over_y = z/y;
			res(b,i) = z_over_y / (z_over_y + 1.0);
		}
	}

	return res;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Features - Frequentist - Bounded alphabet
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double compute_UB_bdd_oracle( const int& n, const Rcpp::NumericVector& pj, const double& b, const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();
	int M = pj.size();
	if(n <= 0)
		throw std::runtime_error("Error in compute_UB_bdd_oracle: n must be strictly positive ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_UB_bdd_oracle: alpha_lev must be in (0,1) ");
	if( b <= 0 )
		throw std::runtime_error("Error in compute_UB_bdd_oracle: b must be strictly positive");

	std::vector<double> log_mp_vec(M,-inf); // vector of log values
	double max{-inf};
	int idx_max{-1};
	for(int j=0; j < M; j++){
		log_mp_vec[j] = b * gsl_log1p( -pj[j] ); // compute log of jth element
		if(log_mp_vec[j] > max){
			max = log_mp_vec[j]; // save max
			idx_max = j; // save position max
		}	
		if(std::isnan(log_mp_vec[j])){
			throw std::runtime_error("Error in compute_UB_bdd_oracle. Get a NaN");
		}
	}
	//mp = std::exp( log_stable_sum(log_mp_vec, TRUE, max, idx_max) );
	double log_mp = log_stable_sum(log_mp_vec, TRUE, max, idx_max);
	return 1.0/( (double)n - b )*( log_mp - std::log(alpha_lev) );
}


// [[Rcpp::export]]
double compute_UB_analytical( const int& n, const Rcpp::IntegerVector& Nj, const int& M, 
							  const double& b, const double& alpha_lev, bool IsRegular)
{
	double inf = std::numeric_limits<double>::infinity();

	if(n <= 0)
		throw std::runtime_error("Error in compute_UB_analytical: n must be strictly positive ");
	if(M <= 0)
		throw std::runtime_error("Error in compute_UB_analytical: M must be strictly positive ");
	if(Nj.size()!= M)
		throw std::runtime_error("Error in compute_UB_analytical: the length of Nj must be equal to M ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_UB_analytical: alpha_lev must be in (0,1) ");
	if( b <= 0 )
		throw std::runtime_error("Error in compute_UB_analytical: b must be strictly positive");

	double a1 = 0.99*alpha_lev; // \alpha in the paper
	double a2 = 0.01*alpha_lev; // \delta in the paper

	std::vector<double> log_mp_vec(M+1,-inf); // vector of log values
	double max{-inf};
	int idx_max{-1};
	for(int j=0; j < M; j++){
		if(Nj[j] < n){
			log_mp_vec[j] = b * gsl_log1p( -(double)Nj[j]/(double)n ); // compute log of jth element	
		}
		if(log_mp_vec[j] > max){
			max = log_mp_vec[j]; // save max
			idx_max = j; // save position max
		}
		//if(std::isnan(log_mp_vec[j])){
			//throw std::runtime_error("Error in compute_UB_analytical. Get a NaN");
		//}
	}
	//mp = std::exp( log_stable_sum(log_mp_vec, TRUE, max, idx_max) );
	if(!IsRegular){
		log_mp_vec[M] = std::log(b) + 0.5*( std::log( (double)M ) - std::log( (double)n ) + std::log(-std::log(a2)) );
		if(log_mp_vec[M] > max){
			max = log_mp_vec[M]; // save max
			idx_max = M; // save position max
		}
	}
	double log_mp = log_stable_sum(log_mp_vec, TRUE, max, idx_max);
	return 1.0/( (double)n - b )*( log_mp - std::log(a1) );
}

// [[Rcpp::export]]
double compute_lmbar_bounded( const int& n, const Rcpp::IntegerVector& Nj, const int& M, 
							               const double& b, const double& alpha_lev, const int& nmax)
{
	double inf = std::numeric_limits<double>::infinity();

	if(n <= 0)
		throw std::runtime_error("Error in compute_mbar_bounded: n must be strictly positive ");
	if(M <= 0)
		throw std::runtime_error("Error in compute_mbar_bounded: M must be strictly positive ");
	if(Nj.size()!= M)
		throw std::runtime_error("Error in compute_mbar_bounded: the length of Nj must be equal to M ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_mbar_bounded: alpha_lev must be in (0,1) ");
	if( b <= 0 )
		throw std::runtime_error("Error in compute_mbar_bounded: b must be strictly positive");

	double a1 = 0.99*alpha_lev; // \alpha in the paper
	double a2 = 0.01*alpha_lev; // \delta in the paper

	std::vector<double> log_mp_vec(M+1,-inf); // vector of log values
	double max{-inf};
	int idx_max{-1};
	for(int j=0; j < M; j++){
		if(Nj[j] < n){
			log_mp_vec[j] = b * gsl_log1p( -(double)Nj[j]/(double)n ); // compute log of jth element	
		}
		if(log_mp_vec[j] > max){
			max = log_mp_vec[j]; // save max
			idx_max = j; // save position max
		}
	}
	
	// define extra term that is the eps(..) term
	log_mp_vec[M] = std::log(b) + 0.5*( std::log( (double)M ) - std::log( (double)n ) + std::log(-std::log(a2/nmax)) );
	if(log_mp_vec[M] > max){
		max = log_mp_vec[M]; // save max
		idx_max = M; // save position max
	}
	double log_mp = log_stable_sum(log_mp_vec, TRUE, max, idx_max);
	return log_mp;
}


// [[Rcpp::export]]
double compute_LB_analytical( const int& n, const Rcpp::IntegerVector& Nj, const double& alpha_lev)
{
	double inf = std::numeric_limits<double>::infinity();

	if(n <= 0)
		throw std::runtime_error("Error in compute_UB_analytical: n must be strictly positive ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_UB_analytical: alpha_lev must be in (0,1) ");

	int n_zeros = sum(Nj == 0);
	if(n_zeros == 0)
		return 0.0;

	// else
	double lb = 1.0/(double)n * std::log((double)n_zeros/alpha_lev);
	return lb;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------
//	Features - Frequentist - Unbounded alphabet
//------------------------------------------------------------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double compute_UB_Unbdd_oracle(const int& n, const double& alpha_lev, const double& S)
{
	if(n <= 0)
		throw std::runtime_error("Error in compute_UB_Unbdd_oracle: n must be strictly positive ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_UB_Unbdd_oracle: alpha_lev must be in (0,1) ");
	if(S < 0)
		throw std::runtime_error("Error in compute_UB_Unbdd_oracle: S must be positive ");

	double Warg = ( S*(double)n )/( -std::log(1.0 - alpha_lev) ); 

	if(Warg < 0)
		throw std::runtime_error("Error in compute_UB_Unbdd_oracle: the argument of the Lambert W function can not be negative");

	// UB = 1/n * W(..)
	double res = gsl_sf_lambert_W0(Warg)/(double)n;
	
	return std::min(res,1.0);				
}

// [[Rcpp::export]]
double compute_UB_intersection(const int& n, const double& alpha_lev, const double& beta_lev, const double& Shat)
{
	if(n <= 0)
		throw std::runtime_error("Error in compute_UB_intersection: n must be strictly positive ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_UB_intersection: alpha_lev must be in (0,1) ");
	if( beta_lev <= 0 )
		throw std::runtime_error("Error in compute_UB_intersection: beta_lev must be strictly positive");
	if( (1.0 - alpha_lev + beta_lev) <= 0 )
		throw std::runtime_error("Error in compute_UB_intersection: 1.0 - alpha_lev + beta_lev must be strictly positive");
	if(Shat < 0)
		throw std::runtime_error("Error in compute_UB_intersection: Shat must be positive ");

	double Warg = (double)n/(-std::log(1.0 - alpha_lev + beta_lev)); 
	Warg *= (  std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) ) + std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) + Shat ) ) *
			(  std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) ) + std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) + Shat ) ) ;

	if(Warg < 0)
		throw std::runtime_error("Error in compute_UB_intersection: the argument of the Lambert W function can not be negative");

	// UB = 1/n * W(..)
	double res = gsl_sf_lambert_W0(Warg)/(double)n;
	
	return std::min(res,1.0);				
}

// [[Rcpp::export]]
double compute_UB_rnorm(const int& n, const double& alpha_lev, const double& beta_lev, const int& r, const double& Shat)
{
	if(n <= 0)
		throw std::runtime_error("Error in compute_UB_rnorm: n must be strictly positive ");
	if(alpha_lev <= 0 || alpha_lev >= 1)
		throw std::runtime_error("Error in compute_UB_rnorm: alpha_lev must be in (0,1) ");
	if( beta_lev <= 0 )
		throw std::runtime_error("Error in compute_UB_rnorm: beta_lev must be strictly positive");
	if( (alpha_lev - beta_lev) <= 0 )
		throw std::runtime_error("Error in compute_UB_rnorm:  alpha_lev + beta_lev must be strictly positive");
	if(Shat < 0)
		throw std::runtime_error("Error in compute_UB_rnorm: Shat must be positive ");
	if(r <= 0)
		throw std::runtime_error("Error in compute_UB_rnorm: r must be strictly positive ");
	
	double Sstar = (  std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) ) + std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) + Shat ) ) *
				   (  std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) ) + std::sqrt( (std::log(1.0/beta_lev))/(2.0*(double)n) + Shat ) ) ;
	
	double res{ 0 };
	res += 1.0/(double)r * ( std::log(Sstar) - std::log(alpha_lev - beta_lev) ) ;
	res += (double)(r-1)/(double)r * ( std::log( (double)(r-1) ) - std::log( (double)(n+r-1) ) );
	res += (double)(n)/(double)r * ( std::log( (double)(n) ) - std::log( (double)(n+r-1) ) );
	res = std::exp(res);

	return std::min(res,1.0);			
}
