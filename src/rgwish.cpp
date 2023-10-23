// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2012 - 2020  Reza Mohammadi                                                   |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  
#include "rgwish.h"
#include "copula.h"
#include "matrix.h"
#include "util.h"
#include <Rcpp.h>
  
  
extern "C" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sampling from Multivariate Normal Distribution
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rmvn_c( double rand_values[], double mus[], double K[], int p )
{
  int dim = p, pxp = dim * dim, one = 1, num_upper = dim*( dim + 1 )/2;
  double alpha = 1.0, beta   = 0.0;
  char upper = 'U';																	
  
  vector<double> rand_normal( dim, 0.0 );
  vector<double> copy_K( pxp, 0.0 );
  vector<double> cov( pxp, 0.0 );
  vector<double> chol_cov( pxp, 0.0 );
  vector<double> chol_cov_compact( num_upper, 0.0 );
  
  memcpy( &copy_K[0], K, sizeof( double ) * pxp );
  inverse( &copy_K[0], &cov[0], &dim );	
  cholesky( &cov[0], &chol_cov[0], &dim );	
  
  // NOTE ON PARALLELIZATION: My testing showed for p~=300, parallelization takes
  // more time than just straight C.  I'm leaving the code for this commented in case
  // we try much larger values for p.
  
  // - - - Sample a vector of Standard Normal Values - - - - - - - - - - - - - - - - - - - - - - - |
  GetRNGstate();
  // #pragma omp parallel
  // {
    // #pragma omp for
      for( int i = 0; i < dim; i++ )
        rand_normal[ i ] = norm_rand();
  // }
  PutRNGstate();
  
  // - - - Get compact matrix for LAPACK code  - - - - - - - - - - - - - - - - - - - - - - - |
  // #pragma omp parallel
  // {
    // #pragma omp for collapse(2)
      for( int i = 0; i < dim; i++ )
        for( int j = 0; j <= i; j++ )
          chol_cov_compact[ j + i ] = chol_cov[ j * dim + i ];	
  // }
  // - - - Calculate the Random Normal Values - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  // rand_values = chol_cov %*% rand_normal
  F77_NAME(dspmv)( &upper, &dim, &alpha, &chol_cov_compact[0], &rand_normal[0], &one, &beta, rand_values, &one FCONE );
  
  // #pragma omp parallel
  // {
    // #pragma omp for
      for( int i = 0; i < dim; i++ )
        rand_values[ i ] += mus[ i ];	
  // }

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sampling from Wishart distribution, in which Ts = chol( solve( Ds ) )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rwish_c( double Ts[], double K[], int *b, int *p )
{
	int dim = *p, pxp = dim * dim, bK = *b;
	double alpha = 1.0, beta   = 0.0;
	char transT  = 'T', transN = 'N', side = 'R', upper = 'U';																	

	vector<double> psi( pxp, 0.0 ); 

	// - - - Sample values in Psi matrix - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	GetRNGstate();
	//#pragma omp parallel for
	for( int i = 0; i < dim; i++ )
		psi[ i * dim + i ] = sqrt( Rf_rgamma( ( bK + dim - i - 1 ) / 2.0, 2.0 ) );
		//psi[i * dim + i] = sqrt( rchisq( bK + dim - i - 1 ) );

	//#pragma omp parallel for
	for( int j = 1; j < dim; j++ )
		for( int i = 0; i < j; i++ )
			psi[ j * dim + i ] = norm_rand();
	PutRNGstate();
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

  // C = psi %*% Ts   I used   psi = psi %*% Ts
	// dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
	F77_NAME(dtrmm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, Ts, &dim, &psi[0], &dim FCONE FCONE FCONE FCONE );

	// This is where K is updated.
	// K = t(C) %*% C 
	// LAPACK function to compute  C := alpha * A * B + beta * C																				
	F77_NAME(dgemm)( &transT, &transN, &dim, &dim, &dim, &alpha, &psi[0], &dim, &psi[0], &dim, &beta, K, &dim FCONE FCONE );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// G is adjacency matrix which has zero in its diagonal // threshold = 1e-8
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

void rgwish_c( double G[], double Ts[], double K[], int *b, int *p, double *threshold, int *failed )
{
  int info, i, j, l, size_node, one = 1, dim = *p, pxp = dim * dim;	
  
  double threshold_c = *threshold;
  double alpha = 1.0, beta  = 0.0;
  int num_iters = 10000;

  char transN  = 'N', uplo  = 'U'; 
  
  rwish_c( Ts, K, b, &dim );
  
  vector<double> sigma_start( pxp ); 
  inverse( K, &sigma_start[0], &dim );
	// Now that I have a randomly selected Wishart, I edit it iteratively so that
	// it ''matches'' my G.
	
	vector<double> sigma( sigma_start ); 
	vector<double> sigma_last( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_start_i( dim ); 

	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used

  double mean_diff = 1.0;
    
  int counter = 0;
	while( (mean_diff > threshold_c) and ( counter < num_iters ) )
	{
	  counter++;
		memcpy( &sigma_last[0], &sigma[0], sizeof( double ) * pxp );
		
		for( i = 0; i < dim; i++ )
		{
			// Count size of node
			size_node = 0;
			for( j = 0; j < dim; j++ ) size_node += (int)G[ j * dim + i ];

			if( size_node > 0 )
			{
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					if( (int)G[ j * dim + i ] )
					{
						sigma_start_N_i[ l ] = sigma_start[ i * dim + j ]; 
						N_i[ l++ ]           = j;
					}else
						beta_star[ j ] = 0.0; 
				}
				// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
				
				sub_matrix( &sigma[0], &sigma_N_i[0], &N_i[0], &size_node, &dim );
					
				// A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
				F77_NAME(dposv)( &uplo, &size_node, &one, &sigma_N_i[0], &size_node, &sigma_start_N_i[0], &size_node, &info FCONE );

				for( j = 0; j < size_node; j++ ) beta_star[ N_i[ j ] ] = sigma_start_N_i[ j ];
				
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &sigma[0], &dim, &beta_star[0], &dim, &beta, &sigma_start_i[0], &dim FCONE FCONE );
				
				for( j = 0; j < i; j++ )
				{
					sigma[ j * dim + i ] = sigma_start_i[ j ];
					sigma[ i * dim + j ] = sigma_start_i[ j ];
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					sigma[ j * dim + i ] = sigma_start_i[ j ];
					sigma[ i * dim + j ] = sigma_start_i[ j ];
				}
			} else {
				for( j = 0; j < i; j++ )
				{
					sigma[ j * dim + i ] = 0.0;
					sigma[ i * dim + j ] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					sigma[ j * dim + i ] = 0.0;
					sigma[ i * dim + j ] = 0.0;
				}
			} 
		}

		mean_diff = fabs( static_cast<double>( sigma[ 0 ] - sigma_last[ 0 ] ) );
		for( i = 1; i < pxp; i++ )
		    mean_diff += fabs( static_cast<double>( sigma[ i ] - sigma_last[ i ] ) );
		mean_diff /= pxp;
	}
	
	if(counter==num_iters){
	  *failed = 1;
	}
	
	
	 inverse( &sigma[0], K, &dim );
}
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Part of function "gnorm"
// which is for calculating Normalizing constant of G-Wishart distribution 
// based on Monto Carlo algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_exp_mc( int G[], int nu[], int *b, double H[], int *check_H, int *mc, int *p, double f_T[] )
{
	int iter, i, j, ij, h, r, mc_iter = *mc, dim = *p, pxp = dim * dim, b_c = *b;
	
	double sumPsi, sumPsiH, sumPsiHi, sumPsiHj;
	double max_numeric_limits_ld = numeric_limits<double>::max() / 1000;
	double min_numeric_limits_ld = numeric_limits<double>::min() * 1000;
	
	vector<double> psi( pxp, 0.0 );      

	GetRNGstate();
	if( *check_H == 1 )
	{ 
		for( iter = 0; iter < mc_iter; iter++ ) 
		{
			for( i = 0; i < dim; i++ )
				psi[ i * dim + i ] = sqrt( Rf_rgamma( ( b_c + nu[ i ] ) / 2.0, 2.0 ) );
				//psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );

			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					//if( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 ); else psi[ij] = 0.0;
					psi[ ij ] = ( G[ ij ] == 1 ) ? norm_rand() : 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if( G[ ij ] == 0 )
					{
						psi[ ij ] = 0.0; // it's not necessary
						if( i > 0 )  
						{
							sumPsi = 0.0;
							//sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] )
							// for( h = 0; h < ( i - 1 ); h++ )
							for( h = 0; h < i; h++ )
							{
								if( sumPsi == R_PosInf ) sumPsi = max_numeric_limits_ld;	
								if( sumPsi == R_NegInf ) sumPsi = min_numeric_limits_ld;	
								sumPsi += ( psi[ i * dim + h ] * psi[ j * dim + h ] );
							}
							
							//psi[i, j] <- - sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] ) / psi[i, i]
							psi[ ij ] = - sumPsi / psi[ i * dim + i ];
						}
						
						if( psi[ ij ] == R_PosInf ) psi[ ij ] = max_numeric_limits_ld;	
						if( psi[ ij ] == R_NegInf ) psi[ ij ] = min_numeric_limits_ld;	
												
						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[ iter ] += ( psi[ ij ] * psi[ ij ] ); 
					}
				}
		
			// checking Inf values
			if( f_T[ iter ] == R_PosInf ) f_T[ iter ] = max_numeric_limits_ld;			
		} 
	}else{
		for( iter = 0; iter < mc_iter; iter++ ) 
		{
			for( i = 0; i < dim; i++ )
				psi[ i * dim + i ] = sqrt( Rf_rgamma( ( b_c + nu[ i ] ) / 2.0, 2.0 ) );
				//psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );

			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					//if( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 ); elsepsi[ij] = 0.0;
					psi[ ij ] = ( G[ ij ] == 1 ) ? norm_rand() : 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if( G[ ij ] == 0 )
					{
						//psi[i, j] = - sum( psi[ i, i : ( j - 1 ) ] * H[ i : ( j - 1 ), j ] )
						sumPsiH = 0.0;
						for( h = i; h < j; h++ )
						{
							if( sumPsiH == R_PosInf ) sumPsiH = max_numeric_limits_ld;	
							if( sumPsiH == R_NegInf ) sumPsiH = min_numeric_limits_ld;	
							sumPsiH += ( psi[ h * dim + i ] * H[ j * dim + h ] ); 
						}
						psi[ ij ] = - sumPsiH;
						
						if( i > 0 )  //if( i > 1 )
							for( r = 0; r < i; r++ ) //for( r in 1 : ( i - 1 ) )
							{
								//sum( psi[ r, r : i ] * H[ r : i, i ] )
								sumPsiHi = 0.0;
								for( h = r; h < i + 1; h++  )
								{
									if( sumPsiHi == R_PosInf ) sumPsiHi = max_numeric_limits_ld;	
									if( sumPsiHi == R_NegInf ) sumPsiHi = min_numeric_limits_ld;	
									sumPsiHi += ( psi[ h * dim + r ] * H[ i * dim + h ] );	
								}
									
								//sum( psi[ r, r : j ] * H[ r : j, j ] ) )
								sumPsiHj = 0.0;
								for( h = r; h < j + 1; h++  )
									sumPsiHj += ( psi[ h * dim + r ] * H[ j * dim + h ] );
								
								//psi[i, j] <- psi[i, j] - ( ( sum( psi[ r, r : i ] * H[ r : i, i ] ) ) * ( sum( psi[ r, r : j ] * H[ r : j, j ] ) ) ) / ( psi[i, i] )
								psi[ ij ] -= ( sumPsiHi * sumPsiHj ) / psi[ i * dim + i ];
							}

						if( psi[ ij ] == R_PosInf ) psi[ ij ] = max_numeric_limits_ld;	
						if( psi[ ij ] == R_NegInf ) psi[ ij ] = min_numeric_limits_ld;	

						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[ iter ] += ( psi[ ij ] * psi[ ij ] ); 
					}
				}
		
			// checking Inf values
			if( f_T[ iter ] == R_PosInf ) f_T[ iter ] = max_numeric_limits_ld;			
		}
	}
	PutRNGstate();	
} 
      
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
     

} // End of Exern C     


using namespace Rcpp;

//' Samples from a multivariate normal distribution.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
Rcpp::NumericVector rmvn_Rcpp( Rcpp::NumericVector mus, 
                            Rcpp::NumericVector K, int p) {
  
  Rcpp::NumericVector rand_values(p);
  
  rmvn_c( REAL(rand_values), REAL(mus), REAL(K), p );
    
  return(rand_values);
}

//' Samples from a Wishart distribution.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
Rcpp::NumericVector rwish_Rcpp( Rcpp::NumericVector Ts, int b, int p) {
  int pxp = p * p;
  Rcpp::NumericVector K(pxp);
  
  rwish_c( REAL(Ts), REAL(K), &b, &p );
  
  return(K);
}

//' Samples from a G-Wishart distribution according to the algorithm by Dobra and Lenkowski.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List rgwish_Rcpp( const Rcpp::NumericVector G, 
                                 const Rcpp::NumericVector D, 
                                 int b, int p, double threshold) {
  int pxp = p * p, failed;
  Rcpp::NumericVector K( pxp );
  vector<double> inv_Ds( pxp ); 
  vector<double> copy_Ds( pxp ); 
  vector<double> Ts_vectorized( pxp );
  vector<double> Ds_vectorized( pxp );
  
  List return_items;
  
  for( int i = 0; i < p; i++ ){
    for( int j = 0; j < p; j++ ){
      Ds_vectorized[ j * p + i ] = D(j * p + i );
      Ts_vectorized[ j * p + i ] = 0.0;
    }
  }
  
  // This inverts our D matrix, and cholesky decomposes.
  get_Ts( &Ds_vectorized[0], &Ts_vectorized[0], &inv_Ds[0], &copy_Ds[0], &p );
  // This draws from the G Wishart.
  rgwish_c( REAL(G), &Ts_vectorized[0], REAL(K), &b, &p, &threshold, &failed );
  
  return_items["failed"] = failed;
  return_items["K"]      = K;
  
  return(return_items);
}

//' Approximates the G wishart normalizing using an MCMC algorithm.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
Rcpp::NumericVector log_wishart_normalizingConstant_mc_Rcpp( const Rcpp::NumericVector G, 
                                                 const Rcpp::NumericVector nu, const int b, const Rcpp::NumericVector H, const int check_H,
                                                 const int mc_iters, const int p) {
    NumericVector f_T(mc_iters);
    int pxp = p * p;

    vector<int> G_vectorized( pxp );
    vector<int> nu_vectorized( p );
    vector<double> f_vectorized( mc_iters );
    vector<double> H_matrix( pxp );

    int temp_b, temp_check_H, temp_mc_iters, temp_p;
    temp_b = b;
    temp_check_H = check_H;
    temp_mc_iters = mc_iters;
    temp_p = p;

    for( int i = 0; i < p; i++ ){
        nu_vectorized[i] = nu(i);
        for( int j = 0; j < p; j++ ){
          G_vectorized[ j * p + i ] = G(j * p + i );
          H_matrix[ j * p + i ] = H(j * p + i );
        }
    }

    log_exp_mc( &G_vectorized[0], &nu_vectorized[0], &temp_b, &H_matrix[0], &temp_check_H, &temp_mc_iters, &temp_p, &f_vectorized[0] );
    for( int i = 0; i < mc_iters; i++ ){
        f_T(i) = f_vectorized[i];   
    }

    return(f_T);
}



