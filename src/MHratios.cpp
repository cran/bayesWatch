#include "matrix.h"
#include "rgwish.h"
#include "copula.h"
#include "math.h"
#include <Rcpp.h>

extern "C" {
  
  void log_MH_Gupdate( double* log_alpha_ij, int selected_edge_i, 
                       int selected_edge_j, double oldG[], double newG[], double oldK[], double newK[],    
                                                                                               int b, int p )
  {
    int one = 1, dim1 = p + 1, dim = p, pxp = dim * dim;
    double neg_one = -1.0, trace_of_difference = 0.0, log_G_wish_ratio, newK_det, oldK_det;

    double double_b          = (double)b;
    double bminustwo_onehalf = 0.5*(double_b) - 1.0;
    int ij                   = selected_edge_j * dim + selected_edge_i;
//    int jj                   = selected_edge_j * dim1;

    vector<double> copy_newK( pxp ); 
    vector<double> copy_oldK( pxp ); 
    memcpy( &copy_newK[0], newK, sizeof( double ) * pxp );
    memcpy( &copy_oldK[0], oldK, sizeof( double ) * pxp );

    vector<double> new_minus_old_K( pxp );
    memcpy( &new_minus_old_K[0], newK, sizeof( double ) * pxp );
    F77_NAME(daxpy)( &pxp, &neg_one, &oldK[0], &one, &new_minus_old_K[0], &one );
    
    // I'll start by getting the ratio of Wisharts (G-Wisharts without their normalizing constants)
    for(int i = 0; i < dim; i++ ){
      trace_of_difference += new_minus_old_K[ dim1 * i ];
    }

    log_determinant(&copy_newK[0], &newK_det, &p );
    log_determinant(&copy_oldK[0], &oldK_det, &p );

    *log_alpha_ij  = newK_det - oldK_det; 
    *log_alpha_ij *=  bminustwo_onehalf;
    *log_alpha_ij -= 0.5 * trace_of_difference;
    // Now, we calculate the ratio of the normalizing constants.
    // nu_star = b + sum( Gf[,i] * Gf[,j] )
    // Note:       This approximation requires that we count the number of paths of length two from vertex
    //             i to vertex j (they're undirected, so also vice versa).  This is equivalent to counting the
    //             number of instances where we have an edge that goes from i to m AND a vertex from m to j.
    //             This is equivalent to the following product of column vectors in the matrix G.
    double nu_star = b;
    for( int k = 0; k < dim; k++ )
      nu_star += oldG[selected_edge_i * dim + k] * oldG[selected_edge_j * dim + k];
    nu_star = 0.5 * nu_star;
    
    // With this nu_star value (see Letac et al 2017), we can calculate the approximation
    // of the ratio of normalizing constants.
    log_G_wish_ratio = log(2.0) + 0.5 * log(M_PI) + lgammaf( nu_star + 0.5 ) - lgammaf( nu_star );
    
    // Depending on if we added or subtracted an edge, this ratio may flip.
    if( (int)newG[ ij ] == 0 ) log_G_wish_ratio = - log_G_wish_ratio;	
    
    *log_alpha_ij += log_G_wish_ratio;
  }
  
  
  void log_MH_mergesplit( double *log_alpha_ij, double oldK[], 
                          double newK[], int b, int p )
  {
    double neg_one = -1.0, trace_of_difference = 0.0, newK_det, oldK_det;
    int one = 1, dim1 = p + 1, dim = p, pxp = dim * dim;
    
    double double_b         = (double)b;
    double bminustwo_onehalf = 0.5*(double_b) - 1.0;
    
    vector<double> copy_newK( pxp ); 
    vector<double> copy_oldK( pxp ); 
    memcpy( &copy_newK[0], newK, sizeof( double ) * pxp );
    memcpy( &copy_oldK[0], oldK, sizeof( double ) * pxp );
    
    vector<double> new_minus_old_K( pxp );
    memcpy( &new_minus_old_K[0], newK, sizeof( double ) * pxp );
    F77_NAME(daxpy)( &pxp, &neg_one, &oldK[0], &one, &new_minus_old_K[0], &one );
    
    // I'll start by getting the ratio of Wisharts (G-Wisharts without their normalizing constants)
    for(int i = 0; i < dim; i++ ){
      trace_of_difference += new_minus_old_K[ dim1 * i ];
    }
    // This should grab the determinants using Lapack code.
    log_determinant(&copy_newK[0], &newK_det, &p );
    log_determinant(&copy_oldK[0], &oldK_det, &p );
    
    *log_alpha_ij  = newK_det - oldK_det; 
    *log_alpha_ij *= bminustwo_onehalf;
    *log_alpha_ij -= 0.5 * trace_of_difference;
  }
  
  
  
  void log_transition_probability_HMM( double *log_prob, double transition_probabilities[], 
                                       double my_states[], int length_of_vector )
  {
    double diff;
    int length_minus_one = length_of_vector - 1, state_index;
    *log_prob            = 0.0;
    
    for( int index = 0; index < length_minus_one; index++ )
    {   
      diff = my_states[index+1] - my_states[index];
      // Note that the counting of my_states starts at 1, but C requires that we start at zero.
      state_index = (int) my_states[index] - 1;
      if(diff > 0.5){
        *log_prob += log(1.0 - transition_probabilities[state_index]);
      } else {
        *log_prob += log(transition_probabilities[state_index]);
      }
    }
  }
  
  
  void select_edge_from_G_prior( double G[], double g_prior[], int *selected_edge_i, int *selected_edge_j, int p )
  {
    int selected_edge, counter;
    int ip, i, j, ij, dim = p; //, pxp = dim * dim, p1 = dim - 1, p2 = dim - 2, p2x2 = p2 * 2
    int qp = dim * ( dim - 1 ) / 2;
    
    // Count size of notes
    vector<int> size_node( dim, 0 );
    for( i = 0; i < dim; i++ )
    {
      ip = i * dim;
      for( j = 0; j < dim; j++ ) size_node[ i ] += (int)G[ ip + j ];
    }
    
    // For finding the index of selected edge 
    vector<int> index_row( qp );
    vector<int> index_col( qp );
    counter = 0;
    for( j = 1; j < dim; j++ ) {
      for( i = 0; i < j; i++ ) {
        ij = g_prior[ j * dim + i ];
        if( ( ij != 0.0 ) )
        {
          index_row[ counter ] = i;
          index_col[ counter ] = j;
          counter++;
        }
      }
    }
    
    int sub_qp = counter;
    GetRNGstate();
    // Murph Note: Since I am grabbing an edge random and uniformly, this should cancel
    //             in the MH ratio.  If I do something more principled with this in the 
    //             future, I'll need to also return a density value in this method.
    // Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
    selected_edge   = static_cast<int>( unif_rand() * sub_qp );
    *selected_edge_i = index_row[ selected_edge ];
    *selected_edge_j = index_col[ selected_edge ];
    
    PutRNGstate();
}
  
} // End of extern "C"

using namespace Rcpp;

//' Calculates MH ratio for a graph update.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
SEXP log_MH_Gupdate_Rcpp(int selected_edge_i,
                         int selected_edge_j, Rcpp::NumericVector oldG, 
                         Rcpp::NumericVector newG, Rcpp::NumericVector oldK, 
                         Rcpp::NumericVector newK, int b, int p ) {
  double log_alpha_ij;

  log_MH_Gupdate( &log_alpha_ij, selected_edge_i, 
                  selected_edge_j, REAL(oldG), REAL(newG), REAL(oldK), REAL(newK), b, p );
  
  return wrap(log_alpha_ij);
}

//' Calcuates the MH ratio for a merge split on regime vector.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
SEXP log_MH_mergesplit_Rcpp( Rcpp::NumericVector oldK, Rcpp::NumericVector newK, int b, int p ) {
  double log_alpha_ij;
  
  log_MH_mergesplit( &log_alpha_ij, REAL(oldK), 
                     REAL(newK), b, p );
    
  return wrap(log_alpha_ij);
}

//' Calculates probability of new regime vector according to the Markov process.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
SEXP log_transition_probability_HMM_Rcpp( Rcpp::NumericVector transition_probabilities, 
                                          Rcpp::NumericVector my_states, int length_of_vector ) {
  double log_prob;
  
  log_transition_probability_HMM( &log_prob, REAL(transition_probabilities), 
                                  REAL(my_states), length_of_vector );

  return wrap(log_prob);
}

//' Selects an element of G to update according to the prior probability of edge inclusion.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
Rcpp::NumericVector select_edge_from_G_prior_Rcpp( Rcpp::NumericVector G, 
                                                   Rcpp::NumericVector g_prior, 
                                                   int p )
{
  int selected_edge_i, selected_edge_j;
  NumericVector edges(2);
  
  select_edge_from_G_prior( REAL(G), REAL(g_prior), &selected_edge_i, &selected_edge_j, p );
  
  // Note that I need to transition between C, which starts indexing at zero, and R,
  // which starts indexing at 1.
  edges[0] = selected_edge_i + 1;
  edges[1] = selected_edge_j + 1;
  
  return(edges);
}





