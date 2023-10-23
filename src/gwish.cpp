// bayesWatch modifies and uses source code released in the supplementary material of the paper "Computational Aspects Related to Inference in Gaussian Graphical Models With the G-Wishart Prior." by Alex Lenkoski and Adrian Dobra.  This code was used for this package, and placed until the GNU license, by direct permission from the authors.  The files that use this code are newgraph.cpp, graph.h, gwish.cpp.

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "graph.h"
#include "Rcpp.h"

#include <iostream>
using namespace std;

// #include <RcppArmadillo.h>
// #include <boost/math/special_functions/erf.hpp>
// using namespace Rcpp;
// // [[Rcpp::depends(RcppArmadillo)]]


#define REPS_MC 1000

void util_es_to_A(int *es, int *A, int p)
{
  int i,j,k;
  k = 0;
  for(i = 0; i < p - 1; i++)
    {
      for(j = i + 1; j < p; j++)
    {
      A[i * p + j] = es[k];
      A[j * p + i] = es[k];
      k++;
    }
    }
  return;
}

void transpose(int p_i, int p_j, double *A, double *At)
{
  int i, j;

  for(i = 0;i < p_i; i++)
    {
      for(j = 0; j < p_j; j++)
    {
      At[j * p_i + i] = A[i * p_j + j];
    }
    }
  return;
}

void get_upper_triangle(int p,double *A_full,double *A_tri)
{
  int i, j, k;
  k = 0;
  for(i = 0; i < p; i++)
    {
      for(j = i; j < p; j++)
    {
      A_tri[k] = A_full[i * p + j];
      k++;
    }
    }
  return;
}

void set_mat_identity(int p, double *A)
{
  int i;

  for(i = 0; i < p * p; i++) A[i] = 0;
  for(i = 0; i < p; i++) A[i * p + i] = 1;
  return;
}

void invert(int p, double *A, double *A_inv)
{
  double *A_copy;
  char uplo = 'U';
  int info;
  int i;
  A_copy = new double[p * p];
  for(i=0; i < p * p; i++)A_copy[i] = A[i];

  set_mat_identity(p, A_inv);

  //-----  Use LAPACK  -------
  F77_NAME(dposv)(&uplo, &p, &p, A_copy, &p, A_inv, &p, &info FCONE );
  //--------------------------
  delete[] A_copy;

  return;
}

double log_det(int p, double* A) 
{

  int N = p;
  char jobz = 'N';
  char uplo = 'L';
  int lda = p;
  double* eigvals = new double[N];
  double* work = new double[5*N];
  int LWORK=5*N;
  int info;

    double* X = new double[p * p];
    memcpy(X, A, p * p * sizeof(double));
    double logdet = 0.0;

    F77_NAME(dsyev)(&jobz,&uplo,&N, X, &lda,eigvals, work, &LWORK,&info FCONE FCONE );

    if(0==info) 
      {
        for(int i = 0;i < p;i++) logdet += log(fabs(eigvals[i]));
      }
    else
      {//This is really unlikely. I'd like to know it happened, but its not worth shutting down
      }
    delete[] eigvals;
    delete[] work;
    delete[] X;
    return logdet;
}

void mult_mats(int p_i, int p_k, int p_j, double *A, double *B, double *C)
{
  double a=1.0;
  double c=0.0;
  char CblasTrans = 'N';
  // Major Note 12/9/201 by Murph: I was only able to link this version of dgemm with Mayo's avaliable modules.  I believe that this has causes a SWAP of the order.
  // Lenkowski uses row major order with his matrices (see usage of mult_mats below), but I'm pretty sure that this version of dgemm assumes col major order.  So, 
  // I transpose the matrices. 
  //---  Use CBLAS  -----
  // 
  // dgemm	(	character 	TRANSA,integer 	M,integer 	N,integer 	K,double precision 	ALPHA,
  //         double precision, dimension(lda,*) 	A,
  //         integer 	LDA,
  //         double precision, dimension(ldb,*) 	B,
  //         integer 	LDB,
  //         double precision 	BETA,
  //         double precision, dimension(ldc,*) 	C,
  //         integer 	LDC 
  // )	
    
    
  F77_NAME(dgemm)(&CblasTrans,&CblasTrans,&p_i,&p_j,&p_k,&a,A,&p_i,B,&p_k,&c,C,&p_i FCONE FCONE);
  //---------------------

  return;
}


void mult_square_mats(int p, double* A,double* B,double* C)
{

  double a = 1.0;
  double c = 0.0;
  char CblasTrans = 'T';

  //-- Use CBLAS -------
  F77_NAME(dgemm)(&CblasTrans,&CblasTrans,&p,&p,&p,&a,A,&p,B,&p,&c,C,&p FCONE FCONE);
  //--------------------
  return;
}

void chol(int p,double *A)
{
  int i,j;
  char uplo = 'U';
  int info;

  //----Use Lapack-----
  F77_NAME(dpotrf)(&uplo,&p,A,&p,&info FCONE );
  //-------------------

  //--This isn't exactly what we want--
  //--So reformat----------------
  for(i = 0; i < p-1;i++)
    {
      for(j = i + 1; j < p; j++)
    {
      A[i * p + j] = A[j * p + i];
      A[j * p + i] = 0;
    }
    }

  return;
}

void get_cond_matrix(int p,int p_clique,int *clique_ID,
             int *V_ID,double *Full, double *Result)
{

  int i,j;

  double *submatrix_cV;
  double *submatrix_Vc,*submatrix_VV,*submatrix_VVinv;
  double *submatrix_cV_VVinv;
  double *TempResult;
  //----  Start: Form Matrix Objects --------
  //-----  Make the two cross parts  ------
  submatrix_cV=new double[p_clique * (p - p_clique)];
  for(i = 0; i < p_clique; i++)
    {
      for(j = 0;j < p-p_clique; j++)
    {
      // submatrix_cV[i * (p-p_clique) + j]=Full[ clique_ID[i] * p + V_ID[j] ];
      submatrix_cV[j * (p_clique) + i]=Full[ clique_ID[i] * p + V_ID[j] ];
    }
    }
  submatrix_Vc=new double[(p-p_clique)*p_clique];
  transpose(p_clique,p-p_clique,submatrix_cV,submatrix_Vc);
  //----------------------------------------------

  //---Now the inverse of the remaining elements----
  submatrix_VV = new double[(p-p_clique) * (p-p_clique)];
  for(i = 0; i < (p-p_clique); i++)
    {
      for(j=0;j<p-p_clique;j++)
    {
      // submatrix_VV[i*(p-p_clique)+j]=Full[ V_ID[i] * p + V_ID[j] ];
      submatrix_VV[j*(p-p_clique)+i]=Full[ V_ID[i] * p + V_ID[j] ];
    }
    }

  submatrix_VVinv=new double[(p-p_clique) * (p-p_clique)];
  invert(p-p_clique,submatrix_VV,submatrix_VVinv);

  //------  Form the conditioning matrix --------
  submatrix_cV_VVinv=new double[p_clique * (p-p_clique)];
  mult_mats(p_clique,p-p_clique,p-p_clique,
        submatrix_cV,
        submatrix_VVinv,
        submatrix_cV_VVinv);

  mult_mats(p_clique,p-p_clique,p_clique,
        submatrix_cV_VVinv,
        submatrix_Vc,
        Result);
  
  // Murph: I had to change the above method to column-major order, while Lenkowski's
  // next methods all require row-major.  To fix this, I will convert the Result matrix 
  // back to row-major order.
  TempResult=new double[p_clique*p_clique];
  for(i = 0 ; i < p_clique; i++)
  {
    for(j = 0; j < p_clique; j++)
    {
      TempResult[i * p_clique + j] = Result[i * p_clique + j];
    }
  }
  transpose(p_clique, p_clique,TempResult,Result);
  
  //---------------------------------------------
  delete[] submatrix_cV;
  submatrix_cV=0;
  delete[] submatrix_Vc;
  submatrix_Vc=0;
  delete[] submatrix_VV;
  submatrix_VV=0;
  delete[] submatrix_VVinv;
  submatrix_VVinv=0;
  delete[] submatrix_cV_VVinv;
  submatrix_cV_VVinv=0;
  delete[] TempResult;
  TempResult = 0;

  return;
}//get_cond_matrix


//This algorithm takes a matrix target and updates
//it according to the edges in the graph thegraph.
void copy_mats(int p_i, int p_j, double *A, double *B)
{
  int i,j;
  for(i = 0 ; i < p_i; i++)
    {
      for(j = 0; j < p_j; j++)
    {
      B[i * p_j + j] = A[i * p_j + j];
    }
    }
  return;
}

void make_sub_mat_dbl(int p, int p_sub, int *sub, double *A, double *B)
{
  int i,j;
  for(i = 0; i < p_sub; i++)
    {
      for(j = 0; j < p_sub; j++)
    {
      B[i * p_sub + j] = A[  sub[i] * p + sub[j] ];
    }
    }

  return;
}

void make_sub_mat_int(int p, int p_sub, int *sub, int **A, int *B)
{
  int i,j;
  for(i = 0; i < p_sub; i++)
    {
      for(j = 0; j < p_sub; j++)
    {
      B[i * p_sub + j] = A[ sub[i] ][ sub[j] ];
    }
    }
  return;
}

void get_complementary_set(int p, int p_clique, int *clique_ID, int *V_ID)
{
  int i,j,k,temp_in_clique;
  k = 0;
  for(i = 0; i < p; i++)
    {
      temp_in_clique = 0;
      for(j = 0; j < p_clique; j++)
    {
      if(i == clique_ID[j])
        {
          temp_in_clique = 1;
        }
    }
      if(temp_in_clique == 0) 
    {
      V_ID[k] = i;
      k++;
    }
    }
  return;
}

// This code taken from something Adrian Dobra sent me
int choose(int n, int m)
{
   if(m > n) return 0;
   if(m == 0) return 1;
   if(m == 1) return n;
   if(m == n) return 1;
   if(m == (n-1)) return n;
   return(choose(n-1,m)+choose(n-1,m-1));
}

double logchoose(int n, int m)
{
  int i;
  double s = 0;
  for(i = 0; i < n;i++)s += log(i + 1);
  for(i = 0; i < m;i++) s -= log(i + 1);
  for(i = 0; i < n - m;i++) s -= log(i + 1);
  return(s);
}
//This is an increment for a combinations operation
//It returns 1 if there is a another combination after the current one
//And zero if this is the final combination
int combinations_increment(int n, int r, int *c)
{
  int i, k;
  for(i = 0; i < r; i++)
    {
      if(c[i] != n - i - 1)
    {
      c[i]++;
      for(k = i - 1; k > -1; k--) c[k] = c[k + 1] + 1;
      return 1;
    }
    }
  return 0;
}

//This initializes a vector to the first combination in a sequence
void combinations_init(int n, int r, int *c)
{
  int i;
  for(i = 0; i < r; i++) c[i] = r - i - 1;
}

int test_add_var(int *A, int p, int *var_list, int list_length, int prop_var)
{
  int is_included = 1;
  int i;
  for(i =0; i < list_length; i++)
    {
      if( A[ var_list[i] * p + prop_var] == 0)
    {
      is_included = 0;
      break;
    }
    }
  return(is_included);
}//test_add_var

int is_subclique(int *var_list, int list_size, int *clmat, int *cldims, int p)
{
  int i,j,k;
  int temp_in;
  int temp_subclique;
  int temp_var_found;
  int ee = p * (p - 1) / 2;
  temp_subclique = 0;
  for(i = 0; i < ee; i++)//go through each potential clique
    {
      if(cldims[i] > 0)
    {
      temp_in = 1;
      //--------- move through each vertex in clique -------------------
      for(j = 0; j < list_size; j++)
        {
          temp_var_found = 0;
          //----- now see if this variable is in this clique --------
          for(k = 0; k < cldims[i]; k++)
        {
          if(clmat[i * p + k] == var_list[j]) temp_var_found = 1;
        }
          if(temp_var_found == 0)temp_in = 0;
          //--------------------------------------------------------
        }
      //---------------------------------------------------------------
      if(temp_in == 1) temp_subclique = 1;
    }
    }
  return(temp_subclique);
}//is_subclique

void add_clique(int *var_list_curr, int curr_length, int *clmat, int *cldims, int p)
{
  int i;
  int pos = -1;
  int ee = p * (p - 1) / 2;

  for(i = 0; i < ee; i++)
    {
      if(cldims[i] == 0) 
    {
      pos = i;
      break;
    }
    }

  cldims[pos] = curr_length;
  for(i = 0; i < curr_length; i++) clmat[pos * p + i] = var_list_curr[i];
  return;
}//add_clique

void list_can_augment(int *A, int p, int *var_list_curr, int curr_length, int *clmat, int *cldims)
{
  int *list_augment;
  int i,cc;
  int start_point;
  int can_add;
  int maximal = 1;
  int sc_check;
  start_point = var_list_curr[ curr_length - 1] + 1;
  for(i = start_point; i < p; i++)
    {

      can_add = test_add_var(A, p, var_list_curr, curr_length, i);
      if(can_add)
    {
      maximal = 0;
      list_augment = new int[curr_length + 1];
      for(cc = 0; cc < curr_length; cc++) list_augment[cc] = var_list_curr[cc];
      list_augment[curr_length] = i;
      list_can_augment(A, p, list_augment, curr_length + 1, clmat, cldims);
      delete[] list_augment;
    }
    }
  if(maximal)
    {
      sc_check = is_subclique(var_list_curr,curr_length, clmat, cldims, p);
      if(!sc_check)
    {
      add_clique(var_list_curr, curr_length, clmat, cldims, p);
    }
  }
}//list_can_augment

int get_cliques(int *A, int p, int *clmat, int *cldims)
{

  int ee = p * (p - 1) /2;
  int i,k;
  int *temp_init_array;
  temp_init_array = new int[1];
  for(i = 0; i < ee * p; i++) clmat[i] = -1;
  for(i = 0; i < ee; i++) cldims[i] = 0;
  for(i = 0; i < p; i++)
    {
      temp_init_array[0] = i;
      list_can_augment(A, p, temp_init_array, 1, clmat, cldims);
    }

  delete[] temp_init_array;
  k = 0;
  for(i = 0; i < ee; i++)
    {
      if(cldims[i] > 0) k++;
    }
  return(k);
}


void IPF_MLE(int *CliqueMat, int *CliqueDims, int totcliques, 
         double *target, int p, double thresh, int *nonconverge_flag)
{

  int i, j;
  int *clique_ID, *V_ID;
  double tempdiff, maxdiff;

  int p_clique;
  int clique_num;
  double *final, *final_last;
  double *submatrix_cc, *submatrix_ccINV;
  double *submatrix_cond;
  int count = 0;
  // int nonconverge_flag = 0;

  //-----------------------------------------
  //Initialize the matrices that we write over
  //each iteration to the identity
  final = new double[p * p];
  final_last = new double[p * p];
  set_mat_identity(p, final);
  //-----------------------------------------

  maxdiff = thresh + 1;

  while(maxdiff > thresh)
    {
      count += 1;
      if(count%50000==0){
      }
    
      if(count > 100000){
        *nonconverge_flag += 1;
        break;
      }
      copy_mats(p, p, final, final_last);
      for(clique_num = 0; clique_num < totcliques; clique_num++)
    {
      //*******Just recording Indices***********
      p_clique = CliqueDims[clique_num];
      clique_ID = new int[p_clique];
      V_ID = new int[p-p_clique];
      for(i = 0; i < p_clique; i++) clique_ID[i] = CliqueMat[clique_num * p + i];
      get_complementary_set(p, p_clique, clique_ID, V_ID);
      //******************************************

      //******  Start: Making Submatrices  ************
      //---First make the part specific to the clique-----
      submatrix_cc = new double[p_clique * p_clique];
      submatrix_ccINV = new double[p_clique * p_clique];
      make_sub_mat_dbl(p, p_clique, clique_ID, target, submatrix_ccINV);
      invert(p_clique, submatrix_ccINV, submatrix_cc);//And form the inverse, which is what we want.
      //--------------------------------------------------

      //-----   Now form the part that's based ------
      //----- on the entries outside of the clique --
      submatrix_cond = new double[p_clique * p_clique];
      get_cond_matrix(p, p_clique, clique_ID, V_ID, final, submatrix_cond);
      //---------------------------------------------

      //--- Add the two matrices  -------------------
      for(i = 0; i < p_clique * p_clique; i++) submatrix_cc[i] = submatrix_cc[i] + submatrix_cond[i];
      //---------------------------------------------
          //********  End Forming Matrices ********************

      //---Update the result matrix for the clique---
      for(i = 0; i < p_clique; i++)
        {
          for(j = 0; j < p_clique; j++)
        {
          final[clique_ID[i] * p + clique_ID[j]] = submatrix_cc[j * p_clique + i];
        }
        }
      //---------------------------------------------
      //----- Clean up from updating this clique-----
      delete[] clique_ID;
      clique_ID = 0;
      delete[] V_ID;
      V_ID = 0;
      delete[] submatrix_cc;
      submatrix_cc = 0;
      delete[] submatrix_ccINV;
      submatrix_ccINV = 0;
      delete[] submatrix_cond;
      submatrix_cond = 0;
      //----------------------------------------------

    }//End Loop through cliques
      //--We've now looped through all cliques in the graph and updated---

      //-----Now test for convergence.----------
      tempdiff = 0;maxdiff = 0;
      for(i = 0; i < p; i++)
    {
      for(j = 0; j < p; j++)
        {
          tempdiff=fabs(final[i*p + j] - final_last[i*p+j]);
          if(tempdiff > maxdiff)maxdiff = tempdiff;
        }
    }
      //-----------------------------------------

    }//End While(maxdiff>thresh)

  //-- Update the supplied array with the inverse----
  invert(p,final,target);
  //----------------------------
  
  
  
  delete[] final; 
  final = 0;
  delete[] final_last;
  final = 0;
  return;
}//IPF_MLE


double gwish_nc_complete(int delta, int p, double *D)
{
  double I;
  double d,c,a,g;
  double dblDelta;//Recasting the inputs makes life easier below
  double dblP;
  int i;

  dblDelta=delta;
  dblP=p;

  //----    Calculate Each Piece --------
  d = ((dblDelta+dblP-1) / 2) * log_det(p,D);

  c = (dblP * (dblDelta + dblP - 1)) / 2 * log(2);

  a = (dblDelta + dblP - 1) / 2;
  g = dblP * (dblP - 1) / 4 * log(M_PI);
  for(i = 0;i < dblP;i++) g += lgamma(a - (double)i / 2);
  I = -d + c + g;
  //------------------------------------
  return(I);
}


//----- Returns Log of C_G---------------
//  The constant term in the expression  
//  Of the normalizing constant for a 
//  non-decomposable G-Wishart Variate
//---------------------------------------
double gwish_logC(int *A, int delta, double *T, int p)
{

  double *nu;//Made double due to casting issues
  double *b;
  int i,j,qq;
  double a;
  double d;
  double g;
  double tot;
  double temp1;
  double dblDelta;

  nu = new double[p];
  b = new double[p];
  dblDelta = delta;

  for(i = 0; i < p; i++)nu[i]=0.0;
  for(i = 0; i < p; i++)b[i]=0.0;
  for(i = 0; i < p; i++)
    {
      for(j = i + 1; j < p; j++) nu[i] += (double)A[i * p + j];
      b[i] = nu[i] + 1.0;
      for(qq = 0; qq < i; qq++)
    {
      b[i] += (double)A[qq * p + i];//This embeds k_i in b_i
    }
    }
  tot=0;
  for(i = 0; i < p; i++)
    {
      temp1 = (dblDelta + 2 * nu[i]) / 2;
      a = temp1 * log(2) + nu[i] / 2 * log(M_PI);
      g = lgamma((dblDelta + nu[i]) / 2);
      d=(dblDelta + b[i] - 1) * log(T[i * p + i]);
      tot += a + g + d;

    }
  delete[] nu;
  delete[] b;

  return(tot);
}//gwish_logC


//----- Returns Laplace Approximation -------
// For the normalizing constant of a 
// G-Wishart Variate

// NOTE:  This already assumes that the scale
// Matrix D has been properly completed
// via the IPF such that
// D^{-1}\in M^+(G)

// ALSO:  The Laplace approximation relies
// On the use of the mode of the G-Wishart
// distribution with params delta and D
// However, in this case
// \hat{K}=(delta-2)*D^{-1}
//  Furthermore, all calculations
//  Involving \hat{K} can be derived
//  From the IPF'd D, thus allowing
//  us to invert one less matrix.
//-----------------------------------------
double gwish_norm_laplace(int p, int *A, int delta, double *D)
{

  int i,j,k,l,ee;
  double *M_1,*M_2,*M_3;
  double *H;
  double temp;
  int q;
  double a,pp,hh;
  int *label_map;
  double dblDelta;

  dblDelta = delta; //convenient for latter.

  //----First Find how many terms will be in the hessian, called ee---
  ee = p;
  for(i = 0;i < p-1; i++)for(j = i + 1; j < p; j++)if(A[i * p + j] == 1) ee++;
  //--------------------------------------------------------------------

  //----  Initialization  -----------
  H = new double[ee * ee];
  M_1 = new double[p * p];
  M_2 = new double[p * p];
  M_3 = new double[p * p];
  label_map=new int[2 * ee];
  //---------------------------------

  //------   Lookup Table   -------------
  //    We need to form a lookup table
  //    That maps between elements
  //    of the matrix K and the Hessian H.
  for(i = 0;i < ee * ee;i++) H[i] = 0;

  for(i = 0; i < p; i++) label_map[i * 2] = label_map[i * 2 + 1] = i;

  q = p;//q is our carrier as we search through the graph for edges

  for(i = 0; i < p-1; i++)
    {
      for(j = i + 1; j < p; j++)
    {
      if(A[i * p + j] == 1)
        {
          label_map[q * 2] = i;
          label_map[q * 2 + 1] = j; 
          q++;
        }
    }
    }
  //--------  End Lookup Table ---------

  //---------Fill in Hessian--------------------
  // Note that the mode of the G-Wishart with scale D
  // is (delta-2)*D^{-1}, and since the Hessian is actually
  // evaluated at the inverse of the mode I never actually
  // form K_hat.
  for(i = 0; i < ee; i++)
    {
      for(j = i; j < ee; j++)
    {

      //-----Form the left-hand side------
      for(l = 0;l < p * p; l++) M_1[l] = D[l];
      for(l = 0; l < p * p; l++) M_2[l] = 0;
      M_2[ label_map[i * 2] * p + label_map[i * 2 + 1] ] = 1;
      M_2[ label_map[i * 2 + 1] * p + label_map[i * 2] ] = 1;
      mult_square_mats(p,M_1,M_2,M_3);
      //----------------------------------

      //---  The left three ------
      for(l = 0;l < p * p; l++) M_1[l] = D[l];
      mult_square_mats(p,M_3,M_1,M_2);
      //--------------------------

      //--- And the final matrix on the right ---
      for(l = 0; l<p * p; l++) M_1[l] = 0;
      M_1[label_map[j * 2] * p + label_map[j * 2 + 1]] = 1;
      M_1[label_map[j * 2 + 1] * p + label_map[j * 2]] = 1;
      mult_square_mats(p,M_2,M_1,M_3);
      //--------------------------------------

      //--- Now take the trace of this product ---
      temp = 0;
      for(k = 0; k < p; k++) temp += M_3[k * p + k];
      //-------------------------------------------

      //-- Fill in the entry of the hessian  ----
      //the delta-2 adusts for using D and not \hat{K}^{-1}=D/(delta-2)
      H[i * ee + j] = -temp / (delta-2) / (delta-2);
      H[j * ee + i] = -temp / (delta-2) / (delta-2);
      //------------------------------------------

    }//j
    }//i
  //----- Finished filling in the Hessian-----------


  //------  Calculate relevant quantities  -------
  //  Note that since \hat{K}=(delta-2)*D^{-1}
  //  The <K,D> term in P[K] is just p*(delta-2)
  //  which is handy.  Also, because of this
  //  log(det(K))=-log(det(D))+p*log(delta-2)
  //  This means that I never actually form the 
  //  matrix K_hat, but instead just use the IPF'd
  //  scale parameter.
  //-----------------------------------------------
  a = ee / 2 * log(2 * M_PI);

  double pp_det = (dblDelta-2) / 2 * ( -log_det(p,D) + p * log(dblDelta - 2) );
  double pp_kern = -p * (dblDelta-2) / 2;
  pp = pp_det + pp_kern;

  for(i = 0; i < ee * ee; i++) H[i] = -(dblDelta-2) / 2 * H[i];    
  hh = .5 * log_det(ee,H);

  //------   Clean Up  --------------------
  delete[] M_1; delete[] M_2; delete[] M_3;
  delete[] H;
  delete[] label_map;
  //---------------------------------------

  return(a - hh + pp);

}

// Murph: again, I don't want the normalizing constants.
double gwish_exact_posterior(int p, int delta, int n,  double *D_post) // double *D_prior,
{
//I_prior,
  double  I_post;
  //--------------------------
  //I_prior = gwish_nc_complete(delta, p, D_prior);
  I_post = gwish_nc_complete(delta+n, p, D_post);
  //--------------------------

  return(I_post );//- I_prior
}

double gwish_calculateLogPosterior(LPGraph graph, double *D_post, int delta,int n,int *nonconverge_flag)
{
  // Comments added by Murph 2021
  double mypost = 0;
  int i, p;
  double *sub_D_post;
  // Note that I am removing the prior calculations.  In Lenkowski and Dobra, they use this method to calculate
  // the marginal posterior -- I just want the posterior normalizing constant, given the kernel of a g wishart posterior.
  //double *sub_D_prior;
  int sub_p;
  int sub_ee;
  int *sub_A;
  int *CliqueMat, *CliqueDims, totcliques;
  p = graph->nVertices;
  //obtain the prime components and the separators
  graph->GetMPSubgraphs();
  // graph->FindCliqueTree();
  for(i=0;i<graph->nCliques;i++)
    {
      sub_p = graph->CliquesDimens[i];
      //----------  Form Sub-Matrices  -----------
      sub_D_post = new double[sub_p*sub_p];
      //sub_D_prior = new double[sub_p*sub_p];

      //make_sub_mat_dbl(p, sub_p, graph->Cliques[i], D_prior, sub_D_prior);
      make_sub_mat_dbl(p, sub_p, graph->Cliques[i], D_post, sub_D_post);
      //-------------------------------------------

      //------------------------------------------------
      //rintf("Clique %d: %d\n",graph->IsClique(graph->Cliques[i],graph->CliquesDimens[i]));
      if(graph->IsClique(graph->Cliques[i], graph->CliquesDimens[i]))
        {
          mypost += gwish_exact_posterior(sub_p, delta, n, sub_D_post);
        }
      else //otherwise do the approximation
        {
      //-------  Form the SubAdjacency Matrix  -------------------
      sub_ee = sub_p * (sub_p - 1) / 2;
      sub_A = new int[sub_p * sub_p];
      make_sub_mat_int(p, sub_p, graph->Cliques[i], graph->Edge, sub_A);
      //----------------------------------------------------------

      //-------  Get the Clique Decomposition --------------------
      CliqueMat = new int[sub_ee * sub_p];
      CliqueDims = new int[sub_ee];	
      totcliques = get_cliques(sub_A, sub_p, CliqueMat, CliqueDims);
      //----------------------------------------------------------

      //-------  IPF the D matrices Accordingly ------------------
      //IPF_MLE(CliqueMat, CliqueDims, totcliques, sub_D_prior, sub_p, .000001);
      IPF_MLE(CliqueMat, CliqueDims, totcliques, sub_D_post, sub_p, .00001, nonconverge_flag);
      //----------------------------------------------------------

      //----  Add the Posterior Normalizing Constant  ------------
      mypost += gwish_norm_laplace(sub_p, sub_A, delta + n, sub_D_post);


      //---- Subtract the Prior Normalizing Constant -------------
          // Murph: again, I don't want the normalizing constants.
      // mypost -= gwish_norm_mc_alt(sub_A, delta, sub_D_prior, sub_p, REPS_MC, stream);
      //------------------------------------------------------------
      delete[] CliqueMat;
      delete[] CliqueDims;
      delete[] sub_A;
    }
      delete[] sub_D_post; //delete[] sub_D_prior;
    }//Loop through the prime components
  for(i = 0; i < graph->nSeparators; i++)
    {
      sub_p = graph->SeparatorsDimens[i];
      //----------  Form Sub-Matrices  -----------
      sub_D_post = new double[sub_p * sub_p];
      // Murph: again, I don't want the normalizing constants.
      // sub_D_prior = new double[sub_p * sub_p];
      //make_sub_mat_dbl(p, sub_p, graph->Separators[i], D_prior,sub_D_prior);
      make_sub_mat_dbl(p, sub_p, graph->Separators[i], D_post,sub_D_post);
      //-------------------------------------------

      //------------------------------------------
      //   GGMs are graphical, hence the 
      //   separators should always be complete
      // Murph: I updated this function to ignore the prior information.
      mypost -= gwish_exact_posterior(sub_p,delta,n,sub_D_post);
      //--------------------------------------------
      // Murph: again, I don't want the normalizing constants.
      delete[] sub_D_post; // delete[] sub_D_prior;
    }//Loop through the separators
  return(mypost);

}//gwish_calculateLogPosterior



void gwish_get_psi_from_K(int p, int delta, double *D, double *K, double *Psi)
{

  double *T;
  double *T_inv;
  double *Phi;
  T = new double[p * p];
  T_inv = new double[p * p];
  Phi = new double[p * p];

  invert(p, D, T);
  chol(p, T);
  invert(p, T, T_inv);
  copy_mats(p, p, K, Phi);
  chol(p, Phi);
  mult_square_mats(p, Phi, T_inv, Psi);

  delete[] T; delete[] T_inv; delete[] Phi;
  return;
}

double gwish_get_f_Tsq(int p, int *A, double *Psi)
{
  double tot = 0;
  int i;
  int j;
  for(i = 0; i < p - 1; i++)
    {
      for(j = i + 1; j < p; j++)
    {
      if(A[i * p + j]==0)
        {
          tot += Psi[i * p + j] * Psi[i * p + j];
        }
    }
    }
  return(tot);
}
int FindDecomposableNeighbors(LPGraph graph,int* decNeighbors)
    {
        int i,j,k;
        int ndecNeighbors = 0;
        int p = graph->nVertices;

        k = 0;
        for(i=0;i<p-1;i++)
        {
            for(j=i+1;j<p;j++)
            {
                //reverse the edge
                graph->Edge[i][j] = 1-graph->Edge[i][j];
                graph->Edge[j][i] = 1-graph->Edge[j][i];

                if(graph->IsDecomposable())
                {
                    decNeighbors[k] = 1;
                    ndecNeighbors++;
                }
                else
                {
                    decNeighbors[k] = 0;
                }
                k++;

                //reverse back the edge
                graph->Edge[i][j] = 1-graph->Edge[i][j];
                graph->Edge[j][i] = 1-graph->Edge[j][i];
            }
        }

        return(ndecNeighbors);
    }


using namespace Rcpp;

//' Calculates the normalizing constant for a G-Wishart with fixed graph structure using a Laplace approximation.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List log_normalizing_g_wishart_posterior_laplace( NumericMatrix graph,
                            NumericMatrix D_post, int Delta, int n, int p) {
  // Method written by Murph to access Lenkowski's and Dobra's functions.
  List return_items;
  int i,j;
  double s;
  int nonconverge_flag = 0;
  // I'll start be making a graph structure.
  LPGraph new_graph = new Graph;
  new_graph->nVertices = p;
  new_graph->InitGraph(p);
  for(i=0;i<p-1;i++)
  {
      for(j=i+1;j<p;j++)
      {
          new_graph->Edge[i][j] = graph(i,j);
          new_graph->Edge[j][i] = graph(i,j);
      }
  }
  s = gwish_calculateLogPosterior(new_graph, REAL(D_post), Delta, n, &nonconverge_flag);
  delete new_graph;
  
  if(nonconverge_flag > 0){
    nonconverge_flag = 1;
  }
  
  return_items["log_posterior"] = s;
  return_items["nonconverge_flag"] = nonconverge_flag;
  return return_items;
  
}

