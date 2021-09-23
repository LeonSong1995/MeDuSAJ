/* Head file 
 * Residual maximum likelihood (REML) to fit the linear mixed model
 * Author: Liyang Song <songliyang@westlake.edu.cn>
 * Adviser: Jian Yang, Xiwei Sun
 * Copyright: Jian Yang 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 *@param X: vector of the fixed component;
 *@param y: vector of the observed value;
 *@param Z: list matrix of the random component;
 *@param maxiter: maximum iterations.
 
 *@example: (in R)
	 y <- abs(rnorm(100,0,1))
	 X <- abs(rnorm(100,0,1))
	 Z1 <- abs(matrix(rnorm(100*2,0,1),ncol=2)) #random component_1
	 Z2<- abs(matrix(rnorm(100*2,0,1),ncol=2))  #random component_2
	 reml(X = X,y = y,Z = list(Z1,Z2),maxiter=1e+3) # return the fixed effect. 
     reml2(X = X,y = y,Z = list(Z1,Z2),maxiter=1e+3) #return the BLUP for random effect. 
 */


#ifndef _MLM_H
#define _MLM_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <math.h>
// #include <omp.h>



using namespace std;
using namespace Eigen;
using namespace Rcpp;

typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;

//REML paramater
vector<eigenMatrix> _A;
eigenMatrix Z; 
eigenMatrix X;
eigenVector y;
eigenMatrix Vi;
eigenMatrix P;
eigenVector varcomp_init;
eigenVector Py;
bool flag_converge = true;
bool flag_inv_Vi = true;
bool flag_inv_P = true;
bool flag_not_itermax = true;
bool flag_EM = false;
vector<double>L_history;

eigenMatrix mat;


//REML function
vector<eigenMatrix> calcu_A(vector<eigenMatrix> &Z, int n, int rindx);
bool calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet,int n, int rindx);
bool calcu_P(eigenMatrix &X,eigenMatrix &Vi, eigenMatrix &P, double &logdet_Xt_Vi_X);
void calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA,int n,int rindx);
bool inverse_H(eigenMatrix &H);
bool ai_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py,eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx);
void em_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx);
double y_center(eigenVector &y,int n);
void constrain_varcmp(eigenVector &y,eigenVector &varcmp,int n, int rindx);
VectorXd reml_iteration(Eigen::VectorXd start, eigenMatrix &X,eigenVector &y, vector<eigenMatrix> &Z, eigenVector &varcmp, int n, int rindx, int maxiter);

#endif
