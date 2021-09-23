#include "INV.hpp"
#include "MLM.h"


//[[Rcpp::export]]
Eigen::MatrixXd ecor(Eigen::MatrixXd mat)
{
	Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
	Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(mat.rows() - 1);
	return (cov);
}



//[[Rcpp::export]]
std::vector<Eigen::MatrixXd> reml(Eigen::VectorXd start, Eigen::MatrixXd &X, Eigen::VectorXd &y, std::vector<Eigen::MatrixXd> &Z, int maxiter)
{
	flag_converge = true;
	int rindx = (1+Z.size());
	int n = X.rows();
	Eigen::VectorXd varcmp(rindx);
	// int test =1;
	varcmp = reml_iteration(start, X, y, Z, varcmp,n, rindx,maxiter);

	eigenMatrix tX_VI_X;
	eigenMatrix tX_VI_y;
	tX_VI_X = X.transpose() * Vi * X;
	tX_VI_y = X.transpose() * Vi * y;

	eigenMatrix b;
	eigenMatrix sd;

	// int converge = 0;
	// int inv_vi = 0;
	// int inv_p = 0;
	// int not_itermax = 0;

	if (flag_converge)
	{
		// converge = 1;

		b = tX_VI_X.inverse() * tX_VI_y;
		sd = tX_VI_X.inverse();
	}

	// if(flag_inv_Vi) inv_vi = 1;
	// if(flag_inv_P) inv_p = 1;
	// if(flag_not_itermax) not_itermax =1;


	L_history.clear();


	std::vector<Eigen::MatrixXd> r;
	r.resize(4);
	r[0] = b;
	r[1] = sd;
	r[2] = varcmp;
	r[3] = Vi;



	// NumericVector out = NumericVector::create(b);
	return (r);

}





//[[Rcpp::export]]
std::vector<Eigen::MatrixXd>  reml2(Eigen::VectorXd start, Eigen::MatrixXd &X, Eigen::VectorXd &y, std::vector<Eigen::MatrixXd> &Z, int maxiter)
{
	flag_converge = true;
	int rindx = (1+Z.size());
	int n = X.rows();
	Eigen::VectorXd varcmp(rindx);
	// int test = 1;
	varcmp = reml_iteration(start, X, y, Z, varcmp,n, rindx,maxiter);


	eigenMatrix tX_VI_X;
	eigenMatrix tX_VI_y;
	tX_VI_X = X.transpose() * Vi * X;
	tX_VI_y = X.transpose() * Vi * y;

	eigenMatrix b;
	eigenMatrix sd;

	// int converge = 0;
	// int inv_vi = 0;
	// int inv_p = 0;
	// int not_itermax = 0;

	if (flag_converge)
	{
		// converge = 1;

		b = tX_VI_X.inverse() * tX_VI_y;
		sd = tX_VI_X.inverse();
	}
	eigenMatrix fix(2,b.size());
	fix.row(0) = b;
	fix.row(1) = sd;

	eigenMatrix Q;
	eigenMatrix temp;
	eigenVector flattened;
	eigenVector flattened_Xita;
	std::vector<Eigen::VectorXd> e;
	std::vector<Eigen::VectorXd> xita;
	e.resize(rindx-1);
	xita.resize(rindx-1);

	//AUP
	// eigenMatrix MRH;
	// double k;
	// MRH = y.transpose() * Q * Z[i];
	// cout << (MRH * MRH.transpose())[0,0] << endl;
	// k = (MRH * MRH.transpose())[0,0];
	// k = sqrt((varcmp[i]*(Z[i].cols()-1))/(MRH * MRH.transpose()));
	// e[i]  = (k * Z[i].transpose() * Q * y)/ sqrt(Z[i].cols());

	eigenMatrix v_tZ_Q;

	if (flag_converge)
	{
		Q = Vi - Vi * X * tX_VI_X.inverse()  * X.transpose() *Vi;
		for (int i=0; i<(rindx-1);i++)
		{
			v_tZ_Q = varcmp[i] * Z[i].transpose() * Q;
			//BLUP
			e[i] = (v_tZ_Q * y)/Z[i].cols();
			//BLUP se
			temp = varcmp[i] *v_tZ_Q * Z[i];
			xita[i] = temp.diagonal().array().sqrt() / Z[i].cols();
		}


		// concatenate vector
		int len = 0;
		for (auto const &v : e) len += v.size();

		flattened.resize(len); flattened_Xita.resize(len);

		int offset = 0;
		int size;
		for (int i=0; i<(rindx-1);i++)
		{
			size = e[i].size();
			flattened.middleRows(offset,size) = e[i];
			flattened_Xita.middleRows(offset,size) = xita[i];
			offset += size;
		}
	}else{
		flattened.resize(1); flattened_Xita.resize(1);
		flattened << 0; flattened_Xita << 0;
	}

	L_history.clear();

	eigenMatrix res(2,flattened.size());
	res.row(0) = flattened;
	res.row(1) = flattened_Xita;

	std::vector<Eigen::MatrixXd> r;
	r.resize(3);
	r[0] = res;
	r[1] = varcmp;
	r[2] = fix;

	return(r);

}




vector<eigenMatrix> calcu_A(vector<eigenMatrix> &Z, int n, int rindx)
{
	vector<eigenMatrix> _A;
	_A.resize(rindx);
	for (int i=0; i<(rindx-1);i++)
	{
		_A[i] = (Z[i] * Z[i].transpose())/Z[i].cols();
	}
	_A[rindx-1] = eigenMatrix::Identity(n, n);
	return _A;

}



bool calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet,int n, int rindx)
{
	int i = 0, j = 0;
	Vi = eigenMatrix::Zero(n, n);
	for (i = 0; i < rindx; i++) {
		Vi += (_A[i]) * prev_varcmp[i];
	}
	INVmethod method_try;
	method_try = INV_LLT;

	int rank = 0;
	bool ret = true;
	if (!SquareMatrixInverse(Vi, logdet, rank, method_try))
	{
		//cout << "warning:the variance-covariance matrix V is invertible, a small positive value is added to the diagonals" << endl;
		double d_buf = Vi.diagonal().mean() * 1e-4;
		for (j = 0; j < n; j++) Vi(j, j) += d_buf;
		if (!SquareMatrixInverse(Vi, logdet, rank, method_try)) {
			std::cout << "Still can't be inversed" << endl;
			ret = false;
		}
	}
	return ret;
}


bool calcu_P(eigenMatrix &X,eigenMatrix &Vi, eigenMatrix &P, double &logdet_Xt_Vi_X)
{
	eigenMatrix Vi_X;
	eigenMatrix Xt_Vi_X_i;
	Vi_X = Vi * X;
	Xt_Vi_X_i = X.transpose() * Vi_X;
	int rank = 0;
	INVmethod method_try;
	method_try = INV_LLT;
	if (!SquareMatrixInverse(Xt_Vi_X_i, logdet_Xt_Vi_X, rank, method_try))return false;

	P = Vi - Vi_X * Xt_Vi_X_i * Vi_X.transpose();
	return true;
}


void calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA,int n,int rindx)
{
	tr_PA.resize(rindx);
	double *d_bufs = new double[n];
	for (int i = 0; i < rindx; i++)
	{
		double d_buf = 0.0;
		memset(d_bufs, 0, n * sizeof(double));
		#pragma omp parallel for
		for (int k = 0; k < n; k++) {
			for (int l = 0; l < n; l++) d_bufs[k] += P(k, l)*(_A[i])(k, l);
		}
		for (int k = 0; k < n; k++) {
			d_buf += d_bufs[k];
		}
		tr_PA(i) = d_buf;
	}
	delete[] d_bufs;
}



bool inverse_H(eigenMatrix &H)
{
	double d_buf = 0.0;
	INVmethod method_try;
	method_try = INV_LLT;
	int rank = 0;
	if (!SquareMatrixInverse(H, d_buf, rank, method_try)) return false;
	else return true;
}



bool ai_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py,eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx,double step)
{

	Py = P * y;
	eigenMatrix APy(n, rindx);
	#pragma omp parallel for
	for (int i = 0; i < rindx; i++) {
		(APy.col(i)) = (_A[i]) * Py;
	}

	// Calculate Hi
	eigenVector R(rindx);
	eigenMatrix Hi(rindx, rindx);
	#pragma omp parallel for
	for (int i = 0; i < rindx; i++) {
		R(i) = (Py.transpose()*(APy.col(i)))(0, 0);
		eigenVector cvec = P * (APy.col(i));
		Hi(i, i) = ((APy.col(i)).transpose() * cvec)(0, 0);
		for (int j = 0; j < i; j++) Hi(j, i) = Hi(i, j) = ((APy.col(j)).transpose() * cvec)(0, 0);
	}
	Hi = 0.5 * Hi;

	// Calcualte tr(PA) and dL
	eigenVector tr_PA;
	calcu_tr_PA(P, tr_PA,n, rindx);
	R = -0.5 * (tr_PA - R);

	// Calculate variance component
	if (!inverse_H(Hi))
	{
		cout << "the information matrix is not invertible.Transfer to EM-REML(slow)" << endl;
		return false;
	}
	eigenVector delta(rindx);
	delta = Hi * R;
	varcmp = prev_varcmp + step *delta;
	return true;
}


void em_reml(eigenVector &y,eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, int n, int rindx)
{
	eigenVector tr_PA;
	calcu_tr_PA(P, tr_PA,n, rindx);
	Py = P * y;

	eigenVector R(rindx);
	#pragma omp parallel for
	for (int i = 0; i < rindx; i++)
	{
		R(i) = (Py.transpose()*(_A[i]) * Py)(0, 0);
		varcmp(i) = prev_varcmp(i) - prev_varcmp(i) * prev_varcmp(i) * (tr_PA(i) - R(i)) / n;
	}
}


double y_center(eigenVector &y,int n)
{
	eigenVector _y_center(n);
	_y_center.setConstant(y.mean());
	_y_center = y - _y_center;
	return(_y_center.squaredNorm());
}


void constrain_varcmp(eigenVector &y,eigenVector &varcmp,int n, int rindx)
{
	double constr_scale = 1e-6;
	for (int i = 0; i < rindx; i++) {if (varcmp[i] < 0) varcmp[i] = constr_scale * y_center(y,n);}
}



VectorXd reml_iteration(Eigen::VectorXd start, eigenMatrix &X,eigenVector &y, vector<eigenMatrix> &Z, eigenVector &varcmp, int n, int rindx,int maxiter)
{
	double logdet = 0.0, logdet_Xt_Vi_X = 0.0, prev_lgL = -1e20, lgL = -1e20, dlogL = 1000.0, step = 0.316;
	int L_size;

	eigenVector prev_varcmp(rindx);
	prev_varcmp = start;

	_A = calcu_A(Z,n, rindx);
	L_history.push_back(lgL);

	for (int iter = 0; iter <= maxiter; iter++)
	{
		if (!calcu_Vi(Vi, prev_varcmp, logdet,n, rindx))
		{
			cout << "REML ERROR!:V matrix is not positive. Switch to multiple linear regression!" << endl;
			flag_converge = false;
			flag_inv_Vi = false;
			return varcmp;
		}
		if (!calcu_P(X,Vi, P, logdet_Xt_Vi_X))
		{
			cout << "REML ERROR!: the X^t * V^-1 * X matrix is not invertible,please check the Signature Genes. Switch to multiple linear regression!" << endl;
			flag_converge = false;
			flag_inv_P = false;
			return varcmp;
		}

		//initialized with EM-REML
		if (iter == 0)
		{
			em_reml(y,P, Py, prev_varcmp, varcmp, n, rindx);
		}
		else
		{
			if (flag_EM) em_reml(y,P, Py, prev_varcmp, varcmp, n, rindx);
			else if (!ai_reml(y,P, Py, prev_varcmp, varcmp,n, rindx,step)) flag_EM = true;

		}

		constrain_varcmp(y,varcmp,n, rindx);

		lgL = -0.5 * (logdet_Xt_Vi_X + logdet + (y.transpose() * Py)(0, 0));


		//change step size
		if (iter > 300){
			L_size = L_history.size();
			for(int i=0;i<L_size-1;i++)
			{
				if(lgL - L_history[i]<1e-30 || lgL - L_history[i]< -1e-30){
					step = step*0.8;
					break;
				}
			}
			L_history.push_back(lgL);
		}


		// cout << lgL << endl;
		dlogL = lgL - prev_lgL;

		//converge
		if ((varcmp - prev_varcmp).squaredNorm() / varcmp.squaredNorm() < 1e-8 && (fabs(dlogL) < 1e-4 || (fabs(dlogL) < 1e-2 && dlogL < 0))) {
			prev_varcmp = varcmp;
			calcu_Vi(Vi, prev_varcmp, logdet,n, rindx);
			flag_converge = true;
			flag_inv_Vi = true;
			flag_inv_P = true;
			flag_not_itermax = true;
			return varcmp;
		}

		//max iter
		if(iter == maxiter)
		{
			calcu_Vi(Vi, prev_varcmp, logdet,n, rindx);
			cout << "Warning: Log-likelihood not converged. Results are not reliable.\nYou can specify the parameter max_iter to allow for more iterations." << endl;
			flag_not_itermax = false;

		}

		prev_varcmp = varcmp;
		prev_lgL = lgL;
	}

	return varcmp;
}
