#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>

enum INVmethod {
	INV_LLT = 1,
	INV_LU = 2,
	INV_QR = 3,
	INV_FQR = 4,
	INV_SVD = 5,
	INV_ERR = 100
};



template<typename MatrixType>
bool SquareMatrixInverse(MatrixType &A, double &logdet, int &rank, INVmethod &method) {
	int n = A.rows();
	bool ret = false;
	switch (method) {
	case INV_LLT:
	{
		Eigen::LLT<MatrixType> llt(A);
		if (!llt.info()) {
			logdet = llt.matrixLLT().diagonal().array().abs().log().sum() * 2;
			A = llt.solve(MatrixType::Identity(n, n));
			method = INV_LLT;
			ret = true;
			break;
		}

		/*if (_LLT(A, logdet)) {
			method = INV_LLT;
			ret = true;
			break;
		}*/
	}
	case INV_LU:
	{
		Eigen::PartialPivLU<MatrixType> lu(A);
		double det = lu.determinant();
		//std::cout << "LU det: " << std::scientific << det << std::endl;
		if (det >= 1e-10 || det <= -1e-10) {
			logdet = lu.matrixLU().diagonal().array().abs().log().sum();
			A = lu.inverse();
			method = INV_LU;
			ret = true;
			break;
		}
	}
	case INV_QR:
	{
		Eigen::HouseholderQR<MatrixType> qr(A);
		double det = qr.absDeterminant();
		//std::cout << "QR det: " << std::scientific << det << std::endl;
		if (det > 1e-16) {
			logdet = qr.logAbsDeterminant();
			A = qr.solve(MatrixType::Identity(n, n));
			method = INV_QR;
			ret = true;
			break;
		}
	}
	case INV_FQR:
		// not necessary 
		// Eigen::HouseholderQR<MatrixType> qr(A);
	{
		Eigen::ColPivHouseholderQR<MatrixType> qr(A);
		if (qr.isInvertible()) {
			logdet = qr.logAbsDeterminant();
			A = qr.inverse();
			method = INV_QR;
			ret = true;
			break;
		}
		else {
			rank = qr.rank();
			// it will be extreme slow or not accurate
			if (n > 50000 || 1.0 * rank / n < 0.99) {
				ret = false;
				break;
			}
		}
	}
	case INV_SVD:
		//Eigen::BDCSVD<MatrixType> svd(A, Eigen::ComputeThinU|Eigen::ComputeThinV);
		;
	default:
		rank = 0;
		ret = false;
		method = INV_ERR;
	}
	return ret;
}



