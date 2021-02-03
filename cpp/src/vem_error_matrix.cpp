#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedRet, typename DerivedY>
void sim::vem_error_matrix(Eigen::MatrixXx<DerivedRet> &g,
	const EigenMatrixList<DerivedY>& Y,
	const Eigen::SparseMatrix<double> & L,
	const EigenVectorList<int>& W_I,
	const EigenVectorList<int>& C_I,
	int d, int k, int n) {

	std::cout << "ysize: " << Y.size() << std::endl;
	int sz = L.cols();

	Eigen::MatrixXd tmp1;
	Eigen::MatrixXd tmp2;
	tmp1.resize(L.rows(), L.rows());
	tmp2.resize(L.rows(), L.cols());
	tmp1.setZero();
	tmp2.setZero();
	
	int dk = d * k;
	int dkn = dk * n;

	auto assembly_idx = [=](int idx, int dkm, int i) {
		int ret;
		if (idx < dkm) {
			ret = W_I[i](idx / dk) * dk + idx % dk;
		}
		else {
			ret = dkn + C_I[i]((idx - dkm) / d) * d + idx % d;
		}
		return ret;
	};

	for (int i = 0; i < Y.size(); ++i) {
		Eigen::MatrixXd tmp = Y[i].transpose() * Y[i];

		int dkm = dk * W_I[i].rows(); // number of associated parts

		// Assembly
		for (int jj = 0; jj < tmp.cols(); ++jj) {
			for (int ii = 0; ii < tmp.rows(); ++ii) {
				int r = assembly_idx(ii, dkm, i);
				int c = assembly_idx(jj, dkm, i);
				tmp1.coeffRef(r, c) += tmp(ii, jj);
			}
		}

		for (int dim = 0; dim < d; ++dim) {
			// Copy weighted polynomial basis matrix
			for (int ii = 0; ii < W_I[i].rows(); ++ii) {
				for (int jj = 0; jj < k; ++jj) {
					int r = W_I[i](ii)*dk + dim * k + jj;
					int c = i * d + dim;
					// Same for each dim, so only copy from one row of Y
					tmp2.coeffRef(r, c) = Y[i](0, ii*dk + jj);
				}
			}

			// Copy weighted deformation origin entries
			for (int ii = 0; ii < C_I[i].rows(); ++ii) {
				int r = dkn + C_I[i](ii)*d + dim;
				int c = i * d + dim;
				// Same for each dim, so only copy from one row of Y
				tmp2.coeffRef(r, c) = Y[i](0, dkm + ii*d);
			}
		}
	}

	g = L.transpose() * tmp1 * L - L.transpose()*tmp2 - tmp2.transpose()*L 
		+ Eigen::MatrixXd::Identity(L.cols(),L.cols());
}