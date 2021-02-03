#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedRet, typename DerivedY, typename DerivedM>
void sim::vem_mass_matrix(Eigen::MatrixXx<DerivedRet> &g,
	const EigenMatrixList<DerivedY>& Y,
	const Eigen::SparseMatrix<double> & L,
	const Eigen::VectorXx<DerivedM>& mass,
	const EigenVectorList<int>& W_I,
	const EigenVectorList<int>& C_I,
	int d, int k, int n) {

	int sz = L.rows();
	g.resize(sz,sz);
	g.setZero();

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
		Eigen::MatrixXd tmp = Y[i].transpose() * Y[i] * mass(i);

		int dkm = dk * W_I[i].rows(); // number of associated parts

		// Assembly
		for (int jj = 0; jj < tmp.cols(); ++jj) {
			for (int ii = 0; ii < tmp.rows(); ++ii) {
				int r = assembly_idx(ii, dkm, i);
				int c = assembly_idx(jj, dkm, i);
				g.coeffRef(r, c) += tmp(ii, jj);
			}
		}
	}
	g = L.transpose() * g * L;
}