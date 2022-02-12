#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_q.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedC, typename DerivedVol, typename DerivedParam>
double sim::vem3dmesh_neohookean_q(const Eigen::VectorXx<DerivedC> &c,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dc,
	const std::vector<Eigen::SparseMatrix<double>, Eigen::aligned_allocator<Eigen::SparseMatrix<double>> > & dF_dc_S,
	int d) {

	double energy = 0;

	for (int i = 0; i < dF_dc.size(); ++i) {
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);

		Eigen::MatrixXd F_flat = dF_dc[i] * dF_dc_S[i] * c;

		Eigen::Matrix3d F = unflatten<3, 3>(F_flat);

		energy += psi_neohookean_F(F, params) *  volume(i);
	}
	
	return energy;
}

#include <iostream>