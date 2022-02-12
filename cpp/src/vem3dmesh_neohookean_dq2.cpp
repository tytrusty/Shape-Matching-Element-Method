#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_dq.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

#include <chrono>

template<typename DerivedRet, typename DerivedC, typename DerivedVol, typename DerivedParam>
void sim::vem3dmesh_neohookean_dq2(Eigen::MatrixXx<DerivedRet> &g,
	const Eigen::VectorXx<DerivedC> &c,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dc,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & W_I,
	int k, int n, int m) {

	g.resize(3 * (k*n + m), 3 * (k*n + m));
	g.setZero();

	double time_0 = 0;
	double time_1 = 0;
	double time_2 = 0;
	double time_3 = 0;
	using namespace std::chrono;

	for (int i = 0; i < dF_dc.size(); ++i) {
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);

		Eigen::Matrix<DerivedC, 9, 9> d2psi_dF2;

		Eigen::Matrix<DerivedC, 9, 1> F_flat; 
		F_flat.setZero();
		for (int j = 0; j < W_I[i].size(); ++j) {
			int idx = W_I[i](j);
			F_flat += dF_dc[i].block(0,3*k*j, 9,3*k) * c.segment(3 * k*idx, 3 * k);
		}
		Eigen::Matrix3d F = unflatten<3, 3>(F_flat);

		d2psi_neohookean_dF2(d2psi_dF2, F, CD);

		Eigen::MatrixXx<DerivedRet> tmp = dF_dc[i].transpose() * d2psi_dF2 * dF_dc[i] * volume(i);

		// Assembly
		int kd = 3 * k;
		for (int jj = 0; jj < tmp.cols(); ++jj) {
			for (int ii = 0; ii < tmp.rows(); ++ii) {
				int r_i = W_I[i](ii / kd);
				int c_j = W_I[i](jj / kd);
				int r_offset = ii % kd;
				int c_offset = jj % kd;
				int r = r_i * kd + r_offset;
				int c = c_j * kd + c_offset;

				g.coeffRef(r, c) += tmp(ii,jj);
			}
		}

	}
}

#include <iostream>