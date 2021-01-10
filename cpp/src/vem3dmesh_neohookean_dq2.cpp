#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_dq2.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedRet, typename DerivedC, typename DerivedDM, typename DerivedVol, typename DerivedParam>
void sim::vem3dmesh_neohookean_dq2(Eigen::MatrixXx<DerivedRet> &g,
	const Eigen::VectorXx<DerivedC> &c,
	const Eigen::MatrixXx<DerivedDM> &dM_dX,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dc,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & W,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & W_S,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & W_I,
	int k, int n) {

	g.resize(3 * (k * n + 1), 3 * (k * n + 1));
	g.setZero();

	for (int i = 0; i < dM_dX.rows(); ++i) {
		
		Eigen::MatrixXx<DerivedDM> dMi_dX = dM_dX.row(i);
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);

		Eigen::Map<const Eigen::MatrixXx<DerivedDM>> dMi_map(dMi_dX.data(), 9, 3*k);
		Eigen::Matrix<DerivedC, 9, 9> d2psi_dF2;

		// TODO move this to separate file
		Eigen::Matrix3d F = unflatten<3, 3>((dMi_map * W[i] * W_S[i] * c).eval());


		d2psi_neohookean_dF2(d2psi_dF2, F, CD);

		Eigen::MatrixXx<DerivedRet> tmp = dF_dc[i].transpose() * d2psi_dF2 * dF_dc[i] * volume(i);
		//Eigen::MatrixXx<DerivedRet> tmp2 = W_S[i].transpose() *  (dF_dc[i].transpose() * d2psi_dF2 * dF_dc[i]) * W_S[i] * volume(i);
		//g += tmp2;

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