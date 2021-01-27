#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_dq.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

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

	for (int i = 0; i < dF_dc.size(); ++i) {
		//std::cout << "i: " << i << std::endl;
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);

		Eigen::Matrix<DerivedC, 9, 9> d2psi_dF2;

		// TODO move this to separate file
		Eigen::Matrix<DerivedC, 9, 1> F_flat; 
		F_flat.setZero();
		for (int j = 0; j < W_I[i].size(); ++j) {
			int idx = W_I[i](j);
			F_flat += dF_dc[i].block(0,3*k*j, 9,3*k) * c.segment(3 * k*idx, 3 * k);
		}
		//std::cout << " did we get fflat? " << F_flat << std::endl;
		Eigen::Matrix3d F = unflatten<3, 3>(F_flat);

		//std::cout << "F: " << F << std::endl;
		d2psi_neohookean_dF2(d2psi_dF2, F, CD);
		//std::cout << "d2psi: " << i << ": \n" << d2psi_dF2 << std::endl;

		Eigen::MatrixXx<DerivedRet> tmp = dF_dc[i].transpose() * d2psi_dF2 * dF_dc[i] * volume(i);
		if (i == 0) {
			//std::cout << "tmp: " << i << ": \n" << tmp << std::endl;
			//std::cout << "dF_dc: " << i << ": \n" << dF_dc[i] << std::endl;

			//std::cout << "tmp [" << tmp.rows() << ", " << tmp.cols() << std::endl;
		}

		// Project each local stiffness matrix to PSD
		/*
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXx<DerivedRet>> es(tmp);
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    for (int i = 0; i < tmp.rows(); ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    tmp = Evec * DiagEval * Evec.transpose();
		*/
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

				//std::cout << "r: " << r << " c: " << c << std::endl;
				g.coeffRef(r, c) += tmp(ii,jj);
			}
		}
	}
}

#include <iostream>