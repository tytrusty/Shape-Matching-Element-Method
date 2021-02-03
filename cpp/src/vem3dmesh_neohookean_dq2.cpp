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
		//std::cout << "i: " << i << std::endl;
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);

		Eigen::Matrix<DerivedC, 9, 9> d2psi_dF2;

		//high_resolution_clock::time_point t0 = high_resolution_clock::now();

		// TODO move this to separate file
		Eigen::Matrix<DerivedC, 9, 1> F_flat; 
		F_flat.setZero();
		for (int j = 0; j < W_I[i].size(); ++j) {
			int idx = W_I[i](j);
			F_flat += dF_dc[i].block(0,3*k*j, 9,3*k) * c.segment(3 * k*idx, 3 * k);
		}
		Eigen::Matrix3d F = unflatten<3, 3>(F_flat);

		//high_resolution_clock::time_point t1 = high_resolution_clock::now();

		//std::cout << "F: " << F << std::endl;
		d2psi_neohookean_dF2(d2psi_dF2, F, CD);
		//std::cout << "d2psi: " << i << ": \n" << d2psi_dF2 << std::endl;

		//high_resolution_clock::time_point t2 = high_resolution_clock::now();

		Eigen::MatrixXx<DerivedRet> tmp = dF_dc[i].transpose() * d2psi_dF2 * dF_dc[i] * volume(i);

		//high_resolution_clock::time_point t3 = high_resolution_clock::now();

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

		//high_resolution_clock::time_point t4 = high_resolution_clock::now();
		//duration<double> time_span0 = duration_cast<duration<double>>(t1 - t0);
		//duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
		//duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
		//duration<double> time_span3 = duration_cast<duration<double>>(t4 - t3);
		//time_0 += time_span0.count();
		//time_1 += time_span1.count();
		//time_2 += time_span2.count();
		//time_3 += time_span3.count();
	}
	//std::cout << " Time F: " << time_0 << " Time Stiffness: " << time_1 << " Time Stiffness mult: " << time_2 << " Time Assemble: " << time_3 << std::endl;

}

#include <iostream>