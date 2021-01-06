#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_dq2.h>
#endif

#ifdef vem_USE_OPENMP
#include <omp.h>
#endif

#include <chrono>

template<typename DerivedRet, typename DerivedA, typename DerivedDF, typename DerivedDM, typename DerivedW, typename DerivedVol, typename DerivedParam>
void sim::vem3dmesh_neohookean_dq2(Eigen::MatrixXx<DerivedRet> &g,
	const Eigen::MatrixXx<DerivedA> &A,
	const Eigen::MatrixXx<DerivedDF> &dF_dq,
	const Eigen::MatrixXx<DerivedDM> &dM_dq,
	const Eigen::MatrixXx<DerivedW> &w,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	int k, int n,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dq_sp,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & E) {

	g.resize(3 * n, 3 * n);
	g.setZero();

	double time_0 = 0;

	double time_1 = 0;
	double time_2 = 0;
	double time_3 = 0;
	using namespace std::chrono;

	//high_resolution_clock::time_point t_begin = high_resolution_clock::now();

	//#pragma omp parallel for
	for (int i = 0; i < w.rows(); ++i) {
		
		//high_resolution_clock::time_point t0 = high_resolution_clock::now();

		// TODO directly pass the block...
		Eigen::MatrixXx<DerivedDF> dFi;  //= dF_dq.row(i);
		Eigen::MatrixXx<DerivedDM> dMi_dq = dM_dq.row(i);
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);

		//high_resolution_clock::time_point t1 = high_resolution_clock::now();

		Eigen::MatrixXx<DerivedA> Ai;
		Ai.resize(3, k);
		Ai.setZero();
		for (int j = 0; j < A.rows(); ++j) {
			Eigen::MatrixXx<DerivedA> wAij = A.row(j).array() * w(i, j);
			Eigen::Map<const Eigen::MatrixXx<DerivedA>> Aij(wAij.data(), 3,k);
			Ai = Ai + Aij;
		}
		
		//high_resolution_clock::time_point t2 = high_resolution_clock::now();

		Eigen::MatrixXx<DerivedRet> tmp;
		vem3d_neohookean_dq2(tmp, Ai, dFi, dMi_dq, volume(i), CD, k, n, dF_dq_sp[i]);
		

		//high_resolution_clock::time_point t3 = high_resolution_clock::now();

		for (unsigned int jj = 0; jj < tmp.cols(); ++jj) {
			for (unsigned int ii = 0; ii < tmp.rows(); ++ii) {
		
				int r = E[i](ii);
				int c = E[i](jj);

				//std::cout << "i: " << i << " rows " << tmp.rows() << " cols: " << tmp.cols() << " r: " << r << " c: " << c << std::endl;
				g.coeffRef(r, c) += tmp(ii,jj);
			}
		}

		/*
		high_resolution_clock::time_point t4 = high_resolution_clock::now();

		duration<double> time_span0 = duration_cast<duration<double>>(t1 - t0);
		duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
		duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
		duration<double> time_span3 = duration_cast<duration<double>>(t4 - t3);
		time_0 += time_span0.count();

		time_1 += time_span1.count();
		time_2 += time_span2.count();
		time_3 += time_span3.count();
		*/
	}
	/*
	high_resolution_clock::time_point t_end = high_resolution_clock::now();
	duration<double> total_span = duration_cast<duration<double>>(t_begin - t_end);
	double total_time = total_span.count();
	std::cout << "Total: " << total_time << " Time access: " << time_0 << " Time blend: " << time_1 << " Time Stiffness: " << time_2 << " Time Assemble: " << time_3 << std::endl;
	*/
}

#include <iostream>