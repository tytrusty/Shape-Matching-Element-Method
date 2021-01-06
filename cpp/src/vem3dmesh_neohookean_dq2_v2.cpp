#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_dq2_v2.h>
#endif

#ifdef vem_USE_OPENMP
#include <omp.h>
#endif

#include <chrono>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>

template<typename DerivedRet, typename DerivedA, typename DerivedDF, typename DerivedDM, typename DerivedW, typename DerivedVol, typename DerivedParam>
void sim::vem3dmesh_neohookean_dq2_v2(Eigen::SparseMatrix<DerivedRet> &g,
	const Eigen::MatrixXx<DerivedA> &A,
	const Eigen::MatrixXx<DerivedDF> &dF_dq,
	const Eigen::MatrixXx<DerivedDM> &dM_dq,
	const Eigen::MatrixXx<DerivedW> &w,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	int k, int n,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dq_sp,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & E) {

	g.reserve(235368);
	Eigen::MatrixXd g2;
	g2.resize(3 * n, 3 * n);
	g2.setZero();
	std::vector< Eigen::Triplet<double> > tripletList;
	tripletList.reserve(2520000);

	double time_0 = 0;
	double time_1 = 0;
	double time_2 = 0;
	double time_3 = 0;
	using namespace std::chrono;

	high_resolution_clock::time_point t_begin = high_resolution_clock::now();

	std::cout << "num threads from eigen!: " << Eigen::nbThreads() << std::endl;
    
    //#pragma omp parallel for
	for (int i = 0; i < w.rows(); ++i) {
		
		high_resolution_clock::time_point t0 = high_resolution_clock::now();

		// TODO directly pass the block...
		Eigen::MatrixXx<DerivedDF> dFi;  //= dF_dq.row(i);
		Eigen::MatrixXx<DerivedDM> dMi_dq = dM_dq.row(i);
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);

		high_resolution_clock::time_point t1 = high_resolution_clock::now();

		Eigen::MatrixXx<DerivedA> Ai;
		Ai.resize(3, k);
		Ai.setZero();
		for (int j = 0; j < A.rows(); ++j) {
			Eigen::MatrixXx<DerivedA> wAij = A.row(j).array() * w(i, j);
			Eigen::Map<const Eigen::MatrixXx<DerivedA>> Aij(wAij.data(), 3,k);
			Ai = Ai + Aij;
		}
		
		high_resolution_clock::time_point t2 = high_resolution_clock::now();

		Eigen::MatrixXx<DerivedRet> tmp;
		vem3d_neohookean_dq2(tmp, Ai, dFi, dMi_dq, volume(i), CD, k, n, dF_dq_sp[i]);

		high_resolution_clock::time_point t3 = high_resolution_clock::now();

		int offset = 120 * 120 * i;
		typedef Eigen::Triplet<double> T;
		for (unsigned int jj = 0; jj < tmp.cols(); ++jj) {
			for (unsigned int ii = 0; ii < tmp.rows(); ++ii) {
		
				int r = E[i](ii);
				int c = E[i](jj);
				tripletList.push_back(T(r, c, tmp(ii, jj)));
			}
		}
		
		high_resolution_clock::time_point t4 = high_resolution_clock::now();

		duration<double> time_span0 = duration_cast<duration<double>>(t1 - t0);
		duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
		duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
		duration<double> time_span3 = duration_cast<duration<double>>(t4 - t3);
		time_0 += time_span0.count();

		time_1 += time_span1.count();
		time_2 += time_span2.count();
		time_3 += time_span3.count();
	}

	std::sort(tripletList.begin(), tripletList.end(), [=](const auto& lhs, const auto& rhs)
	{
		return (3 * n*lhs.row() + lhs.col()) < (3 * n*rhs.row() + rhs.col());
	});
	
	std::cout << "triplet size: " << tripletList.size() << std::endl;

	std::vector< Eigen::Triplet<double> > mainTriplets;
	std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> bins;
	int start = 0;
	for (int i = 1; i < tripletList.size(); ++i) {
		const Eigen::Triplet<double> trip = tripletList[i];
		if ( (i == tripletList.size() - 1) || (trip.col() != tripletList[start].col() || trip.row() != tripletList[start].row())) {
			mainTriplets.push_back(tripletList[start]);
			int new_size = i - start;
			Eigen::VectorXd tmp(new_size); /// excludes the final i values
			if (i == tripletList.size() - 1) new_size += 1;
			
			for (int j = start; j < start + new_size; ++j) {
				tmp(j - start) = tripletList[j].value();
			}
			bins.emplace_back(tmp);
			start = i;
		}
	}
	std::cout << "bins size:! " << bins.size() << std::endl;
	high_resolution_clock::time_point sortb = high_resolution_clock::now();
	# pragma omp parallel for
	for (int i = 0; i < mainTriplets.size(); ++i) {
		double sum = bins[i].sum();
		//double sum = 0;
		//for (int j = 0; j < bins[i].size(); ++j) {
		//	sum += bins[i](j);
		//}
		g2.coeffRef(mainTriplets[i].row(), mainTriplets[i].col()) = sum;
	}


	//thrust::sort_by_key(keys.begin(), keys.end(), values.begin());
	//int C[2520000]; // input keys
	//double D[2520000];
	//thrust::pair<int*, double*> new_end;
	//new_end = thrust::reduce_by_key(keys.begin(), keys.end(), values.begin(), C, D);
	//std::cout << ""
	high_resolution_clock::time_point sorte = high_resolution_clock::now();
	std::cout << "sort time: " << duration_cast<duration<double>>(sorte - sortb).count() << std::endl;
	//for (int i = 0; i < 10; ++i) {
	//	std::cout << "[" << i << " " << tripletList[i].row() << "," << tripletList[i].col() << std::endl;
	//}
	//g.setFromTriplets(tripletList.begin(), tripletList.end());
	high_resolution_clock::time_point t_end = high_resolution_clock::now();
	duration<double> total_span = duration_cast<duration<double>>(t_begin - t_end);
	double total_time = total_span.count();
	std::cout << "Total: " << total_time << " Time access: " << time_0  << " Time blend: " << time_1 << " Time Stiffness: " << time_2 << " Time Assemble: " << time_3 << std::endl;
	
}

#include <iostream>