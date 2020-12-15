#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_dq2.h>
#endif

#ifdef BARTELS_USE_OPENMP
#include <omp.h>
#endif

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

	//std::cout << "num threads from eigen: " << Eigen::nbThreads() << std::endl;
	//#pragma omp parallel for
	for (int i = 0; i < w.rows(); ++i) {
		
		// TODO directly pass the block...
		Eigen::MatrixXx<DerivedDF> dFi = dF_dq.row(i);
		Eigen::MatrixXx<DerivedDM> dMi_dq = dM_dq.row(i);
		Eigen::MatrixXx<DerivedParam> CD = params.row(i);
		
		Eigen::MatrixXx<DerivedA> Ai;
		Ai.resize(3, k);
		Ai.setZero();
		for (int j = 0; j < A.rows(); ++j) {
			Eigen::MatrixXx<DerivedA> wAij = A.row(j).array() * w(i, j);
			Eigen::Map<const Eigen::MatrixXx<DerivedA>> Aij(wAij.data(), 3,k);
			Ai = Ai + Aij;
		}
		
		Eigen::MatrixXx<DerivedRet> tmp;
		vem3d_neohookean_dq2(tmp, Ai, dFi, dMi_dq, volume(i), CD, k, n, dF_dq_sp[i]);
		
        //#pragma omp critical
		//g += tmp;


		for (unsigned int jj = 0; jj < tmp.cols(); ++jj) {
			for (unsigned int ii = 0; ii < tmp.rows(); ++ii) {
				int r = E[i](ii);
				int c = E[i](jj);

				//std::cout << "i: " << i << " rows " << tmp.rows() << " cols: " << tmp.cols() << " r: " << r << " c: " << c << std::endl;
				g.coeffRef(r, c) += tmp(ii,jj);
			}
		}
	}

}

#include <iostream>