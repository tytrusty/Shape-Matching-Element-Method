#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_neohookean_dq.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedRet, typename DerivedExt, typename DerivedM, typename DerivedY>
void sim::vem_ext_force(Eigen::VectorXx<DerivedRet> &g,
	const Eigen::VectorXx<DerivedExt>& ext,
	const Eigen::VectorXx<DerivedM>& mass,
	const EigenMatrixList<DerivedY>& Y,
	const EigenSparseMatrixList<DerivedY>& Y_S) {

	int sz = Y_S[0].cols();
	g.resize(sz);
	g.setZero();

	for (int i = 0; i < Y.size(); ++i) {
		g += Y_S[i].transpose() * (Y[i].transpose() * ext * mass(i));
	}
}