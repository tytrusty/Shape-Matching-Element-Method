#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3d_neohookean_dq2.h>
#endif

#include <chrono>


template<typename GradientType, typename DerivedA, typename DerivedDF, typename DerivedDM, typename DerivedVol, typename ParamType>
void sim::vem3d_neohookean_dq2(Eigen::MatrixXx<GradientType> &g,
	const Eigen::MatrixXx<DerivedA> &A,
	const Eigen::MatrixXx<DerivedDF> &dF_dq,
	const Eigen::MatrixXx<DerivedDM> &dM_dq,
	DerivedVol volume,
	const Eigen::MatrixXx<ParamType> &params,
	int k, int n,
	const Eigen::MatrixXd &dF_dq_sp) {

	Eigen::Matrix<DerivedDF, 9, 9> d2psi_dF2;
	Eigen::Map<const Eigen::MatrixXx<DerivedDM>> dMi(dM_dq.data(), k, 3);


	//grab per element positions
	d2psi_neohookean_dF2(d2psi_dF2, unflatten<3, 3>((A * dMi).eval()), params);

	//auto start = std::chrono::high_resolution_clock::now();
	//Eigen::Map<const Eigen::MatrixXx<DerivedDF>> dFi(dF_dq.data(), 3*3, 3*n);
	//g = dFi.transpose() * d2psi_dF2 * dFi * volume;
	//auto stop = std::chrono::high_resolution_clock::now();
	//std::cout << "DENSE d2spi* dfi: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;
	//start = stop;
	g = dF_dq_sp.transpose() * d2psi_dF2 * dF_dq_sp * volume;
	//stop = std::chrono::high_resolution_clock::now();
	//std::cout << "SPARSE d2spi* dfi: " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;
}