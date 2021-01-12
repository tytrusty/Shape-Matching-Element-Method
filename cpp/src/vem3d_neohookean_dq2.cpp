#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3d_neohookean_dq2.h>
#endif

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

	d2psi_neohookean_dF2(d2psi_dF2, unflatten<3, 3>((A * dMi).eval()), params);
	g = dF_dq_sp.transpose() * d2psi_dF2 * dF_dq_sp * volume;

}