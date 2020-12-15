#ifndef SIM_VEM3D_NEOHOOKEAN_DQ2_H
#define SIM_VEM3D_NEOHOOKEAN_DQ2_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <d2psi_neohookean_dF2.h>

namespace sim {

template<typename GradientType, typename DerivedA, typename DerivedDF, typename DerivedDM, typename DerivedVol, typename ParamType>
void vem3d_neohookean_dq2(Eigen::MatrixXx<GradientType> &g,
	const Eigen::MatrixXx<DerivedA> &A,
	const Eigen::MatrixXx<DerivedDF> &dF_dq,
	const Eigen::MatrixXx<DerivedDM> &dM_dq, 
	DerivedVol volume,
	const Eigen::MatrixXx<ParamType> &params,
	int k, int n,
	const Eigen::MatrixXd &dF_dq_sp);

}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/vem3d_neohookean_dq2.cpp>
#endif

#endif 