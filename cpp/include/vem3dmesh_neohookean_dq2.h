#ifndef SIM_VEM3DMESH_NEOHOOKEAN_DQ2_H
#define SIM_VEM3DMESH_NEOHOOKEAN_DQ2_H

#include <Eigen/Dense>
#include <EigenTypes.h>

//#include <vem3d_neohookean_dq2.h>
#include <d2psi_neohookean_dF2.h> // TODO remove and replace with neohookean dq2


namespace sim {

template<typename DerivedRet, typename DerivedC, typename DerivedDM, typename DerivedVol, typename DerivedParam>
void vem3dmesh_neohookean_dq2(Eigen::MatrixXx<DerivedRet> &g,
	const Eigen::VectorXx<DerivedC> &c,
	const Eigen::MatrixXx<DerivedDM> &dM_dX,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dc,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & W,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & W_S,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & W_I,
	int k, int n);
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/vem3dmesh_neohookean_dq2.cpp>
#endif

#endif