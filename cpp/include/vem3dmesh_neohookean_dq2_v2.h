#ifndef SIM_VEM3DMESH_NEOHOOKEAN_DQ2_V2_H
#define SIM_VEM3DMESH_NEOHOOKEAN_DQ2_V2_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <vem3d_neohookean_dq2.h>

namespace sim {

template<typename DerivedRet, typename DerivedA, typename DerivedDF, typename DerivedDM, typename DerivedW, typename DerivedVol, typename DerivedParam>
void vem3dmesh_neohookean_dq2_v2(Eigen::SparseMatrix<DerivedRet> &g,
	const Eigen::MatrixXx<DerivedA> &A,
	const Eigen::MatrixXx<DerivedDF> &dF_dq,
	const Eigen::MatrixXx<DerivedDM> &dM_dq,
	const Eigen::MatrixXx<DerivedW> &w,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	int k, int n,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dq_sp,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & E);
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/vem3dmesh_neohookean_dq2_v2.cpp>
#endif

#endif