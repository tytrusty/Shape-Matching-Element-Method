#ifndef SIM_VEM3DMESH_NEOHOOKEAN_DQ_H
#define SIM_VEM3DMESH_NEOHOOKEAN_DQ_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <dpsi_neohookean_dF.h> // TODO remove and replace with neohookean dq2


namespace sim {

template<typename DerivedRet, typename DerivedX, typename DerivedC, typename DerivedVol, typename DerivedParam>
void vem3dmesh_neohookean_dq(Eigen::VectorXx<DerivedRet> &g,
	const Eigen::MatrixXx<DerivedX> &x,
	const Eigen::VectorXx<DerivedC> &c,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dc,
	const std::vector<Eigen::SparseMatrix<double>, Eigen::aligned_allocator<Eigen::SparseMatrix<double>> > & dF_dc_S,
	const Eigen::MatrixXd ME,
	const Eigen::MatrixXd L,
	int k, int n, int d, int x0_coms_size, double k_stability);
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/vem3dmesh_neohookean_dq.cpp>
#endif

#endif