#ifndef SIM_VEM3DMESH_NEOHOOKEAN_Q_H
#define SIM_VEM3DMESH_NEOHOOKEAN_Q_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <psi_neohookean_F.h> // TODO remove and replace with neohookean dq2


namespace sim {

template<typename DerivedC, typename DerivedVol, typename DerivedParam>
double vem3dmesh_neohookean_q(const Eigen::VectorXx<DerivedC> &c,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dc,
	const std::vector<Eigen::SparseMatrix<double>, Eigen::aligned_allocator<Eigen::SparseMatrix<double>> > & dF_dc_S,
	int d);
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/vem3dmesh_neohookean_q.cpp>
#endif

#endif