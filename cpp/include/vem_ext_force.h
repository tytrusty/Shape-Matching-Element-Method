#ifndef SIM_VEM_EXT_FORCE_H
#define SIM_VEM_EXT_FORCE_H

#include <Eigen/Dense>
#include <EigenTypes.h>

namespace sim {

template<typename DerivedRet, typename DerivedExt, typename DerivedM, typename DerivedY>
void vem_ext_force(Eigen::VectorXx<DerivedRet> &g,
	const Eigen::VectorXx<DerivedExt>& ext,
	const Eigen::VectorXx<DerivedM>& mass,
	const EigenMatrixList<DerivedY>& Y,
	const EigenSparseMatrixList<DerivedY>& Y_S);
}

#include <../src/vem_ext_force.cpp>

#endif // SIM_VEM_EXT_FORCE_H