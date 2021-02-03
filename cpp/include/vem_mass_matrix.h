#ifndef SIM_VEM_MASS_MATRIX_H
#define SIM_VEM_MASS_MATRIX_H

#include <Eigen/Dense>
#include <EigenTypes.h>

namespace sim {

template<typename DerivedRet, typename DerivedY, typename DerivedM>
void vem_mass_matrix(Eigen::MatrixXx<DerivedRet> &g,
	const EigenMatrixList<DerivedY>& Y,
	const Eigen::SparseMatrix<double> & L,
	const Eigen::VectorXx<DerivedM>& mass,
	const EigenVectorList<int>& W_I,
	const EigenVectorList<int>& C_I,
	int d, int k, int n);
}

#include <../src/vem_mass_matrix.cpp>

#endif // SIM_VEM_MASS_MATRIX_H