#ifndef SIM_VEM3DMESH_POLYNOMIAL_COEFFICIENTS_H
#define SIM_VEM3DMESH_POLYNOMIAL_COEFFICIENTS_H

#include <Eigen/Dense>
#include <EigenTypes.h>
#include <igl/slice.h>

namespace sim {

template<typename DerivedRet, typename DerivedX>
void vem3dmesh_polynomial_coefficients(Eigen::VectorXx<DerivedRet> &c,
	const Eigen::MatrixXx<DerivedX> &x,
	const Eigen::SparseMatrix<double> &L,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & E);
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/vem3dmesh_polynomial_coefficients.cpp>
#endif

#endif