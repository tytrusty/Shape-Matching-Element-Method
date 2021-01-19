#ifndef SIM_VEM3DMESH_SIMULATE_ONE_STEP_H
#define SIM_VEM3DMESH_SIMULATE_ONE_STEP_H

#include <Eigen/Dense>
#include <EigenTypes.h>

// VEM energy, gradient and hessian
#include <vem3dmesh_neohookean_dq2.h>
#include <vem3dmesh_neohookean_q.h>
#include <vem3dmesh_neohookean_dq.h>
#include <vem3dmesh_polynomial_coefficients.h>

// Bartels
#include <optimization_gradient_based.h>

namespace sim {

template<typename DerivedRet, typename DerivedQ, typename DerivedF, typename DerivedVol, typename DerivedParam>
void vem3dmesh_simulate_one_step(Eigen::VectorXx<DerivedRet> &qdot_new,
	const Eigen::VectorXx<DerivedQ> &q,
	const Eigen::VectorXx<DerivedQ> &qdot,
	const Eigen::VectorXx<DerivedF> &f_ext,
	const Eigen::MatrixXd &x_fixed,
	const Eigen::VectorXx<DerivedVol> &volume,
	const Eigen::MatrixXx<DerivedParam> &params,
	const std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & dF_dc,
	const std::vector<Eigen::SparseMatrix<double>, Eigen::aligned_allocator<Eigen::SparseMatrix<double>> > & dF_dc_S,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & W_I,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & E,
	const Eigen::MatrixXd &M,
	const Eigen::MatrixXd &ME,
	const Eigen::MatrixXd &L,
	const Eigen::SparseMatrix<double> &P,
	const Eigen::SparseMatrix<double> &J,
	int k, int n, int d, int x0_coms_size, double k_stability, double dt);
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/vem3dmesh_simulate_one_step.cpp>
#endif

#endif