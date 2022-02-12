#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_simulate_one_step.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedRet, typename DerivedQ, typename DerivedF, typename DerivedVol, typename DerivedParam>
void sim::vem3dmesh_simulate_one_step(Eigen::VectorXx<DerivedRet> &qdot_new,
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
	const Eigen::SparseMatrix<double> &L,
	const Eigen::SparseMatrix<double> &P,
	const Eigen::SparseMatrix<double> &J,
	int k, int n, int d, int x0_coms_size, double k_stability, double dt) {

	Eigen::MatrixXd JtPMPtJ = J.transpose() * P * M * P.transpose() * J;

	auto get_x_new_and_c = [&](Eigen::MatrixXd &x_new, Eigen::VectorXd &c, const Eigen::VectorXd &qdot_new){

		Eigen::VectorXd q_new = q + dt * qdot_new;
		Eigen::VectorXd tmp = P.transpose() * J * q_new;
		// compute new global position
		x_new = Eigen::Map<Eigen::MatrixXd>(tmp.data(), 3, x_fixed.cols()) + x_fixed;
		// solve for polynomial coefficients
		vem3dmesh_polynomial_coefficients(c, x_new, L, E);

	};

	auto value = [&](Eigen::VectorXd qdot_new) -> double{

		// compute new global position
		Eigen::MatrixXd x_new;
		Eigen::VectorXd c;
		get_x_new_and_c(x_new, c, qdot_new);
		Eigen::VectorXd x_new_vec = Eigen::Map<Eigen::VectorXd>(x_new.data(), x_new.size());

		double neohookean_energy = vem3dmesh_neohookean_q(c, volume, params, dF_dc, dF_dc_S, d);

		Eigen::MatrixXd tmp = -qdot_new.transpose() * JtPMPtJ * qdot;

		double energy = 0.5 * qdot_new.transpose() * JtPMPtJ * qdot_new
										+ tmp(0,0)
										+ k_stability * x_new_vec.transpose() * ME * x_new_vec 
										+ 0.5 * neohookean_energy
										- qdot_new.transpose()*J.transpose()*f_ext;
		
		return energy;
  };
	auto gradient = [&](Eigen::VectorXd &g, Eigen::VectorXd &qdot_new) {
		Eigen::MatrixXd x_new;
		Eigen::VectorXd c;
		get_x_new_and_c(x_new, c, qdot_new);

		// compute vem neohookean gradient
		Eigen::VectorXd g_neohookean;
		vem3dmesh_neohookean_dq(g_neohookean, x_new, c, volume, params, dF_dc, dF_dc_S, ME, L, 
														k, n, d, x0_coms_size, k_stability);
		g = JtPMPtJ * (qdot_new - qdot) + J.transpose() * (dt*P*g_neohookean - f_ext); 	

	};
	auto hessian = [&](Eigen::MatrixXd &H, Eigen::VectorXd &qdot_new) {
		Eigen::MatrixXd x_new;
		Eigen::VectorXd c;
		get_x_new_and_c(x_new, c, qdot_new);

		Eigen::MatrixXd K;
		vem3dmesh_neohookean_dq2(K, c, volume, params, dF_dc, W_I, k, n, x0_coms_size);
		K = L.transpose() * (K) * L;

		H = JtPMPtJ + J.transpose() * P * K * P.transpose() * J * dt * dt;
	};

	Eigen::VectorXd tmp_g(qdot.rows());
	Eigen::MatrixXd tmp_H(qdot.rows(), qdot.rows());
	Eigen::VectorXd tmp_d(qdot.rows());

	qdot_new = qdot;
	newtons_method_bisection(qdot_new, value, gradient, hessian, tmp_g, tmp_H, tmp_d, 1e-3, 1, 5);
}

#include <iostream>