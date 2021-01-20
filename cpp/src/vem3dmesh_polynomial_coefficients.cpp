#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_polynomial_coefficients.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedRet, typename DerivedX>
void sim::vem3dmesh_polynomial_coefficients(Eigen::VectorXx<DerivedRet> &c,
	const Eigen::MatrixXx<DerivedX> &x,
	const Eigen::SparseMatrix<double> &L,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & E){

		// std::cout << "start polynomial_coefficients" << std::endl;
		Eigen::MatrixXd B(x.rows(), x.cols());
		B.setZero();
		int col_start = 0;
		Eigen::MatrixXd x_slice;
		// std::cout << "1 polynomial_coefficients" << std::endl;
		for (int i = 0; i < E.size(); i++) {
			igl::slice(x, E[i], 2, x_slice);
			// std::cout << "2 polynomial_coefficients" << std::endl;
			B.block(0, col_start, 3, E[i].rows()) = x_slice;
			// std::cout << "3 polynomial_coefficients" << std::endl;

			col_start += E[i].rows();
			// std::cout << "4 polynomial_coefficients" << std::endl;
		}
		// std::cout << "5 polynomial_coefficients" << std::endl;
		Eigen::VectorXd b = Eigen::Map<const Eigen::VectorXd>(B.data(), B.size());
		// std::cout << "6 polynomial_coefficients" << std::endl;
		
		// Solve for polynomial coefficients (projection operators).
		c = L * b;
		// std::cout << "pass polynomial_coefficients" << std::endl;
}

#include <iostream>