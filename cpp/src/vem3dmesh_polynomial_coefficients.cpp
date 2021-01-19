#ifdef SIM_STATIC_LIBRARY
# include<../include/vem3dmesh_polynomial_coefficients.h>
#endif

#ifdef VEM_USE_OPENMP
#include <omp.h>
#endif

template<typename DerivedRet, typename DerivedX, typename DerivedL>
void sim::vem3dmesh_polynomial_coefficients(Eigen::VectorXx<DerivedRet> &c,
	const Eigen::MatrixXx<DerivedX> &x,
	const Eigen::MatrixXx<DerivedL> &L,
	const std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > & E) {

		Eigen::MatrixXd B(x.rows(), x.cols());
		B.setZero();
		int col_start = 0;
		Eigen::MatrixXd x_slice;
		for (int i = 0; i < E.size(); i++) {
			igl::slice(x, E[i], 2, x_slice);
			B.block(0, col_start, 3, E[i].rows()) = x_slice;

			col_start += E[i].rows();

			// if (i == 0) {
			// 	std::cout << "E[i] = " << E[i].transpose() << std::endl;			
			// 	std::cout << "x_slice = " << x_slice << std::endl;
			// 	std::cout << "B = " << B << std::endl;
			// }
		}
		Eigen::VectorXd b = Eigen::Map<const Eigen::VectorXd>(B.data(), B.size());

		c = L * b;
}

#include <iostream>