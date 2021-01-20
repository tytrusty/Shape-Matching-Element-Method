#include <mex.h>

//matlab hooks from libigl (taken directly from gptoolbox)
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/list_to_matrix.h>

#include <iostream>


//bartels
#include <vem3dmesh_simulate_one_step.h>

template <typename T>
using EigenMatrixList = std::vector<Eigen::MatrixXx<T>, Eigen::aligned_allocator<Eigen::MatrixXx<T>>>;

template <typename T>
using EigenVectorList = std::vector<Eigen::VectorXx<T>, Eigen::aligned_allocator<Eigen::VectorXx<T>>>;

template <typename T>
using EigenSparseMatrixList = std::vector<Eigen::SparseMatrix<T>, Eigen::aligned_allocator<Eigen::SparseMatrix<T>>>;


template <typename DerivedV>
void parse_cell(const mxArray *prhs[], int index, EigenMatrixList<DerivedV>& V) {
	const mxArray *cell = prhs[index];
	const mwSize *dims = mxGetDimensions(prhs[index]);
	const mxArray *cellElement;
	mwIndex jcell;
	for (jcell = 0; jcell < dims[0]; jcell++) {
		cellElement = mxGetCell(cell, jcell);
		Eigen::MatrixXx<DerivedV> element;
		igl::matlab::parse_rhs_double(&cellElement, element);
		
		V.emplace_back(element);
	}
}

template <typename DerivedV>
void parse_cell_index(const mxArray *prhs[], int index, EigenVectorList<DerivedV>& V) {
	const mxArray *cell = prhs[index];
	const mwSize *dims = mxGetDimensions(prhs[index]);
	const mxArray *cellElement;
	mwIndex jcell;
	for (jcell = 0; jcell < dims[0]; jcell++) {
		cellElement = mxGetCell(cell, jcell);
		Eigen::VectorXx<DerivedV> element;
		igl::matlab::parse_rhs_index(&cellElement, element);
		V.emplace_back(element);
	}
}

template <typename DerivedV>
void parse_cell_sparse(const mxArray *prhs[], int index, EigenSparseMatrixList<DerivedV>& V) {
	const mxArray *cell = prhs[index];
	const mwSize *dims = mxGetDimensions(prhs[index]);
	const mxArray *cellElement;
	mwIndex jcell;
	for (jcell = 0; jcell < dims[0]; jcell++) {
		cellElement = mxGetCell(cell, jcell);
		Eigen::SparseMatrix<DerivedV> element;
		igl::matlab::parse_rhs(&cellElement, element);
		V.emplace_back(element);
	}
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

		// std::cout << "mex start" << std::endl;
    /* variable declarations here */
    Eigen::VectorXd qdot_new;

    Eigen::MatrixXd x_fixed, params, M, ME;
    Eigen::VectorXd q, qdot, f_ext, volume, volumes; 

		Eigen::SparseMatrix<double> P, J, L;

		std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > dF_dc;
		std::vector<Eigen::SparseMatrix<double>, Eigen::aligned_allocator<Eigen::SparseMatrix<double>>> dF_dc_S;
		std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > W_I;
		std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > E;

		// std::cout << "1" << std::endl;

		igl::matlab::parse_rhs_double(prhs+0, q);
		// std::cout << "1.1" << std::endl;
    igl::matlab::parse_rhs_double(prhs+1, qdot);
		// std::cout << "1.2" << std::endl;
		igl::matlab::parse_rhs_double(prhs+2, f_ext);
		// std::cout << "1.3" << std::endl;
    igl::matlab::parse_rhs_double(prhs+3, x_fixed);
		// std::cout << "1.4" << std::endl;
    igl::matlab::parse_rhs_double(prhs+4, volumes);
		// std::cout << "1.5" << std::endl;
    igl::matlab::parse_rhs_double(prhs+5, params);
		// std::cout << "2" << std::endl;
		parse_cell(prhs, 6, dF_dc);
		parse_cell_sparse(prhs, 7, dF_dc_S);
		parse_cell_index(prhs, 8, W_I);
		parse_cell_index(prhs, 9, E);
		// std::cout << "3" << std::endl;
		igl::matlab::parse_rhs_double(prhs+10, M);
		igl::matlab::parse_rhs_double(prhs+11, ME);
		igl::matlab::parse_rhs(prhs+12, L);
		igl::matlab::parse_rhs(prhs+13, P);
		igl::matlab::parse_rhs(prhs+14, J);
		// std::cout << "4" << std::endl;

		int k = (int)*mxGetPr(prhs[15]);
		int n = (int)*mxGetPr(prhs[16]);
		int d = (int)*mxGetPr(prhs[17]);
		int x0_coms_size = (int)*mxGetPr(prhs[18]);
		double k_stability = (double)*mxGetPr(prhs[19]);
		double dt = (double)*mxGetPr(prhs[20]);
		// std::cout << "5" << std::endl;

		sim::vem3dmesh_simulate_one_step(qdot_new, q, qdot, f_ext, x_fixed, volumes, params, 
																			dF_dc, dF_dc_S, W_I, E, 
																			M, ME, L, P, J,
																			k, n, d, x0_coms_size, k_stability, dt);

		igl::matlab::prepare_lhs_double(qdot_new, plhs);
}