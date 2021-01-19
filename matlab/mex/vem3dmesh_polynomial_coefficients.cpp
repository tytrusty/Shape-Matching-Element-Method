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


//bartels
#include <vem3dmesh_polynomial_coefficients.h>

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
    /* variable declarations here */
    Eigen::VectorXd c;
    Eigen::MatrixXd x, L;

  	std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > E;

		igl::matlab::parse_rhs_double(prhs+0, x);
    igl::matlab::parse_rhs_double(prhs+1, L);
		parse_cell_index(prhs, 2, E);

		sim::vem3dmesh_polynomial_coefficients(c, x, L, E);

		igl::matlab::prepare_lhs_double(c, plhs);
}