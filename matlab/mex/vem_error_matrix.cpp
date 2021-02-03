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

#include <parse_matlab_cells.h>

#include <vem_error_matrix.h>

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]) {
	/* variable declarations here */
	Eigen::MatrixXd g;

	Eigen::SparseMatrix<double> L;
	EigenMatrixList<double> Y;
	EigenVectorList<int> W_I, C_I;

	sim::parse_cell(prhs, 0, Y);
	igl::matlab::parse_rhs(prhs + 1, L);
	sim::parse_cell_index(prhs, 2, W_I);
	sim::parse_cell_index(prhs, 3, C_I);
	int d = (int)*mxGetPr(prhs[4]);
	int k = (int)*mxGetPr(prhs[5]);
	int n = (int)*mxGetPr(prhs[6]);

	sim::vem_error_matrix(g, Y, L, W_I, C_I, d, k, n);

	igl::matlab::prepare_lhs_double(g, plhs);
}