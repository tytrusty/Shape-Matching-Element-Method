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
#include <vem3dmesh_neohookean_dq2.h>
#include <vem3dmesh_neohookean_dq2_v2.h>

#include <thrust/detail/config.h>

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    /* variable declarations here */
    Eigen::MatrixXd g;

    Eigen::MatrixXd A, dF_dq, dM_dq, w, params;
    Eigen::VectorXd volumes; 

    igl::matlab::parse_rhs_double(prhs+0, A);
    igl::matlab::parse_rhs_double(prhs+1, dF_dq);
    igl::matlab::parse_rhs_double(prhs+2, dM_dq);
    igl::matlab::parse_rhs_double(prhs+3, w);
    igl::matlab::parse_rhs_double(prhs+4, volumes);
    igl::matlab::parse_rhs_double(prhs+5, params);
	int k = (int)*mxGetPr(prhs[6]);
	int n = (int)*mxGetPr(prhs[7]);

	const mwSize *dims;
	const mxArray *cell;
	const mxArray *cellElement;
	mwIndex jcell;
	cell = prhs[8];
	dims = mxGetDimensions(prhs[8]);

	std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > dF_dq_sp;
	std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > E;

	for (jcell = 0; jcell < dims[0]; jcell++) {
		cellElement = mxGetCell(cell, jcell);
		Eigen::MatrixXd element;
		igl::matlab::parse_rhs_double(&cellElement, element);
		dF_dq_sp.emplace_back(element);
	}

	cell = prhs[9];
	dims = mxGetDimensions(prhs[9]);
	for (jcell = 0; jcell < dims[0]; jcell++) {
		cellElement = mxGetCell(cell, jcell);
		Eigen::VectorXi element;
		igl::matlab::parse_rhs_index(&cellElement, element);
		E.emplace_back(element);
	}
	sim::vem3dmesh_neohookean_dq2(g, A, dF_dq, dM_dq, w, volumes, params, k, n, dF_dq_sp, E);

	//Eigen::SparseMatrixd g_sp(3 * n, 3 * n);
	//sim::vem3dmesh_neohookean_dq2_v2(g_sp, A, dF_dq, dM_dq, w, volumes, params, k, n, dF_dq_sp, E);


    igl::matlab::prepare_lhs_double(g, plhs);
}