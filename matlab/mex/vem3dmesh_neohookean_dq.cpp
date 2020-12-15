#include <mex.h>

//matlab hooks from libigl (taken directly from gptoolbox)
#include <mex.h>
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
#include <vem3dmesh_neohookean_dq.h>

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    /* variable declarations here */
    Eigen::VectorXd g;

    Eigen::MatrixXd A, dF_dq, dM_dq;
    Eigen::VectorXd w;
    Eigen::MatrixXd params;
    Eigen::VectorXd volumes; 

    igl::matlab::parse_rhs_double(prhs+0, A);
    igl::matlab::parse_rhs_double(prhs+1, dF_dq);
    igl::matlab::parse_rhs_double(prhs+2, dM_dq);
    igl::matlab::parse_rhs_double(prhs+3, w);
    igl::matlab::parse_rhs_double(prhs+4, volumes);
    igl::matlab::parse_rhs_double(prhs+5, params);
	int k = (int)*mxGetPr(prhs[6]);
	int n = (int)*mxGetPr(prhs[7]);

	sim::vem3dmesh_neohookean_dq(g, A, dF_dq, dM_dq, w, volumes, params, k, n);

    igl::matlab::prepare_lhs_double(g, plhs);
}