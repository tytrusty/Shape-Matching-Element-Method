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

#include <vem_ext_force.h>

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    /* variable declarations here */
    Eigen::VectorXd g;

    Eigen::VectorXd ext_force, mass;

	EigenMatrixList<double> Y;
	EigenSparseMatrixList<double> Y_S;

	igl::matlab::parse_rhs_double(prhs+0, ext_force);
	igl::matlab::parse_rhs_double(prhs+1, mass);
	sim::parse_cell(prhs, 2, Y);
	sim::parse_cell_sparse(prhs, 3, Y_S);
	
	//sim::vem3dmesh_neohookean_dq(g, x, c, volumes, params, dF_dc, dF_dc_S, ME, L, k, n, d, x0_coms_size, k_stability);
	sim::vem_ext_force(g, ext_force, mass, Y, Y_S);

	igl::matlab::prepare_lhs_double(g, plhs);
}