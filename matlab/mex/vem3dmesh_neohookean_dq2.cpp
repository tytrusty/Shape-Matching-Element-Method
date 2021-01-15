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

template <typename T>
using EigenMatrixList = std::vector<Eigen::MatrixXx<T>, Eigen::aligned_allocator<Eigen::MatrixXx<T>>>;

template <typename T>
using EigenVectorList = std::vector<Eigen::VectorXx<T>, Eigen::aligned_allocator<Eigen::VectorXx<T>>>;


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

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    /* variable declarations here */
    Eigen::MatrixXd g;

    Eigen::MatrixXd params;
    Eigen::VectorXd c, volumes; 

	std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > dF_dc;
	std::vector<Eigen::VectorXi, Eigen::aligned_allocator<Eigen::VectorXi> > W_I;

    igl::matlab::parse_rhs_double(prhs+0, c);
    igl::matlab::parse_rhs_double(prhs+1, volumes);
    igl::matlab::parse_rhs_double(prhs+2, params);
	parse_cell(prhs, 3, dF_dc);
	parse_cell_index(prhs, 4, W_I);

	int k = (int)*mxGetPr(prhs[5]);
	int n = (int)*mxGetPr(prhs[6]);
	int m = (int)*mxGetPr(prhs[7]);

	sim::vem3dmesh_neohookean_dq2(g, c, volumes, params, dF_dc, W_I, k, n, m);

    igl::matlab::prepare_lhs_double(g, plhs);
}