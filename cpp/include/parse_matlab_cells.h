#ifndef SIM_PARSE_MATLAB_CELLS
#define SIM_PARSE_MATLAB_CELLS

#include <Eigen/Dense>
#include <EigenTypes.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>

namespace sim {

	// Parse cells composed of dense matrices
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

	// Parse cells composed of dense vectors
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

	// Parse cells composed of sparse matrices
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

}

#endif // SIM_PARSE_MATLAB_CELLS