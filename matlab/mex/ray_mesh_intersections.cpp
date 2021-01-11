#include <igl/embree/EmbreeIntersector.h>
#include <igl/Hit.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#include <igl/matlab/MexStream.h>
#include <igl/parallel_for.h>

#include <Eigen/Core>
#include <limits>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;
  using namespace igl::embree;
  // Embree uses singles
  MatrixXf V,source,dir;
  MatrixXi F;
  mexErrMsgTxt(nrhs==4,"nrhs should == 4");

  parse_rhs_double(prhs+0,source);
  parse_rhs_double(prhs+1,dir);
  parse_rhs_double(prhs+2,V);
  parse_rhs_index( prhs+3,F);
  mexErrMsgTxt(source.cols()==3,"source must be #source by 3");
  mexErrMsgTxt(dir.cols()==3,"dir must be #dir by 3");
  mexErrMsgTxt(source.rows() == dir.rows(), "#source must equal #dir");
  mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");

  const int n = source.rows();

  int k = 16;
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> t(n, k);
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> ids(n, k);
  Eigen::Matrix<int, Eigen::Dynamic, 1> nhits(n);

  EmbreeIntersector ei;
  ei.init(V,F,true);

  //igl::parallel_for(n,[&](const int si)
  //for (int si = 0; si < n; si++)
  igl::parallel_for(n, [&](const int si)
  {
    Eigen::Vector3f s = source.row(si);
    Eigen::Vector3f d = dir.row(si);
    const float tnear = 1e-4f;
	std::vector<igl::Hit> hit;
	int num_rays = 0;
    if(!ei.intersectRay(s,d,hit,num_rays,tnear))
    {
	  nhits(si) = hit.size();
	  //std::cout << "si: " << si << " nrays: " << num_rays << " hitsize: "<< hit.size() << std::endl;
	  for (int i = 0; i < hit.size(); ++i) {
	    ids(si, i) = hit[i].id;
		t(si, i) = hit[i].t;
	  }
    }else
    {
	  nhits(si) = 0;
    }
  }
  );

  prepare_lhs_double(t, plhs + 0);
  prepare_lhs_double(nhits, plhs + 1);

  if (nlhs > 2)
  {
    prepare_lhs_index(ids, plhs + 2);
  }
}
