#ifndef SIM_OPTIMIZATION_GRADIENT_BASED_H
#define SIM_OPTIMIZATION_GRADIENT_BASED_H

#include <BartelsTypes.h>
#include <EigenTypes.h>
#include <Eigen/SVD> 
#include <linesearch_backtracking_bisection.h>

//gradient-based optimization schemes
namespace sim {
    
    //1st-order
    //unconstrained gradient descent with simple bisection, backtracking
     
    //2nd-order
    //unconstrained newtons method with simple bisection, backtracking
    //x - initial guess for solver, contains final value on exit
    //f - f(x) returns the value for the energy at x
    //g - g(tmp_g, x)  returns the value of the gradient in tmp_g
    //H - H(tmp_H, x) returns the value of the hessian in tmp_H
    //Solver - linear solver to use for search direction calculation
    //tmp_g - work space for computing the gradient during newton's method
    //tmp_H - work space for computing the gradient during newton's method
    //tmp_d - work space for storing the search direction during newton's method
    //gradient_tol - consider converged if gradient is less than this tolerance
    //max_iterations - max number of newton iterations
    //max_iterations_ls - max number of line search iterations
    //alpha - initial percentage of step size for line search (default 1.0)
    //c - sufficient decrease toloerance (default 1e-8)
    //p -  
    //opt_call - opt_call(auto &x, auto &g, auto &H, SolverExitStatus &searchstatus) optimization callback function, default does nothing
    //ls_call -  ls_call(auto &x) line search callback function (default does nothing)
    template <class Energy, class Gradient, class Hessian,
                int GRowsAtCompileTime,
                int GColsAtCompileTime,
                int GOptions = 0,
                int GMaxRowsAtCompileTime=GRowsAtCompileTime,
                int GMaxColsAtCompileTime=GColsAtCompileTime,
                typename Scalar, 
                class OptimizationCallback = decltype(default_optimization_callback),
                class LineSearchCallback = decltype(default_linesearch_callback) >
    inline SolverExitStatus newtons_method_bisection(Eigen::VectorXd &x,
                                                 Energy &f, 
                                                 Gradient &g, 
                                                 Hessian &H, 
                                                 Eigen::Matrix<Scalar, GRowsAtCompileTime, GColsAtCompileTime, GOptions, GMaxRowsAtCompileTime, GMaxColsAtCompileTime> &tmp_g,
                                                 Eigen::Matrix<Scalar, GRowsAtCompileTime, GRowsAtCompileTime, GOptions, GMaxRowsAtCompileTime, GMaxRowsAtCompileTime> &tmp_H,
                                                 Eigen::Matrix<Scalar, GRowsAtCompileTime, GColsAtCompileTime, GOptions, GMaxRowsAtCompileTime, GMaxColsAtCompileTime> &tmp_d,
                                                 Scalar gradient_tol = 1e-3,
                                                 unsigned int max_iterations = 100,
                                                 unsigned int max_iterations_ls = 100,
                                                 Scalar alpha = 1.0,
                                                 Scalar c = 1e-8,
                                                 Scalar p = 0.5,
                                                 const OptimizationCallback opt_call = default_optimization_callback,
                                                 const LineSearchCallback ls_call = default_linesearch_callback) 
                                                 
    {
        unsigned int iteration_count = 0;

        if(x.rows() == 0) {
            std::cout<<"Empty point sent to newtons search, quitting.\n";
            exit(1);
        }
        
        if(tmp_g.rows() != x.rows()) {
            tmp_g.resize(x.rows(),1);
        }   

        
        if((tmp_H.rows() != x.rows()) || (tmp_H.cols() != x.rows()) ) {
            tmp_H.resize(x.rows(), x.rows()); 
        }

        do {

            //compute gradient
            g(tmp_g, x);

            //check for convergence and return if converged 
            if(tmp_g.transpose()*tmp_g < gradient_tol*gradient_tol) 
                return SolverExitStatus::CONVERGED;

            //compute hessian
            H(tmp_H, x);

            // std::cout << "before solve " << std::endl;

            // std::cout << "H = " << tmp_H.block(0, 0, 10, 10) << std::endl;
            // std::cout << "g = " << tmp_g << std::endl;
             
            //compute search direction
            tmp_d = tmp_H.colPivHouseholderQr().solve(-tmp_g);
            // std::cout << "tmp_d:  " << tmp_d.block(0, 0, 10, 1) << std::endl;
            // std::cout << "after solve " << std::endl;
            double relative_error = (tmp_H*tmp_d - (-tmp_g)).norm() / (-tmp_g).norm(); // norm() is L2 norm
            // std::cout << "relative_error = " << relative_error << std::endl;
            if (relative_error > 1.0) {
                std::cout<<"Computing search direction failed in newtons method \n";
                exit(1);
            }

            // if(!search_direction_newton_sparse(tmp_d, tmp_g, tmp_H, solver))
            // {
            //     std::cout<<"Computing search direction failed in newtons method \n";
            //     exit(1);
            // }

            //line search 
            SolverExitStatus ls_status = linesearch_backtracking_bisection(x, tmp_d, f, tmp_g, 
                                                                           max_iterations_ls, 
                                                                           alpha, c, p,
                                                                           ls_call);

            opt_call(x, tmp_g, tmp_H, ls_status);

            iteration_count++;

        } while(iteration_count <= max_iterations);

        return SolverExitStatus::MAX_ITERATIONS_REACHED;
        
    }

    
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/optimization_gradient_based.cpp>
#endif

#endif
