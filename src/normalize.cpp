/*
C++ codes for data normalization adapted from Seurat:
https://github.com/satijalab/seurat/blob/master/src/data_manipulation.cpp 
*/

#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> CosineNormSparse(Eigen::SparseMatrix<double> data, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  Eigen::SparseMatrix<double> data2 = data.cwiseProduct(data);
  Eigen::ArrayXd col2Sums = data2.transpose() * Eigen::VectorXd::Ones(data2.cols());
  col2Sums = col2Sums.sqrt();
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      data.coeffRef(it.row(), it.col()) = double(it.value()) / col2Sums[k];
    }
  }
  return data;
}

// [[Rcpp::export]]
Eigen::MatrixXd CosineNorm(Eigen::MatrixXd data, bool display_progress = true){
  Progress p(data.cols(), display_progress);
  Eigen::VectorXd colSums = data.colwise().sum();
  Eigen::MatrixXd data2 = data.cwiseProduct(data);
  Eigen::ArrayXd col2Sums = data2.transpose() * Eigen::VectorXd::Ones(data2.cols());
  col2Sums = col2Sums.sqrt();
  for (int k=0; k < data.cols(); ++k){
    p.increment();
    Eigen::ArrayXd c = data.col(k).array();
    data.col(k) = c / col2Sums[k];
  }
  return data;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> LogNormSparse(Eigen::SparseMatrix<double> data, int scale_factor, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      it.valueRef() = log1p(double(it.value()) / colSums[k] * scale_factor);
    }
  }
  return data;
}

// [[Rcpp::export]]
Eigen::MatrixXd LogNorm(Eigen::MatrixXd data, int scale_factor, bool display_progress = true){
  Progress p(data.cols(), display_progress);
  Eigen::VectorXd colSums = data.colwise().sum();
  for(int k=0; k < data.cols(); ++k){
    p.increment();
    Eigen::ArrayXd c = data.col(k).array();
    data.col(k) = log1p(c / colSums[k] * scale_factor);
  }
  return data;
}

// [[Rcpp::export]]
Eigen::MatrixXd FastRowScale(Eigen::MatrixXd mat, bool scale = true, bool center = true, double scale_max = 10, bool display_progress = true){
  Progress p(mat.rows(), display_progress);
  Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.rows(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.row(i).array();
    double rowMean = r.mean();
    double rowSdev = 1;
    if(scale == true){
      if(center == true){
        rowSdev = sqrt((r - rowMean).square().sum() / (mat.cols() - 1));
      }
      else{
        rowSdev = sqrt(r.square().sum() / (mat.cols() - 1));
      }
    }
    if(center == false){
      rowMean = 0;
    }
    scaled_mat.row(i) = (r - rowMean) / rowSdev;
    for(int s=0; s<scaled_mat.row(i).size(); ++s){
      if(scaled_mat(i, s) > scale_max){
        scaled_mat(i, s) = scale_max;
      }
    }
  }
  return scaled_mat;
}

// [[Rcpp::export]]
Eigen::MatrixXd FastSparseRowScale(Eigen::SparseMatrix<double> mat, bool scale = true, bool center = true, double scale_max = 10, bool display_progress = true){
  mat = mat.transpose();
  Progress p(mat.outerSize(), display_progress);
  Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colMean = 0;
    double colSdev = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      colMean += it.value();
    }
    colMean = colMean / mat.rows();
    if (scale == true){
      int nnZero = 0;
      if(center == true){
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          nnZero += 1;
          colSdev += pow((it.value() - colMean), 2);
        }
        colSdev += pow(colMean, 2) * (mat.rows() - nnZero);
      }
      else{
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          colSdev += pow(it.value(), 2);
        }
      }
      colSdev = sqrt(colSdev / (mat.rows() - 1));
    }
    else{
      colSdev = 1;
    }
    if(center == false){
      colMean = 0;
    }
    Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
    scaled_mat.col(k) = (col.array() - colMean) / colSdev;
    for(int s=0; s<scaled_mat.col(k).size(); ++s){
      if(scaled_mat(s,k) > scale_max){
        scaled_mat(s,k) = scale_max;
      }
    }
  }
  return scaled_mat.transpose();
}

// [[Rcpp::export]]
NumericVector vector_clip(NumericVector x, double min, double max){
  return clamp(min, x, max) ;
}
