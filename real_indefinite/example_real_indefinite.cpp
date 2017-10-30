// Copyright (c) 2017 Evan S Weinberg
// A reference piece of code which computes matrix elements
// of a reference sparse first derivative operator to fill a dense
// Eigen matrix, then computes the spectrum and prints
// the Eigenvalues.

// This code is for real, indefinite matrices.
// This code lives on github at github.com/weinbe2/eigens-with-eigen/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

// Borrow dense matrix eigenvalue routines.
#include <Eigen/Dense>

using namespace std; 
using namespace Eigen;


// This is just for convenience. By default, all matrices
// in Eigen are column major. You can replace "Dynamic"
// with a specific number to template just one size.
// You can also ust use "MatrixXd".
typedef Matrix<double, Dynamic, Dynamic, ColMajor> dMatrix;
// Likewise, we can use "MatrixXcd".
typedef Matrix<std::complex<double>, Dynamic, Dynamic, ColMajor> cMatrix;
// We need both real and complex because, while the matrix itself is real,
// the eigenvalues and eigenvectors can generically be complex.


// Reference 1-D central difference operator
void central_diff_1d(double* out, double* in, const int L, const double m2);

int main(int argc, char** argv)
{  
  double *in_real;
  double *out_real;
  complex<double> *out_cplx; // needed for eigenvectors

  // Set output precision to be long.
  cout << setprecision(10);

  // Basic information about the lattice.
  const int length = 8;
  const double mass = 0.001;

  // Print the basic info.
  std::cout << "1D central difference operator, length " << length << ", mass " << mass << ", periodic boundary conditions.\n";
  std::cout << "Change the length and mass by modifying the source.\n";
  
  // Allocate.
  in_real = new double[length];
  out_real = new double[length];
  out_cplx = new complex<double>[length];

  // Zero out.
  for (int i = 0; i < length; i++)
  {
    in_real[i] = out_real[i] = 0.0;
  }

  ///////////////////////////
  // REAL, INDEFINITE CASE //
  ///////////////////////////

  std::cout << "Real, Indefinite case.\n\n";

  // Allocate a sufficiently gigantic matrix.
  dMatrix mat_real = dMatrix::Zero(length, length);

  // Form matrix elements. This is where it's important that
  // dMatrix is column major.
  for (int i = 0; i < length; i++)
  {
    // Set a point on the rhs for a matrix element.
    // If appropriate, zero out the previous point.
    if (i > 0)
    {
      in_real[i-1] = 0.0;
    }
    in_real[i] = 1.0;

    // Zero out the "out" vector. I defined "laplace_1d" to
    // not require this, but I put this here for generality.
    for (int j = 0; j < length; j++)
      out_real[j] = 0.0;

    // PUT YOUR MAT-VEC HERE.
    central_diff_1d(out_real, in_real, length, mass);

    // Copy your output into the right memory location.
    // If your data layout supports it, you can also pass
    // "mptr" directly as your "output vector" when you call
    // your mat-vec.
    double* mptr = &(mat_real(i*length));
    
    for (int j = 0; j < length; j++)
    {
      mptr[j] = out_real[j];
    }
  }

  // We've now formed the dense matrix. We print it here
  // as a sanity check if it's small enough.
  if (length <= 16)
  {
    std::cout << mat_real << "\n";
  }

  // Get the eigenvalues and eigenvectors.
  EigenSolver< dMatrix > eigsolve_real_indef(length);
  eigsolve_real_indef.compute(mat_real);

  // Remark: if you only want the eigenvalues, you can call
  // eigsolve_real.compute(mat_real, EigenvaluesOnly);

  // Print the eigenvalues.
  // Note that the eigenvalues are of type "cMatrix" as typedef'd above or
  // MatrixXcd generically.
  cMatrix evals = eigsolve_real_indef.eigenvalues();
  std::cout << "The eigenvalues are:\n" << evals << "\n\n";
  // You can also index individual eigenvalues as "evals(i)", where
  // "i" zero-indexes the eigenvalues.

  // Print the eigenvectors if the matrix is small enough.
  // As a remark, this also shows you how to access the eigenvectors. 
  if (length <= 16)
  {
    for (int i = 0; i < length; i++)
    {
      // You can use "VectorXcd" as the type instead.
      cMatrix evec = eigsolve_real_indef.eigenvectors().col(i);
      
      // Print the eigenvector.
      std::cout << "Eigenvector " << i << " equals:\n" << evec << "\n\n";

      // You can also copy the eigenvector into another array as such:
      for (int j = 0; j < length; j++)
      {
        out_cplx[j] = evec(j);
      }
    }
  }

  // Clean up.
  delete[] in_real;
  delete[] out_real;
  delete[] out_cplx;

  return 0;
}



// Reference 1-D central difference function.
void central_diff_1d(double* out, double* in, const int L, const double mass)
{
  // Central difference operator with periodic boundary conditions.
  for (int i = 0; i < L; i++)
  {
    out[i] = mass*in[i] + in[(i+1)%L] - in[(i-1+L)%L];
  }

  // Done.
  return;
}

