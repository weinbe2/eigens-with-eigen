// Copyright (c) 2017 Evan S Weinberg
// A reference piece of code which computes matrix elements
// of a reference sparse gauged Laplace operator to fill a dense
// Eigen matrix, then computes the spectrum and prints
// the Eigenvalues.

// This code is for complex, Hermitian matrices.
// This code lives on github at github.com/weinbe2/eigens-with-eigen/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <complex>
#include <random>

// Borrow dense matrix eigenvalue routines.
#include <Eigen/Dense>

using namespace std; 
using namespace Eigen;


// This is just for convenience. By default, all matrices
// in Eigen are column major. You can replace "Dynamic"
// with a specific number to template just one size.
// You can also ust use "MatrixXcd".
typedef Matrix<complex<double>, Dynamic, Dynamic, ColMajor> cMatrix;

// Reference 1-D Laplace function.
void laplace_1d_gauge(complex<double>* out, complex<double>* in, const int L, const double m2, complex<double>* field);

int main(int argc, char** argv)
{  
  complex<double> *in_cplx;
  complex<double> *out_cplx;

  // Set output precision to be long.
  cout << setprecision(10);

  // Basic information about the lattice.
  const int length = 4;
  const double m_sq = 0.001;

  // Print the basic info.
  std::cout << "1D Laplace operator, length " << length << ", mass squared " << m_sq << ", periodic boundary conditions.\n";
  std::cout << "Change the length and mass by modifying the source.\n";

  // Allocate a random phase field.
  complex<double>* field = new complex<double>[length];
  std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
  std::uniform_real_distribution<> dist(-3.1415926535, 3.1415926535);
  for (int i = 0; i < length; i++)
  {
    // As a sanity check, you can set these all equal to 1.0.
    field[i] = exp(complex<double>(0.0, dist(generator)));
  }
  
  // Allocate.
  in_cplx = new complex<double>[length];
  out_cplx = new complex<double>[length];

  // Zero out.
  for (int i = 0; i < length; i++)
  {
    in_cplx[i] = out_cplx[i] = 0.0;
  }

  /////////////////////////////
  // COMPLEX, HERMITIAN CASE //
  /////////////////////////////

  std::cout << "Complex, Hermitian case.\n\n";

  // Allocate a sufficiently gigantic matrix.
  cMatrix mat_cplx = cMatrix::Zero(length, length);

  // Form matrix elements. This is where it's important that
  // dMatrix is column major.
  for (int i = 0; i < length; i++)
  {
    // Set a point on the rhs for a matrix element.
    // If appropriate, zero out the previous point.
    if (i > 0)
    {
      in_cplx[i-1] = 0.0;
    }
    in_cplx[i] = 1.0;

    // Zero out the "out" vector. I defined "laplace_1d_gauge" to
    // not require this, but I put this here for generality.
    for (int j = 0; j < length; j++)
      out_cplx[j] = 0.0;

    // PUT YOUR MAT-VEC HERE.
    laplace_1d_gauge(out_cplx, in_cplx, length, m_sq, field);

    // Copy your output into the right memory location.
    // If your data layout supports it, you can also pass
    // "mptr" directly as your "output vector" when you call
    // your mat-vec.
    complex<double>* mptr = &(mat_cplx(i*length));
    
    for (int j = 0; j < length; j++)
    {
      mptr[j] = out_cplx[j];
    }
  }

  // We've now formed the dense matrix. We print it here
  // as a sanity check if it's small enough.
  if (length <= 16)
  {
    std::cout << mat_cplx << "\n";
  }

  // Get the eigenvalues and eigenvectors.
  SelfAdjointEigenSolver< cMatrix > eigsolve_cplx(length);
  eigsolve_cplx.compute(mat_cplx);

  // Remark: if you only want the eigenvalues, you can call
  // eigsolve_cplx.compute(mat_cplx, EigenvaluesOnly);

  // Print the eigenvalues.
  cMatrix evals = eigsolve_cplx.eigenvalues();
  std::cout << "The eigenvalues are:\n" << evals << "\n\n";
  // You can also index individual eigenvalues as "evals(i)", where
  // "i" zero-indexes the eigenvalues.

  // Print the eigenvectors if the matrix is small enough.
  // As a remark, this also shows you how to access the eigenvectors. 
  if (length <= 16)
  {
    for (int i = 0; i < length; i++)
    {
      // You can use "VectorXd" as the type instead.
      cMatrix evec = eigsolve_cplx.eigenvectors().col(i);
      
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
  delete[] field;
  delete[] in_cplx;
  delete[] out_cplx;

  return 0;
}



// Reference 1-D Laplace function.
void laplace_1d_gauge(complex<double>* out, complex<double>* in, const int L, const double m2, complex<double>* field)
{
  // Loop over all sites, applying fields.
  for (int i = 0; i < L; i++)
  {
    out[i] = (2+m2)*in[i] - conj(field[(i-1+L)%L])*in[(i-1+L)%L] - field[i]*in[(i+1)%L];
  }

  // Done.
  return;
}

