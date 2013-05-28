#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "../Blitz.h"
#include <blitz/array.h>
#include "../MatrixOps/MatrixOps.h"

using namespace blitz;

void FitEllipse(Array<double,2> &A)
{
  const int M = A.rows();
  const int D = A.cols();

  Array<double,1> avg(D);
  avg = 0;
  for (int k = 0; k < M; k++)
    for (int d = 0; d < D; d++)
      avg(d) += A(k,d);
  avg = avg/M;

  Array<double,2> C(D,D); // Will be the covariance matrix.
  C = 0;                  // Eigenvectors will be major & minor axes. Eigenvalues will be lengths of axes, squared.
  for (int k = 0; k < M; k++)
    for (int d1 = 0; d1 < D; d1++)
      for (int d2 = 0; d2 < D; d2++)
        C(d1,d2) += (A(k,d1) - avg(d1)) * (A(k,d2) - avg(d2));
  C = C/(M-1);

  Array<double,1> Vals(D);
  Array<double,2> Vectors(D,D);
  SymmEigenPairs(C,D,Vals,Vectors);

  cout << avg << C << Vals << Vectors << endl;
}

int main(int argc, char* argv[])
{
  string CoordFile = argv[1];
  std::ifstream ifs(CoordFile.c_str());

  const int M = 512;
  const int D = 3;

  Array<TinyVector<double,D>,1> A(M);
  std::string line;
  int i = 0;
  while(std::getline(ifs,line)) {
    istringstream iss(line);
    TinyVector<double,D> a;
    int j = 0;
    do {
        string sub;
        iss >> sub;
        double x = atof(sub.c_str());
        if (x!=0)
          a(j) = x;
        j++;
    } while (iss);
    A(i) = a;
    i++;
  }


  //for (int k = 0; k < M; k++)
  //  cout << k << " " << A(k) << endl;

  Array<double,2> B(M,D);
  for (int k = 0; k < M; k++)
    for (int d = 0; d < D; d++)
      B(k,d) = A(k)(d);
  FitEllipse(B);

  return 0;
}
