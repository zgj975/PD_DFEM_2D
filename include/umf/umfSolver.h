#ifndef DLUT_SAE_PERIDYNAMIC_UMF_SOLVER_HXX_20190903
#define DLUT_SAE_PERIDYNAMIC_UMF_SOLVER_HXX_20190903

#include <iostream>
#include <vector>
#include <map>
#include "umfpack.h"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

#include "Eigen/Dense"
#include "Eigen/SparseLU"
#include "Eigen/Eigenvalues"
#include "pd_base_toolkit.hxx"

using namespace DLUT::SAE::PERIDYNAMIC;
using namespace Eigen;
using namespace std;

void umf_solver(const vector< map<int, double> >& A, vector<double>& x, const vector<double>& b)
{
	vector<int> Ap;
	vector<int> Ai;
	vector<double> Ax;

	int n = (int)(A.size());
	Ap.resize(n + 1);
	Ap[0] = 0;

	int num = 0;
	for (const map<int, double>& col_datas : A)
	{
		for (const pair<int, double>& j_v : col_datas)
		{
			if (j_v.second > ERR_VALUE || j_v.second < -ERR_VALUE)
			{
				num++;				
			}
		}
	}
	Ai.resize(num);
	Ax.resize(num);

	int k = 0;
	for (int i = 0; i < (int)(A.size()); ++i)
	{
		Ap[i + 1] = Ap[i];
		for (const pair<int, double>& j_v : A[i])
		{
			if (j_v.second > ERR_VALUE || j_v.second < -ERR_VALUE)
			{
				Ax[k] = j_v.second;

				Ai[k] = (int)(j_v.first);
				k++;
				Ap[i + 1]++;
			}
		}
	}

	double* null = (double*)NULL;
	void* Symbolic, * Numeric;

	umfpack_di_symbolic(n, n, &Ap[0], &Ai[0], &Ax[0], &Symbolic, null, null);
	umfpack_di_numeric(&Ap[0], &Ai[0], &Ax[0], Symbolic, &Numeric, null, null);
	umfpack_di_free_symbolic(&Symbolic);
	umfpack_di_solve(UMFPACK_A, &Ap[0], &Ai[0], &Ax[0], &x[0], &b[0], Numeric, null, null);
	umfpack_di_free_numeric(&Numeric);
}

#endif