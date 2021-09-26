#ifndef DLUT_SAE_PERIDYNAMIC_FEM_SHAPE_FUNCTIONS_20210528
#define DLUT_SAE_PERIDYNAMIC_FEM_SHAPE_FUNCTIONS_20210528

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			//	Beam单元形函数，其中Chi范围为[0,1]
			Eigen::Matrix<double, 6, 12>		N_SF_BEAM(double L, double Chi)
			{
				Eigen::Matrix<double, 6, 12> Res;
				Res.setZero();
				Res(0, 0) = 1 - Chi;
				Res(0, 6) = Chi;

				Res(1, 1) = 1 - 3 * Chi * Chi + 2 * Chi * Chi * Chi;
				Res(1, 5) = (Chi - 2 * Chi * Chi + Chi * Chi * Chi) * L;
				Res(1, 7) = 3 * Chi * Chi - 2 * Chi * Chi * Chi;
				Res(1, 11) = (Chi * Chi * Chi - Chi * Chi) * L;

				Res(2, 2) = 1 - 3 * Chi * Chi + 2 * Chi * Chi * Chi;
				Res(2, 4) = -(Chi - 2 * Chi * Chi + Chi * Chi * Chi) * L;
				Res(2, 8) = 3 * Chi * Chi - 2 * Chi * Chi * Chi;
				Res(2, 10) = -(Chi * Chi * Chi - Chi * Chi) * L;

				Res(3, 3) = 1 - Chi;
				Res(3, 9) = Chi;

				Res(4, 2) = -(-6 * Chi + 6 * Chi * Chi) / L;
				Res(4, 4) = (1 - 4 * Chi + 3 * Chi * Chi);
				Res(4, 8) = -(6 * Chi - 6 * Chi * Chi) / L;
				Res(4, 10) = (3 * Chi * Chi - 2 * Chi);

				Res(5, 1) = (-6 * Chi + 6 * Chi * Chi) / L;
				Res(5, 5) = (1 - 4 * Chi + 3 * Chi * Chi);
				Res(5, 7) = (6 * Chi - 6 * Chi * Chi) / L;
				Res(5, 11) = (3 * Chi * Chi - 2 * Chi);

				return Res;
			}

			//	Rod单元形函数，其中Chi范围为[0,1]
			Eigen::Matrix<double, 3, 6>			N_SF_ROD(double Chi)
			{
				Eigen::Matrix<double, 3, 6> Res;
				Res.setZero();
				Res(0, 0) = 1 - Chi;
				Res(0, 3) = Chi;
				Res(1, 1) = 1 - Chi;
				Res(1, 4) = Chi;
				Res(2, 2) = 1 - Chi;
				Res(2, 5) = Chi;

				return Res;
			}

			//	矩形单元Plane形函数，其中s和t范围均为[-1,1]
			Eigen::Matrix<double, 3, 12>		M_SF_RECTANGLE_PLANE(double s, double t)
			{
				Eigen::Matrix<double, 3, 12> Res;
				Res.setZero();
				double N1 = (1 - s) * (1 - t) / 4.0;
				double N2 = (1 + s) * (1 - t) / 4.0;
				double N3 = (1 + s) * (1 + t) / 4.0;
				double N4 = (1 - s) * (1 + t) / 4.0;

				Res(0, 0) = Res(1, 1) = Res(2, 2) = N1;
				Res(0, 3) = Res(1, 4) = Res(2, 5) = N2;
				Res(0, 6) = Res(1, 7) = Res(2, 8) = N3;
				Res(0, 9) = Res(1, 10) = Res(2, 11) = N4;

				return Res;
			}

			//	矩形单元Shell形函数，其中s和t范围均为[-1,1]
			Eigen::Matrix<double, 6, 6>		N_Ipq(double p, double q, double pI, double qI, double a, double b)
			{
				Eigen::Matrix<double, 6, 6>	 Res;
				Res.setZero();

				Res(0, 0) = (1 + p * pI) * (1 + q * qI) / 4.0;
				Res(1, 1) = Res(0, 0);

				Res(2, 2) = -((1 + pI * p) * (1 + qI * q) * (p * p + q * q - pI * p - qI * q - 2)) / 8.0;
				Res(2, 3) = b * qI * ((1 + pI * p) * (1 + qI * q) * (1 + qI * q) * (qI * q - 1)) / 8.0;
				Res(2, 4) = -a * pI * ((1 + pI * p) * (1 + pI * p) * (1 + qI * q) * (pI * p - 1)) / 8.0;

				Res(3, 2) = (-qI * (pI * p + 1) * (3 * q * q + p * p - pI * p - 3)) / (b * 8.0);
				Res(3, 3) = (b * (3 * qI * q - 1) * (1 + pI * p) * (1 + qI * q)) / (b * 8.0);
				Res(3, 4) = (-a * pI * qI * (1 + pI * p) * (1 + pI * p) * (pI * p - 1)) / (b * 8.0);

				Res(4, 2) = (-pI * (qI * q + 1) * (3 * p * p + q * q - qI * q - 3)) / (-a * 8.0);
				Res(4, 3) = (b * pI * qI * (qI * q + 1) * (qI * q + 1) * (qI * q - 1)) / (-a * 8.0);
				Res(4, 4) = (-a * (3 * pI * p - 1) * (pI * p + 1) * (qI * q + 1)) / (-a * 8.0);

				Res(5, 5) = Res(0, 0);

				return Res;
			}	
			Eigen::Matrix<double, 6, 24>	M_SF_RECTANGLE_SHELL(double s, double t, double a, double b)
			{
				Eigen::Matrix<double, 6, 24> Res;
				Res.setZero();

				double pI[4] = { -1, 1, 1, -1 };
				double qI[4] = { -1, -1, 1, 1 };

				Res.block(0, 0, 6, 6) = N_Ipq(s, t, pI[0], qI[0], a, b);
				Res.block(0, 6, 6, 6) = N_Ipq(s, t, pI[1], qI[1], a, b);
				Res.block(0, 12, 6, 6) = N_Ipq(s, t, pI[2], qI[2], a, b);
				Res.block(0, 18, 6, 6) = N_Ipq(s, t, pI[3], qI[3], a, b);

				return Res;
			}
		}
	}
}

#endif