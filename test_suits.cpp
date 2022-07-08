#include <iostream>
#include <vector>
#include "mmaOpt.h"
#include "Eigen/Eigen"
#include "matlab_utils.h"


/*This script is the "beam problem" from the MMA paper of Krister Svanberg.
**
**minimize 0.0624*(x(1) + x(2) + x(3) + x(4) + x(5))
**subject to 61 / (x(1) ^ 3) + 37 / (x(2) ^ 3) + 19 / (x(3) ^ 3) + 7 / (x(4) ^ 3) + 1 / (x(5) ^ 3) = < 1,
**	1 = < x(j) = < 10, for j = 1, .., 5.
*/
std::tuple<double, std::vector<double>, std::vector<double>, std::vector<double>> beam2(const std::vector<double>& xval) {
	double f0 = 0;
	std::vector<double> dfdx;
	for (int i = 0; i < xval.size(); i++) {
		f0 += 0.0624 * xval[i];
		dfdx.push_back(0.0624);
	}

	std::vector<double> gval;
	double g0 = 0;
	g0 += 61 / (pow(xval[0], 3)) +
		37 / (pow(xval[1], 3)) +
		19 / (pow(xval[2], 3)) +
		7 / (pow(xval[3], 3)) +
		1 / (pow(xval[4], 3)) - 1;
	gval.emplace_back(g0);
	std::vector<double> dgdx;
	dgdx.push_back(-3 * 61 / pow(xval[0], 4));
	dgdx.push_back(-3 * 37 / pow(xval[1], 4));
	dgdx.push_back(-3 * 19 / pow(xval[2], 4));
	dgdx.push_back(-3 * 7  / pow(xval[3], 4));
	dgdx.push_back(-3 * 1  / pow(xval[4], 4));

	return std::make_tuple(f0, dfdx, gval, dgdx);
}

std::tuple<double, std::vector<double>, std::vector<double>, std::vector<double>> Problem1(const std::vector<double>& xval) {
	std::vector<double> w({ 0.6020,0.2630,0.6541,0.6892,0.7482,0.4505,0.0838,0.2290,0.9133,0.1524 });
	double f0val = 0;
	for (int i = 0; i < w.size(); i++) f0val += w[i] * xval[i];
	std::vector<double> df0dx(xval.size(), 0);
	//df0dx[1] = 1;
	df0dx = w;
	std::vector<double> gval(2, 0);

	for (int i = 0; i < xval.size(); i++) { gval[0] += pow(xval[i], 2); }
	gval[0] -= pow(1.5, 2);

	for (int i = 0; i < xval.size(); i++) {
		if (i == 0) {
			gval[1] += pow(xval[0] - 2, 2);
		} else {
			gval[1] += pow(xval[i], 2);
		}
	}
	gval[1] -= pow(1.5, 2);

	std::vector<double> dgdx(xval.size() * 2, 0);
	for (int i = 0; i < xval.size(); i++) {
		if (i == 0) {
			dgdx[i] = 2 * xval[0];
			dgdx[i + xval.size()] = 2 * (xval[0] - 2);
		}
		else {
			dgdx[i] = 2 * xval[i];
			dgdx[i + xval.size()] = 2 * xval[i];
		}
	}
	return std::make_tuple(f0val, df0dx, gval, dgdx);
}

std::tuple<double, std::vector<double>, std::vector<double>, std::vector<double>> Problem2(
	const std::vector<double>& xval,
	const Eigen::Matrix<double, -1, -1>& weight,
	const Eigen::Matrix<double, -1, -1>& centers
) {
	int m = centers.cols();
	int n = weight.rows();
	Eigen::Matrix<double, -1, 1> w;
	w = weight.col(0);
	Eigen::Matrix<double, -1, 1> x(n, 1);
	std::copy(xval.begin(), xval.end(), x.data());
	double f0val = w.dot(x);
	std::vector<double> df0(w.begin(), w.end());
	
	std::vector<Eigen::Matrix<double, -1, 1>> sample_centers(m);
	for (int i = 0; i < sample_centers.size(); i++) {
		sample_centers[i] = centers.col(i);
	}

	std::vector<double> gval(m);
	for (int i = 0; i < m; i++) {
		gval[i] = (x - sample_centers[i]).squaredNorm() - pow(35, 2);
	}
	std::vector<double> dgdx(xval.size() * gval.size());
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			dgdx[i * n + j] = 2 * (x[j] - sample_centers[i][j]);
		}
	}

	return std::make_tuple(f0val, df0, gval, dgdx);
}


void test_mma(void) {
	// beam2
	if (0)
	{
		std::vector<double> xval(5, 5);
		std::vector<double> xmin(5, 1);
		std::vector<double> xmax(5, 10);
		std::vector<double> low(5), upp(5);
		std::vector<double> xold1(xval), xold2(xval);
		double a0 = 1;
		std::vector<double> a(1, 0);
		std::vector<double> c(1, 1000);
		std::vector<double> d(1, 1);
		double move = 1;
		double xch = 1;
		int itn = 0;
		while (itn++ <100 && xch >1e-5) {
			auto[f0, dfdx, gval, dgdx] = beam2(xval);
			mmasub(1, 5, itn,
				xval.data(), xmin.data(), xmax.data(), xold1.data(), xold2.data(),
				f0, dfdx.data(), gval.data(), dgdx.data(),
				low.data(), upp.data(), a0, a.data(), c.data(), d.data(), move);

			std::cout << "\033[32mxval = ";
			for (int i = 0; i < xval.size(); i++) { std::cout << xval[i] << " "; }
			xch = (Eigen::VectorXd::Map(xval.data(), 5) - Eigen::VectorXd::Map(xold1.data(), 5)).norm();
			std::cout << "\033[34mxch = " << xch;
			std::cout << "\033[0m" << std::endl;
		}
	} 
	else if (0) {
		std::vector<double> xval(10, 0);
		xval[0] = 1.15;
		std::vector<double> xmin(10, -4);
		std::vector<double> xmax(10, 8);
		std::vector<double> low(10), upp(10);
		std::vector<double> xold1(xval), xold2(xval);
		double a0 = 1;
		std::vector<double> a(2, 0);
		std::vector<double> c(2, 1000);
		std::vector<double> d(2, 1);
		double move = 1;
		double xch = 1;
		int itn = 0;
		while (itn++<100 && xch >1e-5) {
			std::cout << "\n\033[33m" << "* * * * * * * * * * * * * iter " << itn << " * * * * * * * * * * * * * \033[0m" << std::endl;
			auto[f0, dfdx, gval, dgdx] = Problem1(xval);
			mmasub(2, 10, itn,
				xval.data(), xmin.data(), xmax.data(), xold1.data(), xold2.data(),
				f0, dfdx.data(), gval.data(), dgdx.data(),
				low.data(), upp.data(), a0, a.data(), c.data(), d.data(), move);
			std::cout << "\033[32mxval = ";
			for (int i = 0; i < xval.size(); i++) { printf("%6.4lf ", xval[i]); }
			xch = (Eigen::VectorXd::Map(xval.data(), 10) - Eigen::VectorXd::Map(xold1.data(), 10)).norm();
			std::cout << "\033[34mxch = " << xch;
			std::cout << "\033[0m" << std::endl;
		}
	}
	else if (1) {
		Eigen::Matrix<double, -1, -1> weightvector = loadMatrixFromFile(
			"weightvector",
			"../data/problem2.mat"
		);
		Eigen::Matrix<double, -1, -1> centers = loadMatrixFromFile(
			"centers",
			"../data/problem2.mat"
		);
		//Eigen::Matrix<double, -1, -1> nstep = loadMatrixFromFile(
		//	"nstep",
		//	"D:/sources/GCMMA-MMA-code-1.5/nstep.mat"
		//);
		//nlineSearchStep.resize(nstep.cols());
		//for (int i = 0; i < nstep.cols(); i++) { nlineSearchStep[i] = nstep(0, i) + 0.5; }
		eigen2ConnectedMatlab("weightvector", weightvector);
		eigen2ConnectedMatlab("centers", centers);
		int m = centers.cols();
		int n = centers.rows();
		std::vector<double> xval(n, 0);
		std::vector<double> xmin(n, -4);
		std::vector<double> xmax(n, 8);
		std::vector<double> low(n), upp(n);
		std::vector<double> xold1(xval), xold2(xval);
		double a0 = 1;
		std::vector<double> a(m, 0);
		std::vector<double> c(m, 1000);
		std::vector<double> d(m, 1);
		double move = 1;
		double xch = 1;
		int itn = 0;
		while (itn++<100 && xch >1e-5) {
			std::cout << "\n\033[33m" << "* * * * * * * * * * * * * iter " << itn << " * * * * * * * * * * * * * \033[0m" << std::endl;
			auto[f0, dfdx, gval, dgdx] = Problem2(xval, weightvector, centers);
			mmasub(m, n, itn,
				xval.data(), xmin.data(), xmax.data(), xold1.data(), xold2.data(),
				f0, dfdx.data(), gval.data(), dgdx.data(),
				low.data(), upp.data(), a0, a.data(), c.data(), d.data(), move);
			std::cout << "\033[32mxval = ";
			for (int i = 0; i < 10; i++) { printf("%6.4lf ", xval[i]); }
			xch = (Eigen::VectorXd::Map(xval.data(), n) - Eigen::VectorXd::Map(xold1.data(), n)).norm();
			std::cout << "\033[34mxch = " << xch;
			std::cout << "\033[0m" << std::endl;
		}
	}
}

int main(int argc, char** argv) {
	//test_gVector();
	test_mma();
	return 0;
}
