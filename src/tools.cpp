#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// This code is more or less the same as explained on the conferences.
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,const vector<VectorXd> &ground_truth) {

	//TODO:
	/* Calculate the RMSE here.
	*/
	//from udacity 
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	if (estimations.size() != ground_truth.size() && estimations.size() == 0) {
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		// ... your code here
		VectorXd r;
		r = estimations[i] - ground_truth[i];
		r = r.array()*r.array();
		rmse = rmse + r;
	}

	//calculate the mean
	// ... your code here
	rmse = rmse / estimations.size();

	//calculate the squared root
	// ... your code here
	rmse = rmse.array().sqrt();
	//return the result
	return rmse;
}

// This code is more or less the same as explained on the conferences.
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3, 4);

	if (x_state.size() != 4) {
		cout << "ERROR - CalculateJacobian () - The state vector must have size 4." << endl;
		return Hj;
	}
	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	double c1 = px * px + py * py;
	double c2 = sqrt(c1);
	double c3 = (c1*c2);

	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "ERROR - CalculateJacobian () - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix for Radar measurement update. To Calculate the P and X
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;

	return Hj;
}