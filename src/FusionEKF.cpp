#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
* Constructor.
*/
FusionEKF::FusionEKF() {
	is_initialized_ = false;

	previous_timestamp_ = 0;

	// initializing matrices
	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	Hj_ = MatrixXd(3, 4);

	//measurement covariance matrix - laser constant
	R_laser_ << 0.0225, 0,
		        0, 0.0225;

	//measurement covariance matrix - radar constant
	R_radar_ << 0.09, 0, 0,
				0, 0.0009, 0,
				0, 0, 0.09;

	//measurement matrix constant
	H_laser_ << 1, 0, 0, 0,
				0, 1, 0, 0;
	/**
	TODO:
	* Finish initializing the FusionEKF.
	* Set the process and measurement noises
	*/
	// Initializing P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1000, 0,
				0, 0, 0, 1000;

	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
				0, 1, 0, 1,
				0, 0, 1, 0,
				0, 0, 0, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

	/*****************************************************************************
	*  Initialization stage
	****************************************************************************/
	if (!is_initialized_) {
		/**
		TODO:
		* Initialize the state ekf_.x_ with the first measurement.
		* Create the covariance matrix.
		* Remember: you'll need to convert radar from polar to cartesian coordinates.
		*/
		// first measurement with four 1 values.
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << 1, 1, 1, 1;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			cout << "EKF : First measurement RADAR" << endl;
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			double rho = measurement_pack.raw_measurements_[0]; // range
			double phi = measurement_pack.raw_measurements_[1]; // bearing
			double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho. Coordinates convertion from polar to cartesian
			double x;
			double y;
			double vx;
			double vy;

			x = rho * cos(phi);
			if (x < 0.0001) {
				x = 0.0001;
			}
			y = rho * sin(phi);
			if (y < 0.0001) {
				y = 0.0001;
			}
			vx = rho_dot * cos(phi);
			vy = rho_dot * sin(phi);
			ekf_.x_ << x, y, vx, vy;
			cout << "EKF Radar initialization finished" << endl;
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			// Laser no need to change measurement matrix into cartesian
			cout << "EKF : First measurement LASER" << endl;
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
		}

		// Saving first timestamp in seconds for calculating interval
		previous_timestamp_ = measurement_pack.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		cout << "EKF Laser initialization finished" << endl;
		return;
	}

	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	/**
	TODO:
	* Update the state transition matrix F according to the new elapsed time.
	- Time is measured in seconds.
	* Update the process noise covariance matrix.
	* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	*/
	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = measurement_pack.timestamp_;

	// State transition matrix update
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	// Noise covariance matrix computation
	// Noise values from the task

	double noise_ax = 9.0;
	double noise_ay = 9.0;

	cout << "calculate dt4,dt3,dt2 " << endl;
	double dt4 = dt * dt* dt* dt;
	double dt3 = dt * dt* dt;
	double dt2 = dt * dt;
	cout << "calculate ekf_.Q_ " << endl;
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << (dt4 * noise_ax) / 4, 0, (dt3 * noise_ax) / 2, 0,
		0, (dt4 * noise_ay) / 4, 0, (dt3 * noise_ay) / 2,
		(dt3 * noise_ax) / 2, 0, dt2 * noise_ax, 0,
		0, (dt3 * noise_ay) / 2, 0, dt2 * noise_ay;
	cout << "start predict step " << endl;
	ekf_.Predict();

	/*****************************************************************************
	*  Update
	****************************************************************************/

	/**
	TODO:
	* Use the sensor type to perform the update step.
	* Update the state and covariance matrices.
	*/

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.R_ = R_radar_;
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	}
	else {
		// Laser updates
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}