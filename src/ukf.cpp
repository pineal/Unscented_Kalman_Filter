#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	this->use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	this->use_radar_ = true;

	// initial state vector
	this->x_ = VectorXd(5);

	// initial covariance matrix
	this->P_ = MatrixXd(5, 5);
	P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;


	// Process noise standard deviation longitudinal acceleration in m/s^2
	this->std_a_ = 10;

	// Process noise standard deviation yaw acceleration in rad/s^2
	this->std_yawdd_ = 30;

	//DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
	// Laser measurement noise standard deviation position1 in m
	this->std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	this->std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	this->std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	this->std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	this->std_radrd_ = 0.3;

	//DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
	this->is_initialized_ = false;
	this->time_us_ = 0.0;
	//set state dimension
	this->n_x_ = 5;

	// set augmented dimension
	this->n_aug_ = n_x_ + 2;

	this->Xsig_pred_ = MatrixXd(n_x_, n_aug_ * 2 + 1);

	// set sigma point spreading parameter
	this->lambda_ = 3 - n_x_;

	//set weights
	this->weights_ = VectorXd(n_aug_ * 2 + 1);
	weights_.fill(0.5 / (n_aug_ + lambda_));
	weights_(0) = lambda_ / (n_aug_ + lambda_);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	// update time stamp
	bool isRADAR = meas_package.sensor_type_ == MeasurementPackage::RADAR? true : false;

	// initialize in the first time
	if (!is_initialized_) {
		if (isRADAR) {
			double rho = meas_package.raw_measurements_[0];
			double phi = meas_package.raw_measurements_[1];
			double rho_dot = meas_package.raw_measurements_[2];
			double x = rho * cos(phi);
			double y = rho * sin(phi);
			double vx = rho_dot * cos(phi);
			double vy = rho_dot * sin(phi);
			double v = sqrt(vx * vx + vy * vy);
			this->x_ << x, y, v, 0, 0;
		}
		else {
			x_ << 	meas_package.raw_measurements_[0], 
					meas_package.raw_measurements_[1],
					0,
					0,
					0;
		}
		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		return;
	}

	// prediction UKF
	double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;
	Prediction(dt);

	// update UKF 
	if (isRADAR) {
		UpdateRadar(meas_package);
	} 
	else {
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	int n_sig = 2 * n_aug_ + 1;
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_,n_x_) = P_;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

	// generate sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig);
	Xsig_aug.col(0) = x_aug;

	MatrixXd A = P_aug.llt().matrixL();
	double gamma = sqrt(lambda_ + n_aug_);
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + gamma * A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - gamma * A.col(i);
	}
	// predict sigma points
	for (int i = 0; i< n_sig; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0,i);
		double p_y = Xsig_aug(1,i);
		double v = Xsig_aug(2,i);
		double yaw = Xsig_aug(3,i);
		double yawd = Xsig_aug(4,i);
		double nu_a = Xsig_aug(5,i);
		double nu_yawdd = Xsig_aug(6,i);

		//avoid division by zero
    if (fabs(yawd) > COV_THREAD) {
        Xsig_pred_(0,i) = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw)) +
                          (0.5 * delta_t * delta_t * cos(yawd) * nu_a);
        Xsig_pred_(1,i) = p_y + v/yawd * (-cos(yaw + yawd * delta_t) + cos(yaw)) +
                          (0.5 * delta_t * delta_t * sin(yawd) * nu_a);
    }
    else {
        Xsig_pred_(0,i) = p_x + v*delta_t*cos(yaw) +
                          (0.5 * delta_t * delta_t * cos(yawd) * nu_a);
        Xsig_pred_(1,i) = p_y + v*delta_t*sin(yaw) +
                          (0.5 * delta_t * delta_t * sin(yawd) * nu_a);
    }

    Xsig_pred_(2,i) = v + nu_a * delta_t;
    Xsig_pred_(3,i) = yaw + yawd * delta_t + 0.5 * nu_yawdd * delta_t * delta_t;
    Xsig_pred_(4,i) = yawd + nu_yawdd * delta_t;

	}

	//predicted state mean
	x_ = Xsig_pred_ * weights_;

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < n_sig; i++) {

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//angle normalization
		Tools::Normalize(x_diff, 3);
		P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

	VectorXd z = meas_package.raw_measurements_;

	int n_z = 2;
	
	MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_aug_ * 2 + 1);

	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_ * std_laspx_, 0,
	  0,std_laspy_ * std_laspy_; 

	S += R;

	// Update State
	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		Tools::Normalize(z_diff, 1);

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		Tools::Normalize(x_diff, 3);

		Tc = Tc + this->weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	Tools::Normalize(z_diff, 1);

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	VectorXd z = meas_package.raw_measurements_;

	int n_z = 3;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v   = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                         //r
		Zsig(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);     //r_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		Tools::Normalize(z_diff, 1);

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R <<  std_radr_*std_radr_, 0, 0,
	  0, std_radphi_*std_radphi_, 0,
	  0, 0, std_radrd_*std_radrd_;
	S = S + R;



	// build cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		Tools::Normalize(z_diff, 1);
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		Tools::Normalize(x_diff, 3);
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	Tools::Normalize(z_diff, 1);
	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();


	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
