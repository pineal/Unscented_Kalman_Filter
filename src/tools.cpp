#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	if (estimations.size() != ground_truth.size() || estimations.empty()) {
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	for (size_t i = 0; i < estimations.size(); ++i) {
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}
	
	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();
	cout << "RMSE: [ "<< rmse[0] << ", "<< rmse[1] << ", "<< rmse[2] << ", "<< rmse[3]<< " ]"<< endl;
	return rmse;
}

void Tools::Normalize(double & v) {
	v -= static_cast<int>(v / (2 * M_PI)) / (2 * M_PI);
};