#pragma once

#include <Eigen/Dense>

namespace polyfem::mesh
{
	void apply_slim(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::MatrixXd &V_new);
}