#pragma once

#include <memory>
#include <Eigen/Core>

namespace polyfem::solver
{
	/** This parameterize a function f : x -> y
	 * and provides the chain rule with respect to previous gradients
	 */
	class Parametrization
	{
	public:
		Parametrization() {}
		virtual ~Parametrization() {}

		virtual Eigen::VectorXd inverse_eval(const Eigen::VectorXd &y);

		virtual int size(const int x_size) const = 0; // just for verification
		virtual Eigen::VectorXd eval(const Eigen::VectorXd &x) const = 0;
		virtual Eigen::VectorXd apply_jacobian(const Eigen::VectorXd &grad_full, const Eigen::VectorXd &x) const = 0;
	};

	class IndexedParametrization : public Parametrization
	{
	public:
		IndexedParametrization() {}
		virtual ~IndexedParametrization() {}

		void set_output_indexing(const Eigen::VectorXi &output_indexing) { output_indexing_ = output_indexing; }
		Eigen::VectorXi get_output_indexing(const Eigen::VectorXd &x) const;

	protected:
		Eigen::VectorXi output_indexing_;
	};

	class CompositeParametrization : public IndexedParametrization
	{
	public:
		CompositeParametrization() {}
		CompositeParametrization(const std::vector<std::shared_ptr<Parametrization>> &parametrizations) : parametrizations_(parametrizations) {}
		virtual ~CompositeParametrization() {}

		Eigen::VectorXd inverse_eval(const Eigen::VectorXd &y) override;

		int size(const int x_size) const override;
		Eigen::VectorXd eval(const Eigen::VectorXd &x) const override;
		Eigen::VectorXd apply_jacobian(const Eigen::VectorXd &grad_full, const Eigen::VectorXd &x) const override;

	private:
		std::vector<std::shared_ptr<Parametrization>> parametrizations_;
	};
} // namespace polyfem::solver