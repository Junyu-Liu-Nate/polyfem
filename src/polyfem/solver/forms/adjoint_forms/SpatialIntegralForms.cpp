#include "SpatialIntegralForms.hpp"
#include <polyfem/utils/CubicHermiteSplineParametrization.hpp>
#include <polyfem/io/Evaluator.hpp>

using namespace polyfem::utils;

namespace polyfem::solver
{
	double SpatialIntegralForm::value_unweighted(const Eigen::VectorXd &x) const
	{
		assert(time_step_ < state_.diff_cached.size());
		return AdjointTools::integrate_objective(state_, get_integral_functional(), state_.diff_cached[time_step_].u, ids_, spatial_integral_type_, time_step_);
	}

	void SpatialIntegralForm::compute_partial_gradient_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
	{
		assert(time_step_ < state_.diff_cached.size());
		gradv.setZero(x.size());
		for (const auto &param_map : variable_to_simulations_)
		{
			const auto &parametrization = param_map->get_parametrization();
			const auto &state = param_map->get_state();
			const auto &param_type = param_map->get_parameter_type();

			if (&state != &state_)
				continue;

			Eigen::VectorXd term;
			if (param_type == ParameterType::Shape)
				AdjointTools::compute_shape_derivative_functional_term(state_, state_.diff_cached[time_step_].u, get_integral_functional(), ids_, spatial_integral_type_, term, time_step_);
			else if (param_type == ParameterType::MacroStrain)
				AdjointTools::compute_macro_strain_derivative_functional_term(state_, state_.diff_cached[time_step_].u, get_integral_functional(), ids_, spatial_integral_type_, term, time_step_);
      
			if (term.size() > 0)
				gradv += parametrization.apply_jacobian(term, x);
		}
	}

	Eigen::VectorXd SpatialIntegralForm::compute_adjoint_rhs_unweighted_step(const Eigen::VectorXd &x, const State &state)
	{
		if (&state != &state_)
			return Eigen::VectorXd::Zero(state.ndof(), state.diff_cached.size());

		assert(time_step_ < state_.diff_cached.size());

		Eigen::VectorXd rhs;
		AdjointTools::dJ_du_step(state, get_integral_functional(), state.diff_cached[time_step_].u, ids_, spatial_integral_type_, time_step_, rhs);

		return rhs;
	}

	// TODO: call local assemblers instead
	IntegrableFunctional StressForm::get_integral_functional() const
	{
		IntegrableFunctional j;

		const std::string formulation = state_.formulation();
		const int power = in_power_;

		j.set_j([formulation, power](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
			val.setZero(grad_u.rows(), 1);
			Eigen::MatrixXd grad_u_q, stress;
			for (int q = 0; q < grad_u.rows(); q++)
			{
				if (formulation == "Laplacian")
				{
					stress = grad_u.row(q);
				}
				else if (formulation == "LinearElasticity")
				{
					vector2matrix(grad_u.row(q), grad_u_q);
					stress = mu(q) * (grad_u_q + grad_u_q.transpose()) + lambda(q) * grad_u_q.trace() * Eigen::MatrixXd::Identity(grad_u_q.rows(), grad_u_q.cols());
				}
				else if (formulation == "NeoHookean")
				{
					vector2matrix(grad_u.row(q), grad_u_q);
					Eigen::MatrixXd def_grad = Eigen::MatrixXd::Identity(grad_u_q.rows(), grad_u_q.cols()) + grad_u_q;
					Eigen::MatrixXd FmT = def_grad.inverse().transpose();
					stress = mu(q) * (def_grad - FmT) + lambda(q) * std::log(def_grad.determinant()) * FmT;
				}
				else
					log_and_throw_error("Unknown formulation!");
				val(q) = pow(stress.squaredNorm(), power / 2.);
			}
		});

		j.set_dj_dgradu([formulation, power](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
			val.setZero(grad_u.rows(), grad_u.cols());
			const int dim = sqrt(grad_u.cols());
			const int actual_dim = (formulation == "Laplacian") ? 1 : dim;
			Eigen::MatrixXd grad_u_q, stress, stress_dstress;
			for (int q = 0; q < grad_u.rows(); q++)
			{
				if (formulation == "Laplacian")
				{
					stress = grad_u.row(q);
					stress_dstress = 2 * stress;
				}
				else if (formulation == "LinearElasticity")
				{
					vector2matrix(grad_u.row(q), grad_u_q);
					stress = mu(q) * (grad_u_q + grad_u_q.transpose()) + lambda(q) * grad_u_q.trace() * Eigen::MatrixXd::Identity(grad_u_q.rows(), grad_u_q.cols());
					stress_dstress = mu(q) * (stress + stress.transpose()) + lambda(q) * stress.trace() * Eigen::MatrixXd::Identity(stress.rows(), stress.cols());
				}
				else if (formulation == "NeoHookean")
				{
					vector2matrix(grad_u.row(q), grad_u_q);
					Eigen::MatrixXd def_grad = Eigen::MatrixXd::Identity(grad_u_q.rows(), grad_u_q.cols()) + grad_u_q;
					Eigen::MatrixXd FmT = def_grad.inverse().transpose();
					stress = mu(q) * (def_grad - FmT) + lambda(q) * std::log(def_grad.determinant()) * FmT;
					stress_dstress = mu(q) * stress + FmT * stress.transpose() * FmT * (mu(q) - lambda(q) * std::log(def_grad.determinant())) + (lambda(q) * (FmT.array() * stress.array()).sum()) * FmT;
				}
				else
					log_and_throw_error("Unknown formulation!");

				const double coef = power * pow(stress.squaredNorm(), power / 2. - 1.);
				for (int i = 0; i < actual_dim; i++)
					for (int l = 0; l < dim; l++)
						val(q, i * dim + l) = coef * stress_dstress(i, l);
			}
		});

		return j;
	}

  	void StressForm::compute_partial_gradient_unweighted(const Eigen::VectorXd &x, Eigen::VectorXd &gradv) const
  {
    	SpatialIntegralForm::compute_partial_gradient_unweighted(x, gradv);
		for (const auto &param_map : variable_to_simulations_)
		{
			const auto &parametrization = param_map->get_parametrization();
			const auto &state = param_map->get_state();
			const auto &param_type = param_map->get_parameter_type();

			if (&state != &state_)
				continue;

			Eigen::VectorXd term;
			if (param_type == ParameterType::Material)
				log_and_throw_error("Doesn't support stress derivative wrt. material!");

			if (term.size() > 0)
				gradv += parametrization.apply_jacobian(term, x);
		}
  }

	IntegrableFunctional PositionForm::get_integral_functional() const
	{
		IntegrableFunctional j;
		const int dim = dim_;

		j.set_j([dim](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
			val = u.col(dim) + pts.col(dim);
		});

		j.set_dj_du([dim](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
			val.setZero(u.rows(), u.cols());
			val.col(dim).setOnes();
		});

		j.set_dj_dx([dim](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
			val.setZero(pts.rows(), pts.cols());
			val.col(dim).setOnes();
		});

		return j;
	}

	IntegrableFunctional TargetForm::get_integral_functional() const
	{
		IntegrableFunctional j;
		if (target_state_)
		{
			assert(target_state_->diff_cached.size() > 0);

			auto j_func = [this](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
				val.setZero(u.rows(), 1);
				const int e = params["elem"];
				int e_ref;
				if (auto search = e_to_ref_e_.find(e); search != e_to_ref_e_.end())
					e_ref = search->second;
				else
					e_ref = e;
				const auto &gbase_ref = target_state_->geom_bases()[e_ref];

				Eigen::MatrixXd pts_ref;
				gbase_ref.eval_geom_mapping(local_pts, pts_ref);

				Eigen::MatrixXd u_ref, grad_u_ref;
				const Eigen::MatrixXd &sol_ref = target_state_->problem->is_time_dependent() ? target_state_->diff_cached[params["step"].get<int>()].u : target_state_->diff_cached[0].u;
				io::Evaluator::interpolate_at_local_vals(*(target_state_->mesh), target_state_->problem->is_scalar(), target_state_->bases, target_state_->geom_bases(), e_ref, local_pts, sol_ref, u_ref, grad_u_ref);

				for (int q = 0; q < u.rows(); q++)
				{
					val(q) = ((u_ref.row(q) + pts_ref.row(q)) - (u.row(q) + pts.row(q))).squaredNorm();
				}
			};

			auto djdu_func = [this](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
				val.setZero(u.rows(), u.cols());
				const int e = params["elem"];
				int e_ref;
				if (auto search = e_to_ref_e_.find(e); search != e_to_ref_e_.end())
					e_ref = search->second;
				else
					e_ref = e;
				const auto &gbase_ref = target_state_->geom_bases()[e_ref];

				Eigen::MatrixXd pts_ref;
				gbase_ref.eval_geom_mapping(local_pts, pts_ref);

				Eigen::MatrixXd u_ref, grad_u_ref;
				const Eigen::MatrixXd &sol_ref = target_state_->problem->is_time_dependent() ? target_state_->diff_cached[params["step"].get<int>()].u : target_state_->diff_cached[0].u;
				io::Evaluator::interpolate_at_local_vals(*(target_state_->mesh), target_state_->problem->is_scalar(), target_state_->bases, target_state_->geom_bases(), e_ref, local_pts, sol_ref, u_ref, grad_u_ref);

				for (int q = 0; q < u.rows(); q++)
				{
					auto x = (u.row(q) + pts.row(q)) - (u_ref.row(q) + pts_ref.row(q));
					val.row(q) = 2 * x;
				}
			};

			j.set_j(j_func);
			j.set_dj_du(djdu_func);
			j.set_dj_dx(djdu_func); // only used for shape derivative
		}
		else if (have_target_func)
		{
			auto j_func = [this](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
				val.setZero(u.rows(), 1);
				for (int q = 0; q < u.rows(); q++)
				{
					Eigen::VectorXd x = u.row(q) + pts.row(q);
					val(q) = target_func(x(0), x(1), x.size() == 2 ? 0 : x(2), 0, params["elem"]);
				}
			};

			auto djdu_func = [this](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
				val.setZero(u.rows(), u.cols());
				for (int q = 0; q < u.rows(); q++)
				{
					Eigen::VectorXd x = u.row(q) + pts.row(q);
					for (int d = 0; d < val.cols(); d++)
						val(q, d) = target_func_grad[d](x(0), x(1), x.size() == 2 ? 0 : x(2), 0, params["elem"]);
				}
			};

			j.set_j(j_func);
			j.set_dj_du(djdu_func);
			j.set_dj_dx(djdu_func); // only used for shape derivative
		}
		else // error wrt. a constant displacement
		{
			if (target_disp.size() == state_.mesh->dimension())
			{
				auto j_func = [this](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
					val.setZero(u.rows(), 1);

					for (int q = 0; q < u.rows(); q++)
					{
						Eigen::VectorXd err = u.row(q) - this->target_disp.transpose();
						for (int d = 0; d < active_dimension_mask.size(); d++)
							if (!active_dimension_mask[d])
								err(d) = 0;
						val(q) = err.squaredNorm();
					}
				};
				auto djdu_func = [this](const Eigen::MatrixXd &local_pts, const Eigen::MatrixXd &pts, const Eigen::MatrixXd &u, const Eigen::MatrixXd &grad_u, const Eigen::MatrixXd &lambda, const Eigen::MatrixXd &mu, const json &params, Eigen::MatrixXd &val) {
					val.setZero(u.rows(), u.cols());

					for (int q = 0; q < u.rows(); q++)
					{
						Eigen::VectorXd err = u.row(q) - this->target_disp.transpose();
						for (int d = 0; d < active_dimension_mask.size(); d++)
							if (!active_dimension_mask[d])
								err(d) = 0;
						val.row(q) = 2 * err;
					}
				};

				j.set_j(j_func);
				j.set_dj_du(djdu_func);
			}
			else
				log_and_throw_error("Only constant target displacement is supported!");
		}

		return j;
	}

	void TargetForm::set_reference(const std::shared_ptr<const State> &target_state, const std::set<int> &reference_cached_body_ids)
	{
		target_state_ = target_state;

		std::map<int, std::vector<int>> ref_interested_body_id_to_e;
		int ref_count = 0;
		for (int e = 0; e < target_state_->bases.size(); ++e)
		{
			int body_id = target_state_->mesh->get_body_id(e);
			if (reference_cached_body_ids.size() > 0 && reference_cached_body_ids.count(body_id) == 0)
				continue;
			if (ref_interested_body_id_to_e.find(body_id) != ref_interested_body_id_to_e.end())
				ref_interested_body_id_to_e[body_id].push_back(e);
			else
				ref_interested_body_id_to_e[body_id] = {e};
			ref_count++;
		}

		std::map<int, std::vector<int>> interested_body_id_to_e;
		int count = 0;
		for (int e = 0; e < state_.bases.size(); ++e)
		{
			int body_id = state_.mesh->get_body_id(e);
			if (reference_cached_body_ids.size() > 0 && reference_cached_body_ids.count(body_id) == 0)
				continue;
			if (interested_body_id_to_e.find(body_id) != interested_body_id_to_e.end())
				interested_body_id_to_e[body_id].push_back(e);
			else
				interested_body_id_to_e[body_id] = {e};
			count++;
		}

		if (count != ref_count)
			logger().error("Number of interested elements in the reference and optimization examples do not match! {} {}", count, ref_count);
		else
			logger().trace("Found {} matching elements.", count);

		for (const auto &kv : interested_body_id_to_e)
		{
			for (int i = 0; i < kv.second.size(); ++i)
			{
				e_to_ref_e_[kv.second[i]] = ref_interested_body_id_to_e[kv.first][i];
			}
		}
	}

	void TargetForm::set_reference(const json &func, const json &grad_func)
	{
		target_func.init(func);
		for (size_t k = 0; k < grad_func.size(); k++)
			target_func_grad[k].init(grad_func[k]);
		have_target_func = true;
	}

} // namespace polyfem::solver