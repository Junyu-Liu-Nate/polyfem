#include "OptimizationProblem.hpp"
#include "Optimizations.hpp"
#include "ShapeProblem.hpp"
#include <polyfem/assembler/RhsAssembler.hpp>
#include "TopologyOptimizationProblem.hpp"
#include "MaterialProblem.hpp"
#include "InitialConditionProblem.hpp"
#include "LBFGSBSolver.hpp"
#include "LBFGSSolver.hpp"
#include "BFGSSolver.hpp"
#include "MMASolver.hpp"
#include "GradientDescentSolver.hpp"
#include <polyfem/utils/SplineParam.hpp>

#include <map>

namespace polyfem
{
	double cross2(double x, double y)
	{
		x = abs(x);
		y = abs(y);
		if (x > y)
			std::swap(x, y);

		if (x < 0.1)
			return 0.05;
		return 0.95;
	}

	double cross3(double x, double y, double z)
	{
		x = abs(x);
		y = abs(y);
		z = abs(z);
		if (x > y)
			std::swap(x, y);
		if (y > z)
			std::swap(y, z);
		if (x > y)
			std::swap(x, y);

		if (y < 0.2)
			return 0.001;
		return 1;
	}

	template <typename ProblemType>
	std::shared_ptr<cppoptlib::NonlinearSolver<ProblemType>> make_nl_solver(const json &solver_params)
	{
		const std::string name = solver_params.contains("solver") ? solver_params["solver"].template get<std::string>() : "lbfgs";
		if (name == "GradientDescent" || name == "gradientdescent" || name == "gradient")
		{
			return std::make_shared<cppoptlib::GradientDescentSolver<ProblemType>>(
				solver_params);
		}
		else if (name == "lbfgs" || name == "LBFGS" || name == "L-BFGS")
		{
			return std::make_shared<cppoptlib::LBFGSSolver<ProblemType>>(
				solver_params);
		}
		else if (name == "bfgs" || name == "BFGS" || name == "BFGS")
		{
			return std::make_shared<cppoptlib::BFGSSolver<ProblemType>>(
				solver_params);
		}
		else if (name == "lbfgsb" || name == "LBFGSB" || name == "L-BFGS-B")
		{
			return std::make_shared<cppoptlib::LBFGSBSolver<ProblemType>>(
				solver_params);
		}
		else if (name == "mma" || name == "MMA")
		{
			return std::make_shared<cppoptlib::MMASolver<ProblemType>>(
				solver_params);
		}
		else
		{
			throw std::invalid_argument(fmt::format("invalid nonlinear solver type: {}", name));
		}
	}

	double matrix_dot(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) { return (A.array() * B.array()).sum(); }

	void initial_condition_optimization(State &state, const std::shared_ptr<CompositeFunctional> j)
	{
		const auto &opt_params = state.args["optimization"];
		const auto &opt_nl_params = state.args["solver"]["optimization_nonlinear"];

		std::shared_ptr<InitialConditionProblem> initial_problem = std::make_shared<InitialConditionProblem>(state, j);
		std::shared_ptr<cppoptlib::NonlinearSolver<InitialConditionProblem>> nlsolver = make_nl_solver<InitialConditionProblem>(opt_nl_params);
		nlsolver->setLineSearch(opt_nl_params["line_search"]["method"]);

		json initial_params;
		for (const auto &param : opt_params["parameters"])
		{
			if (param["type"] == "initial")
			{
				initial_params = param;
				break;
			}
		}

		// fix certain object
		std::set<int> optimize_body_ids;
		if (initial_params.contains("volume_selection"))
		{
			for (int i : initial_params["volume_selection"])
				optimize_body_ids.insert(i);
		}
		else
			logger().info("No optimization body specified, optimize initial condition of every mesh...");

		// to get the initial velocity
		{
			const auto &gbases = state.iso_parametric() ? state.bases : state.geom_bases;
			json rhs_solver_params = state.args["solver"]["linear"];
			if (!rhs_solver_params.contains("Pardiso"))
				rhs_solver_params["Pardiso"] = {};
			rhs_solver_params["Pardiso"]["mtype"] = -2; // matrix type for Pardiso (2 = SPD)

			state.step_data.rhs_assembler = std::make_shared<assembler::RhsAssembler>(
				state.assembler, *state.mesh, state.obstacle, state.input_dirichelt,
				state.n_bases, state.problem->is_scalar() ? 1 : state.mesh->dimension(),
				state.bases, gbases, state.ass_vals_cache,
				state.formulation(), *state.problem,
				state.args["space"]["advanced"]["bc_method"],
				state.args["solver"]["linear"]["solver"], state.args["solver"]["linear"]["precond"], rhs_solver_params);

			assembler::RhsAssembler &rhs_assembler = *state.step_data.rhs_assembler;
			rhs_assembler.initial_solution(state.initial_sol_update);
			rhs_assembler.initial_velocity(state.initial_vel_update);
		}
		auto initial_guess_vel = state.initial_vel_update, initial_guess_sol = state.initial_sol_update;

		const int dim = state.mesh->dimension();
		const int dof = state.n_bases;

		std::map<int, std::array<int, 2>> body_id_map; // from body_id to {node_id, index}
		int n = 0;
		for (int e = 0; e < state.bases.size(); e++)
		{
			const int body_id = state.mesh->get_body_id(e);
			if (!body_id_map.count(body_id) && (optimize_body_ids.count(body_id) || optimize_body_ids.size() == 0))
			{
				body_id_map[body_id] = {{state.bases[e].bases[0].global()[0].index, n}};
				n++;
			}
		}
		logger().info("{} objects found, each object has a constant initial velocity and position...", body_id_map.size());

		// by default optimize for initial velocity
		if (!initial_params.contains("restriction"))
			initial_params["restriction"] = "velocity";

		if (initial_params["restriction"].get<std::string>() == "velocity")
		{
			initial_problem->x_to_param = [&](const InitialConditionProblem::TVector &x, Eigen::MatrixXd &init_sol, Eigen::MatrixXd &init_vel) {
				init_sol = initial_guess_sol;
				init_vel = initial_guess_vel;
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					if (!body_id_map.count(body_id))
						continue;
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
							for (int d = 0; d < dim; d++)
								init_vel(g.index * dim + d) = x(body_id_map[body_id][1] * dim + d);
				}
				std::cout << "initial velocity: " << std::setprecision(16) << x.transpose() << "\n";
			};

			initial_problem->param_to_x = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size());
				for (auto i : body_id_map)
					for (int d = 0; d < dim; d++)
						x(i.second[1] * dim + d) = init_vel(i.second[0] * dim + d);
			};

			initial_problem->dparam_to_dx = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size());
				std::vector<bool> visited(dof, false);
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					if (!body_id_map.count(body_id))
						continue;
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
						{
							if (!visited[g.index])
								visited[g.index] = true;
							else
								continue;
							for (int d = 0; d < dim; d++)
								x(dim * body_id_map[body_id][1] + d) += init_vel(g.index * dim + d);
						}
				}
			};
		}
		else if (initial_params["restriction"].get<std::string>() == "velocityXZ")
		{
			initial_problem->x_to_param = [&](const InitialConditionProblem::TVector &x, Eigen::MatrixXd &init_sol, Eigen::MatrixXd &init_vel) {
				init_sol = initial_guess_sol;
				init_vel = initial_guess_vel;
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					if (!body_id_map.count(body_id))
						continue;
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
							for (int d = 0; d < dim; d++)
								init_vel(g.index * dim + d) = x(body_id_map[body_id][1] * dim + d);
				}
				logger().debug("initial velocity: {}", x.transpose());
			};

			initial_problem->param_to_x = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size());
				for (auto i : body_id_map)
					for (int d = 0; d < dim; d++)
						x(i.second[1] * dim + d) = init_vel(i.second[0] * dim + d);
				logger().debug("initial velocity: {}", x.transpose());
			};

			initial_problem->dparam_to_dx = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size());
				std::vector<bool> visited(dof, false);
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					if (!body_id_map.count(body_id))
						continue;
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
						{
							if (!visited[g.index])
								visited[g.index] = true;
							else
								continue;
							for (int d = 0; d < dim; d++)
							{
								if (d == 1)
									continue;
								x(dim * body_id_map[body_id][1] + d) += init_vel(g.index * dim + d);
							}
						}
				}
			};
		}
		else if (initial_params["restriction"].get<std::string>() == "position")
		{
			initial_problem->x_to_param = [&](const InitialConditionProblem::TVector &x, Eigen::MatrixXd &init_sol, Eigen::MatrixXd &init_vel) {
				init_sol.setZero(dof * dim, 1);
				init_vel.setZero(dof * dim, 1);
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					if (!body_id_map.count(body_id))
						continue;
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
							for (int d = 0; d < dim; d++)
								init_sol(g.index * dim + d) = x(body_id_map[body_id][1] * dim + d);
				}
				logger().debug("initial solution: {}", x.transpose());
			};

			initial_problem->param_to_x = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size());
				for (auto i : body_id_map)
				{
					for (int d = 0; d < dim; d++)
						x(i.second[1] * dim + d) = init_sol(i.second[0] * dim + d);
				}
				logger().debug("initial solution: {}", x.transpose());
			};

			initial_problem->dparam_to_dx = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size());
				std::vector<bool> visited(dof, false);
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					if (!body_id_map.count(body_id))
						continue;
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
						{
							if (!visited[g.index])
								visited[g.index] = true;
							else
								continue;
							for (int d = 0; d < dim; d++)
								x(dim * body_id_map[body_id][1] + d) += init_sol(g.index * dim + d);
						}
				}
			};
		}
		else
		{
			initial_problem->x_to_param = [&](const InitialConditionProblem::TVector &x, Eigen::MatrixXd &init_sol, Eigen::MatrixXd &init_vel) {
				init_sol = initial_guess_sol;
				init_vel = initial_guess_vel;
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
							for (int d = 0; d < dim; d++)
							{
								init_sol(g.index * dim + d) = x(body_id_map[body_id][1] * dim + d);
								init_vel(g.index * dim + d) = x(body_id_map[body_id][1] * dim + d + dim * body_id_map.size());
							}
				}
				logger().debug("initial velocity: {}", x.tail(x.size() / 2).transpose());
				logger().debug("initial position: {}", x.head(x.size() / 2).transpose());
			};

			initial_problem->param_to_x = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size() * 2);
				for (auto i : body_id_map)
					for (int d = 0; d < dim; d++)
					{
						x(i.second[1] * dim + d) = init_sol(i.second[0] * dim + d);
						x(i.second[1] * dim + d + dim * body_id_map.size()) = init_vel(i.second[0] * dim + d);
					}
				logger().debug("initial velocity: {}", x.tail(x.size() / 2).transpose());
				logger().debug("initial position: {}", x.head(x.size() / 2).transpose());
			};

			initial_problem->dparam_to_dx = [&](InitialConditionProblem::TVector &x, const Eigen::MatrixXd &init_sol, const Eigen::MatrixXd &init_vel) {
				x.setZero(dim * body_id_map.size() * 2);
				std::vector<bool> visited(dof, false);
				for (int e = 0; e < state.bases.size(); e++)
				{
					const int body_id = state.mesh->get_body_id(e);
					for (auto &bs : state.bases[e].bases)
						for (auto &g : bs.global())
						{
							if (!visited[g.index])
								visited[g.index] = true;
							else
								continue;
							for (int d = 0; d < dim; d++)
							{
								x(dim * body_id_map[body_id][1] + d) += init_sol(g.index * dim + d);
								x(dim * body_id_map[body_id][1] + d + dim * body_id_map.size()) += init_vel(g.index * dim + d);
							}
						}
				}
			};
		}
		Eigen::VectorXd x;
		initial_problem->param_to_x(x, state.initial_sol_update, state.initial_vel_update);

		nlsolver->minimize(*initial_problem, x);

		json solver_info;
		nlsolver->getInfo(solver_info);
		std::cout << solver_info << std::endl;
	}

	void material_optimization(State &state, const std::shared_ptr<CompositeFunctional> j)
	{
		const auto &opt_params = state.args["optimization"];
		const auto &opt_nl_params = state.args["solver"]["optimization_nonlinear"];
		json material_params;

		state.args["output"]["paraview"]["options"]["material"] = true;

		for (const auto &param : opt_params["parameters"])
		{
			if (param["type"] == "material")
			{
				material_params = param;
				break;
			}
		}

		std::shared_ptr<MaterialProblem> material_problem = std::make_shared<MaterialProblem>(state, j);
		std::shared_ptr<cppoptlib::NonlinearSolver<MaterialProblem>> nlsolver = make_nl_solver<MaterialProblem>(opt_nl_params);
		nlsolver->setLineSearch(opt_nl_params["line_search"]["method"]);

		// fix certain object
		std::set<int> optimize_body_ids;
		if (material_params.contains("volume_selection"))
		{
			for (int i : material_params["volume_selection"])
				optimize_body_ids.insert(i);
		}
		else
			logger().info("No optimization body specified, optimize material of every mesh...");

		const int dim = state.mesh->dimension();
		const int dof = state.bases.size();

		std::map<int, std::array<int, 2>> body_id_map; // from body_id to {elem_id, index}
		int n = 0;
		for (int e = 0; e < dof; e++)
		{
			const int body_id = state.mesh->get_body_id(e);
			if (!body_id_map.count(body_id) && (optimize_body_ids.count(body_id) || optimize_body_ids.size() == 0))
			{
				body_id_map[body_id] = {{e, n}};
				n++;
			}
		}

		// constraints on optimization
		if (material_params.contains("restriction"))
		{
			if (material_params["restriction"].get<std::string>() == "constant")
			{
				logger().info("{} objects found, each object has constant material parameters...", body_id_map.size());

				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					auto cur_lambdas = state.assembler.lame_params().lambda_mat_;
					auto cur_mus = state.assembler.lame_params().mu_mat_;

					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);

						cur_lambdas(e) = x(body_id_map[body_id][1] * 2 + 0);
						cur_mus(e) = x(body_id_map[body_id][1] * 2 + 1);
					}
					state.assembler.update_lame_params(cur_lambdas, cur_mus);
				};
				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;
					x.setZero(2 * body_id_map.size());
					for (auto i : body_id_map)
					{
						x(i.second[1] * 2 + 0) = cur_lambdas(i.second[0]);
						x(i.second[1] * 2 + 1) = cur_mus(i.second[0]);
					}
					logger().debug("material: {}", x.transpose());
				};
				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					dx.setZero(2 * body_id_map.size());
					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);
						dx(body_id_map[body_id][1] * 2 + 0) += dparams(e);
						dx(body_id_map[body_id][1] * 2 + 1) += dparams(e + dof);
					}
				};
			}
			else if (material_params["restriction"].get<std::string>() == "log")
			{
				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					auto cur_lambdas = state.assembler.lame_params().lambda_mat_;
					auto cur_mus = state.assembler.lame_params().mu_mat_;
					cur_mus = x.tail(dof).array().exp().matrix();
					cur_lambdas = x.head(dof).array().exp().matrix();
					state.assembler.update_lame_params(cur_lambdas, cur_mus);
				};

				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					x.resize(2 * dof);
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;

					x.tail(dof) = cur_mus.array().log().matrix();
					x.head(dof) = cur_lambdas.array().log().matrix();
				};

				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					dx.resize(dof * 2);
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;

					dx = dparams.head(2 * dof);
					dx.tail(dof) *= cur_mus;
					dx.head(dof) *= cur_lambdas;
				};
			}
			else if (material_params["restriction"].get<std::string>() == "constant_log")
			{
				logger().info("{} objects found, each object has constant material parameters...", body_id_map.size());

				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					auto cur_lambdas = state.assembler.lame_params().lambda_mat_;
					auto cur_mus = state.assembler.lame_params().mu_mat_;

					cur_lambdas.setConstant(dof, 1, std::exp(x(0)));
					cur_mus.setConstant(dof, 1, std::exp(x(1)));

					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);
						if (!body_id_map.count(body_id))
							continue;
						cur_lambdas(e) = std::exp(x(body_id_map[body_id][1] * 2 + 0));
						cur_mus(e) = std::exp(x(body_id_map[body_id][1] * 2 + 1));
					}
					Eigen::VectorXd x_display = Eigen::VectorXd::Zero(x.size());
					for (int i = 0; i < x.size(); i++)
						x_display(i) = std::exp(x(i));
					state.assembler.update_lame_params(cur_lambdas, cur_mus);
					logger().debug("material: {}", x_display.transpose());
				};
				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;
					x.setZero(2 * body_id_map.size());
					for (auto i : body_id_map)
					{
						x(i.second[1] * 2 + 0) = std::log(cur_lambdas(i.second[0]));
						x(i.second[1] * 2 + 1) = std::log(cur_mus(i.second[0]));
					}
					Eigen::VectorXd x_display = Eigen::VectorXd::Zero(x.size());
					for (int i = 0; i < x.size(); i++)
						x_display(i) = std::exp(x(i));
					logger().debug("material: {}", x_display.transpose());
				};
				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;
					dx.setZero(2 * body_id_map.size());
					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);
						if (!body_id_map.count(body_id))
							continue;
						dx(body_id_map[body_id][1] * 2 + 0) += dparams(e) * cur_lambdas(e);
						dx(body_id_map[body_id][1] * 2 + 1) += dparams(e + dof) * cur_mus(e);
					}
				};
			}
			else if (material_params["restriction"].get<std::string>() == "constant_log_friction")
			{
				logger().info("{} objects found, each object has constant material parameters...", body_id_map.size());

				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					auto cur_lambdas = state.assembler.lame_params().lambda_mat_;
					auto cur_mus = state.assembler.lame_params().mu_mat_;

					cur_lambdas.setConstant(dof, 1, std::exp(x(0)));
					cur_mus.setConstant(dof, 1, std::exp(x(1)));

					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);

						cur_lambdas(e) = std::exp(x(body_id_map[body_id][1] * 2 + 0));
						cur_mus(e) = std::exp(x(body_id_map[body_id][1] * 2 + 1));
					}
					state.assembler.update_lame_params(cur_lambdas, cur_mus);

					state.args["contact"]["friction_coefficient"] = std::exp(x(x.size() - 1));

					Eigen::VectorXd x_display = Eigen::VectorXd::Zero(x.size() - 1);
					for (int i = 0; i < x_display.size(); i++)
						x_display(i) = std::exp(x(i));
					logger().debug("material: {}", x_display.transpose());
					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
				};
				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;
					x.setZero(2 * body_id_map.size() + 1);
					for (auto i : body_id_map)
					{
						x(i.second[1] * 2 + 0) = std::log(cur_lambdas(i.second[0]));
						x(i.second[1] * 2 + 1) = std::log(cur_mus(i.second[0]));
					}
					x(x.size() - 1) = std::log(state.args["contact"]["friction_coefficient"].get<double>());

					Eigen::VectorXd x_display = Eigen::VectorXd::Zero(x.size() - 1);
					for (int i = 0; i < x_display.size(); i++)
						x_display(i) = std::exp(x(i));
					logger().debug("material: {}", x_display.transpose());
					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
				};
				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;
					dx.setZero(2 * body_id_map.size() + 1);
					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);
						dx(body_id_map[body_id][1] * 2 + 0) += dparams(e) * cur_lambdas(e);
						dx(body_id_map[body_id][1] * 2 + 1) += dparams(e + dof) * cur_mus(e);
					}
					dx(dx.size() - 1) = dparams(2 * dof) * state.args["contact"]["friction_coefficient"].get<double>();
				};
			}
			else if (material_params["restriction"].get<std::string>() == "constant_log_friction_damping")
			{
				logger().info("{} objects found, each object has constant material parameters...", body_id_map.size());

				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					auto cur_lambdas = state.assembler.lame_params().lambda_mat_;
					auto cur_mus = state.assembler.lame_params().mu_mat_;

					cur_lambdas.setConstant(dof, 1, std::exp(x(0)));
					cur_mus.setConstant(dof, 1, std::exp(x(1)));

					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);

						cur_lambdas(e) = std::exp(x(body_id_map[body_id][1] * 2 + 0));
						cur_mus(e) = std::exp(x(body_id_map[body_id][1] * 2 + 1));
					}
					state.assembler.update_lame_params(cur_lambdas, cur_mus);

					state.args["contact"]["friction_coefficient"] = std::exp(x(x.size() - 3));
					state.damping_assembler.local_assembler().set_params(std::exp(x(x.size() - 2)), std::exp(x(x.size() - 1)));

					Eigen::VectorXd x_display = Eigen::VectorXd::Zero(x.size() - 3);
					for (int i = 0; i < x_display.size(); i++)
						x_display(i) = std::exp(x(i));
					logger().debug("material: {}", x_display.transpose());
					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
					logger().debug("psi = {}, phi = {}", state.damping_assembler.local_assembler().get_psi(), state.damping_assembler.local_assembler().get_phi());
				};
				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;
					x.setZero(2 * body_id_map.size() + 3);
					for (auto i : body_id_map)
					{
						x(i.second[1] * 2 + 0) = std::log(cur_lambdas(i.second[0]));
						x(i.second[1] * 2 + 1) = std::log(cur_mus(i.second[0]));
					}
					x(x.size() - 3) = std::log(state.args["contact"]["friction_coefficient"].get<double>());
					x(x.size() - 2) = std::log(state.damping_assembler.local_assembler().get_psi());
					x(x.size() - 1) = std::log(state.damping_assembler.local_assembler().get_phi());

					Eigen::VectorXd x_display = Eigen::VectorXd::Zero(x.size() - 3);
					for (int i = 0; i < x_display.size(); i++)
						x_display(i) = std::exp(x(i));
					logger().debug("material: {}", x_display.transpose());
					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
					logger().debug("psi = {}, phi = {}", state.damping_assembler.local_assembler().get_psi(), state.damping_assembler.local_assembler().get_phi());
				};
				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					const auto &cur_lambdas = state.assembler.lame_params().lambda_mat_;
					const auto &cur_mus = state.assembler.lame_params().mu_mat_;
					dx.setZero(2 * body_id_map.size() + 3);
					for (int e = 0; e < dof; e++)
					{
						const int body_id = state.mesh->get_body_id(e);
						dx(body_id_map[body_id][1] * 2 + 0) += dparams(e) * cur_lambdas(e);
						dx(body_id_map[body_id][1] * 2 + 1) += dparams(e + dof) * cur_mus(e);
					}
					dx.tail(3) = dparams.tail(3);
					dx(dx.size() - 3) = dparams(dparams.size() - 3) * state.args["contact"]["friction_coefficient"].get<double>();
					dx(dx.size() - 2) = dparams(dparams.size() - 2) * state.damping_assembler.local_assembler().get_psi();
					dx(dx.size() - 1) = dparams(dparams.size() - 1) * state.damping_assembler.local_assembler().get_phi();
				};
			}
			else if (material_params["restriction"].get<std::string>() == "friction_damping")
			{
				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					state.args["contact"]["friction_coefficient"] = x(x.size() - 3);
					state.damping_assembler.local_assembler().set_params(x(x.size() - 2), x(x.size() - 1));

					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
					logger().debug("psi = {}, phi = {}", state.damping_assembler.local_assembler().get_psi(), state.damping_assembler.local_assembler().get_phi());
				};
				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					x.setZero(3);

					x(x.size() - 3) = state.args["contact"]["friction_coefficient"].get<double>();
					x(x.size() - 2) = state.damping_assembler.local_assembler().get_psi();
					x(x.size() - 1) = state.damping_assembler.local_assembler().get_phi();

					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
					logger().debug("psi = {}, phi = {}", state.damping_assembler.local_assembler().get_psi(), state.damping_assembler.local_assembler().get_phi());
				};
				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					dx.setZero(3);

					dx.tail(3) = dparams.tail(3);

					logger().debug("opt grad: {}", dx.transpose());
				};
			}
			else if (material_params["restriction"].get<std::string>() == "friction")
			{
				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					state.args["contact"]["friction_coefficient"] = x(0);

					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
				};
				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					x.setZero(1);

					x(0) = state.args["contact"]["friction_coefficient"].get<double>();

					logger().debug("friction coeff = {}", state.args["contact"]["friction_coefficient"].get<double>());
				};
				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					dx.setZero(1);

					dx(0) = dparams(dparams.size() - 3);

					logger().debug("opt grad: {}", dx.transpose());
				};
			}
			else if (material_params["restriction"].get<std::string>() == "damping")
			{
				material_problem->x_to_param = [&](const MaterialProblem::TVector &x, State &state) {
					state.damping_assembler.local_assembler().set_params(x(x.size() - 2), x(x.size() - 1));
					logger().debug("psi = {}, phi = {}", state.damping_assembler.local_assembler().get_psi(), state.damping_assembler.local_assembler().get_phi());
				};
				material_problem->param_to_x = [&](MaterialProblem::TVector &x, State &state) {
					x.setZero(2);

					x(x.size() - 2) = state.damping_assembler.local_assembler().get_psi();
					x(x.size() - 1) = state.damping_assembler.local_assembler().get_phi();

					logger().debug("psi = {}, phi = {}", state.damping_assembler.local_assembler().get_psi(), state.damping_assembler.local_assembler().get_phi());
				};
				material_problem->dparam_to_dx = [&](MaterialProblem::TVector &dx, const Eigen::VectorXd &dparams, State &state) {
					dx.setZero(2);

					dx.tail(2) = dparams.tail(2);

					logger().debug("opt grad: {}", dx.transpose());
				};
			}
		}

		Eigen::VectorXd x;
		material_problem->param_to_x(x, state);

		nlsolver->minimize(*material_problem, x);

		json solver_info;
		nlsolver->getInfo(solver_info);
		std::cout << solver_info << std::endl;
	}

	void shape_optimization(State &state, const std::shared_ptr<CompositeFunctional> j)
	{
		const auto &opt_params = state.args["optimization"];
		const auto &opt_nl_params = state.args["solver"]["optimization_nonlinear"];

		std::shared_ptr<ShapeProblem> shape_problem = std::make_shared<ShapeProblem>(state, j);
		std::shared_ptr<cppoptlib::NonlinearSolver<ShapeProblem>> nlsolver = make_nl_solver<ShapeProblem>(opt_nl_params);
		nlsolver->setLineSearch(opt_nl_params["line_search"]["method"]);

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		state.get_vf(V, F);

		const auto &shape_params = shape_problem->get_shape_params();
		if (shape_params.contains("restriction"))
		{
			if (shape_params["restriction"].get<std::string>() == "cubic_hermite_spline" and shape_params.contains("spline_specification"))
			{
				// Assume there is one spline with id 10.
				const auto &spline_params = shape_params["spline_specification"];
				const int boundary_id = spline_params[0]["id"].get<int>();
				auto control_points = spline_params[0]["control_point"];
				auto tangents = spline_params[0]["tangent"];
				const int sampling = spline_params[0]["sampling"].get<int>();
				std::map<int, Eigen::MatrixXd> control_point, tangent;
				Eigen::MatrixXd c(2, 2), t(2, 2);
				for (int i = 0; i < 2; ++i)
				{
					for (int j = 0; j < 2; ++j)
					{
						c(i, j) = control_points[i][j];
						t(i, j) = tangents[i][j];
					}
				}
				logger().trace("Given tangents are: {}", t);
				control_point = {{boundary_id, c}};
				tangent = {{boundary_id, t}};
				SplineParam spline_param(control_point, tangent, shape_problem->optimization_boundary_to_node, V, sampling);
				shape_problem->param_to_x = [spline_param](ShapeProblem::TVector &x, const Eigen::MatrixXd &V) {
					std::map<int, Eigen::MatrixXd> control_point, tangent;
					spline_param.get_parameters(V, control_point, tangent);
					x.setZero(2 * tangent.size() + 2);
					int index = 0;
					int last_id = -1;
					for (const auto &kv : tangent)
					{
						x.segment(index, 2) = kv.second.row(0);
						index += 2;
						last_id = kv.first;
					}
					x.segment(index, 2) = tangent.at(last_id).row(1);
				};
				shape_problem->x_to_param = [control_point, tangent, spline_param](const ShapeProblem::TVector &x, const Eigen::MatrixXd &V_prev, Eigen::MatrixXd &V_) {
					std::map<int, Eigen::MatrixXd> new_tangent;
					int index = 0;
					for (const auto &kv : tangent)
					{
						Eigen::MatrixXd tangent_matrix(2, 2);
						tangent_matrix.row(0) = x.segment(index, 2);
						tangent_matrix.row(1) = x.segment(index + 2, 2);
						new_tangent[kv.first] = tangent_matrix;
						index += 2;
					}
					spline_param.reparametrize(control_point, new_tangent, V_prev, V_);
				};
				shape_problem->dparam_to_dx = [tangent, spline_param](ShapeProblem::TVector &grad_x, const ShapeProblem::TVector &grad_v) {
					grad_x.setZero(2 * tangent.size() + 2);
					int index = 0;
					for (const auto &kv : tangent)
					{
						Eigen::VectorXd grad_control_point, grad_tangent;
						spline_param.derivative_wrt_params(grad_v, kv.first, grad_control_point, grad_tangent);
						grad_x.segment(index, 4) += 0.5 * grad_tangent;
						index += 2;
					}
					grad_x.segment(0, 2) *= 2;
					grad_x.segment(index, 2) *= 2;
				};
			}
		}

		Eigen::VectorXd x;
		shape_problem->param_to_x(x, V);
		nlsolver->minimize(*shape_problem, x);

		json solver_info;
		nlsolver->getInfo(solver_info);
		std::cout << solver_info << std::endl;
	}

	void topology_optimization(State &state, const std::shared_ptr<CompositeFunctional> j)
	{
		const auto &opt_params = state.args["optimization"];
		const auto &opt_nl_params = state.args["solver"]["optimization_nonlinear"];

		std::shared_ptr<TopologyOptimizationProblem> top_opt = std::make_shared<TopologyOptimizationProblem>(state, j);
		std::shared_ptr<cppoptlib::NonlinearSolver<TopologyOptimizationProblem>> nlsolver = make_nl_solver<TopologyOptimizationProblem>(opt_nl_params);
		nlsolver->setLineSearch(opt_nl_params["line_search"]["method"]);

		Eigen::MatrixXd density_mat = state.assembler.lame_params().density_mat_;
		if (density_mat.size() != state.bases.size())
			density_mat.setZero(state.bases.size(), 1);
		for (const auto &param : opt_params["parameters"])
		{
			if (param["type"] == "topology")
			{
				if (param.contains("initial"))
					density_mat.setConstant(param["initial"]);
				else
				{
					Eigen::MatrixXd barycenters;
					if (state.mesh->is_volume())
					{
						state.mesh->cell_barycenters(barycenters);
						for (int e = 0; e < state.bases.size(); e++)
						{
							density_mat(e) = cross3(barycenters(e,0), barycenters(e,1), barycenters(e,2));
						}
					}
					else
					{
						state.mesh->face_barycenters(barycenters);
						for (int e = 0; e < state.bases.size(); e++)
						{
							density_mat(e) = cross2(barycenters(e,0), barycenters(e,1));
						}
					}
					// density_mat.setOnes();
				}
				
				if (param.contains("power"))
					state.assembler.update_lame_params_density(top_opt->apply_filter(density_mat), param["power"]);
				else
					state.assembler.update_lame_params_density(top_opt->apply_filter(density_mat));
				break;
			}
		}

		Eigen::VectorXd x = density_mat;
		nlsolver->minimize(*top_opt, x);

		json solver_info;
		nlsolver->getInfo(solver_info);
		std::cout << solver_info << std::endl;
	}
} // namespace polyfem