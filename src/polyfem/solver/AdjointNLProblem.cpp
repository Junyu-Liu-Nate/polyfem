#include "AdjointNLProblem.hpp"

#include <polyfem/solver/forms/adjoint_forms/AdjointForm.hpp>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/MaybeParallelFor.hpp>
#include <polyfem/utils/Timer.hpp>
#include <polyfem/io/OBJWriter.hpp>
#include <polyfem/State.hpp>
#include <igl/boundary_facets.h>
#include <igl/writeOBJ.h>

#include <polyfem/solver/NLProblem.hpp>

#include <list>
#include <stack>

namespace polyfem::solver
{
	namespace
	{
		using namespace std;
		// Class to represent a graph
		class Graph
		{
			int V; // No. of vertices'

			// adjacency lists
			vector<list<int>> adj;

			// A function used by topologicalSort
			void topologicalSortUtil(int v, vector<bool> &visited, stack<int> &Stack);

		public:
			Graph(int V); // Constructor

			// function to add an edge to graph
			void addEdge(int v, int w);

			// prints a Topological Sort of the complete graph
			vector<int> topologicalSort();
		};

		Graph::Graph(int V)
		{
			this->V = V;
			adj.resize(V);
		}

		void Graph::addEdge(int v, int w)
		{
			adj[v].push_back(w); // Add w to v’s list.
		}

		// A recursive function used by topologicalSort
		void Graph::topologicalSortUtil(int v, vector<bool> &visited,
										stack<int> &Stack)
		{
			// Mark the current node as visited.
			visited[v] = true;

			// Recur for all the vertices adjacent to this vertex
			list<int>::iterator i;
			for (i = adj[v].begin(); i != adj[v].end(); ++i)
				if (!visited[*i])
					topologicalSortUtil(*i, visited, Stack);

			// Push current vertex to stack which stores result
			Stack.push(v);
		}

		// The function to do Topological Sort. It uses recursive
		// topologicalSortUtil()
		vector<int> Graph::topologicalSort()
		{
			stack<int> Stack;

			// Mark all the vertices as not visited
			vector<bool> visited(V, false);

			// Call the recursive helper function to store Topological
			// Sort starting from all vertices one by one
			for (int i = 0; i < V; i++)
				if (visited[i] == false)
					topologicalSortUtil(i, visited, Stack);

			// Print contents of stack
			vector<int> sorted;
			while (Stack.empty() == false)
			{
				sorted.push_back(Stack.top());
				Stack.pop();
			}

			return sorted;
		}
	} // namespace

	AdjointNLProblem::AdjointNLProblem(std::shared_ptr<AdjointForm> form, const std::vector<std::shared_ptr<VariableToSimulation>> &variables_to_simulation, const std::vector<std::shared_ptr<State>> &all_states, const json &args)
		: FullNLProblem({form}),
		  form_(form),
		  variables_to_simulation_(variables_to_simulation),
		  all_states_(all_states),
		  save_freq(args["output"]["save_frequency"]),
		  solve_in_parallel(args["solver"]["advanced"]["solve_in_parallel"])
	{
		cur_grad.setZero(0);

		solve_in_order.clear();
		{
			Graph G(all_states.size());
			for (int k = 0; k < all_states.size(); k++)
			{
				auto &arg = args["states"][k];
				if (arg["initial_guess"].get<int>() >= 0)
					G.addEdge(arg["initial_guess"].get<int>(), k);
			}

			solve_in_order = G.topologicalSort();
		}

		active_state_mask.assign(all_states_.size(), false);
		for (int i = 0; i < all_states_.size(); i++)
		{
			for (const auto &v2sim : variables_to_simulation_)
			{
				for (const auto &state : v2sim->get_states())
				{
					if (all_states_[i].get() == state.get())
					{
						active_state_mask[i] = true;
						break;
					}
				}
			}
		}
	}

	AdjointNLProblem::AdjointNLProblem(std::shared_ptr<AdjointForm> form, const std::vector<std::shared_ptr<AdjointForm>> stopping_conditions, const std::vector<std::shared_ptr<VariableToSimulation>> &variables_to_simulation, const std::vector<std::shared_ptr<State>> &all_states, const json &args) : AdjointNLProblem(form, variables_to_simulation, all_states, args)
	{
		stopping_conditions_ = stopping_conditions;
	}

	void AdjointNLProblem::hessian(const Eigen::VectorXd &x, StiffnessMatrix &hessian)
	{
		log_and_throw_adjoint_error("Hessian not supported!");
	}

	double AdjointNLProblem::value(const Eigen::VectorXd &x)
	{
		return form_->value(x);
	}

	void AdjointNLProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &gradv)
	{
		if (cur_grad.size() == x.size())
			gradv = cur_grad;
		else
		{
			gradv.setZero(x.size());

			{
				POLYFEM_SCOPED_TIMER("adjoint solve");
				for (int i = 0; i < all_states_.size(); i++)
					all_states_[i]->solve_adjoint_cached(form_->compute_adjoint_rhs(x, *all_states_[i])); // caches inside state
			}

			{
				POLYFEM_SCOPED_TIMER("gradient assembly");
				form_->first_derivative(x, gradv);
			}

			cur_grad = gradv;
		}
	}

	bool AdjointNLProblem::is_step_valid(const Eigen::VectorXd &x0, const Eigen::VectorXd &x1) const
	{
		return form_->is_step_valid(x0, x1);
	}

	bool AdjointNLProblem::is_step_collision_free(const Eigen::VectorXd &x0, const Eigen::VectorXd &x1) const
	{
		return form_->is_step_collision_free(x0, x1);
	}

	double AdjointNLProblem::max_step_size(const Eigen::VectorXd &x0, const Eigen::VectorXd &x1) const
	{
		return form_->max_step_size(x0, x1);
	}

	void AdjointNLProblem::line_search_begin(const Eigen::VectorXd &x0, const Eigen::VectorXd &x1)
	{
		form_->line_search_begin(x0, x1);
	}

	void AdjointNLProblem::line_search_end()
	{
		form_->line_search_end();
	}

	void AdjointNLProblem::post_step(const polysolve::nonlinear::PostStepData &data)
	{
		save_to_file(data.iter_num, data.x);

		form_->post_step(data);
	}

	void AdjointNLProblem::save_to_file(const int iter_num, const Eigen::VectorXd &x0)
	{
		int id = 0;
		if (iter_num % save_freq != 0)
			return;
		adjoint_logger().info("Saving iteration {}", iter_num);
		for (const auto &state : all_states_)
		{
			bool save_vtu = true;
			bool save_rest_mesh = true;

			std::string vis_mesh_path = state->resolve_output_path(fmt::format("opt_state_{:d}_iter_{:d}.vtu", id, iter_num));
			std::string rest_mesh_path = state->resolve_output_path(fmt::format("opt_state_{:d}_iter_{:d}.obj", id, iter_num));
			id++;

			if (!save_vtu)
				continue;
			adjoint_logger().debug("Save final vtu to file {} ...", vis_mesh_path);

			double tend = state->args.value("tend", 1.0);
			double dt = 1;
			if (!state->args["time"].is_null())
				dt = state->args["time"]["dt"];

			Eigen::MatrixXd sol = state->diff_cached.u(-1);

			state->out_geom.save_vtu(
				vis_mesh_path,
				*state,
				sol,
				Eigen::MatrixXd::Zero(state->n_pressure_bases, 1),
				tend, dt,
				io::OutGeometryData::ExportOptions(state->args, state->mesh->is_linear(), state->problem->is_scalar(), state->solve_export_to_file),
				state->is_contact_enabled(),
				state->solution_frames);

			if (!save_rest_mesh)
				continue;
			adjoint_logger().debug("Save rest mesh to file {} ...", rest_mesh_path);

			// If shape opt, save rest meshes as well
			Eigen::MatrixXd V;
			Eigen::MatrixXi F;
			state->get_vertices(V);
			state->get_elements(F);
			if (state->mesh->dimension() == 3)
				F = igl::boundary_facets<Eigen::MatrixXi, Eigen::MatrixXi>(F);

			io::OBJWriter::write(rest_mesh_path, V, F);
		}
	}

	void AdjointNLProblem::solution_changed_no_solve(const Eigen::VectorXd &newX)
	{
		bool need_rebuild_basis = false;

		// update to new parameter and check if the new parameter is valid to solve
		for (const auto &v : variables_to_simulation_)
		{
			v->update(newX);
			if (v->get_parameter_type() == ParameterType::Shape)
				need_rebuild_basis = true;
		}

		if (need_rebuild_basis)
		{
			for (const auto &state : all_states_)
				state->build_basis();
		}

		form_->solution_changed(newX);
	}

	void AdjointNLProblem::solution_changed(const Eigen::VectorXd &newX)
	{
		bool need_rebuild_basis = false;

		// update to new parameter and check if the new parameter is valid to solve
		for (const auto &v : variables_to_simulation_)
		{
			v->update(newX);
			if (v->get_parameter_type() == ParameterType::Shape)
				need_rebuild_basis = true;
		}

		if (need_rebuild_basis)
		{
			for (const auto &state : all_states_)
				state->build_basis();
		}

		// solve PDE
		solve_pde();

		form_->solution_changed(newX);
	}

	void AdjointNLProblem::solve_pde()
	{
		if (solve_in_parallel)
		{
			adjoint_logger().info("Run simulations in parallel...");

			utils::maybe_parallel_for(all_states_.size(), [&](int start, int end, int thread_id) {
				for (int i = start; i < end; i++)
				{
					auto state = all_states_[i];
					if (active_state_mask[i] || state->diff_cached.size() == 0)
					{
						state->assemble_rhs();
						state->assemble_mass_mat();
						Eigen::MatrixXd sol, pressure; // solution is also cached in state
						state->solve_problem(sol, pressure);
					}
				}
			});
		}
		else
		{
			Eigen::MatrixXd sol, pressure; // solution is also cached in state
			for (int i : solve_in_order)
			{
				auto state = all_states_[i];
				if (active_state_mask[i] || state->diff_cached.size() == 0)
				{
					state->assemble_rhs();
					state->assemble_mass_mat();

					state->solve_problem(sol, pressure);
				}
			}
		}

		cur_grad.resize(0);
	}

	bool AdjointNLProblem::stop(const TVector &x)
	{
		if (stopping_conditions_.size() == 0)
			return false;

		for (auto &obj : stopping_conditions_)
		{
			obj->solution_changed(x);
			if (obj->value(x) > 0)
				return false;
		}
		return true;
	}

} // namespace polyfem::solver