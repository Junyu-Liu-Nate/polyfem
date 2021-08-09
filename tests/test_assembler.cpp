#include <polyfem/State.hpp>

#include <catch.hpp>
#include <iostream>

using namespace polyfem;

TEST_CASE("hessian_lin", "[assembler]")
{
	const std::string path = POLYFEM_DATA_DIR;
	json in_args = json({});
	in_args["mesh"] = path + "/circle2.msh";
	in_args["force_linear_geometry"] = true;
	in_args["problem"] = "ElasticExact";

	State state;
	state.init_logger("", 6, false);
	state.init(in_args);
	state.load_mesh();

	state.compute_mesh_stats();
	state.build_basis();

	state.assemble_stiffness_mat();

	SpareMatrixCache mat_cache;
	StiffnessMatrix hessian;
	Eigen::MatrixXd disp(state.n_bases * 2, 1);
	disp.setZero();

	for (int rand = 0; rand < 10; ++rand)
	{
		state.assembler.assemble_energy_hessian(
			"LinearElasticity", false, state.n_bases, false,
			state.bases, state.bases, state.ass_vals_cache, disp, mat_cache, hessian);

		const StiffnessMatrix tmp = state.stiffness - hessian;
		const auto val = Approx(0).margin(1e-8);

		for (int k = 0; k < tmp.outerSize(); ++k)
		{
			for (StiffnessMatrix::InnerIterator it(tmp, k); it; ++it)
			{
				REQUIRE(it.value() == val);
			}
		}

		disp.setRandom();
	}
}