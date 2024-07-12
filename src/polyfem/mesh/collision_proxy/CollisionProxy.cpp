#include "CollisionProxy.hpp"

#include <polyfem/mesh/collision_proxy/UpsampleMesh.hpp>
#include <polyfem/mesh/MeshUtils.hpp>
#include <polyfem/mesh/GeometryReader.hpp>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/MatrixUtils.hpp>

#include <SimpleBVH/BVH.hpp>
#include <igl/edges.h>
#include <igl/barycentric_coordinates.h>
#include <h5pp/h5pp.h>
// #include <fcpw/fcpw.h>

namespace polyfem::mesh
{
	namespace
	{
		template <typename T>
		using RowMajorMatrixX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

		/// @brief Convert from 2D barycentric coordinates (BCs) of a triangle to 3D BCs in a tet.
		/// @param uv 2D BCs of a triangular face
		/// @param local_fid which face do the 2D BCs coorespond to
		/// @return 3D BCs in the tet
		Eigen::MatrixXd uv_to_uvw(const Eigen::MatrixXd &uv, const int local_fid)
		{
			assert(uv.cols() == 2);

			Eigen::MatrixXd uvw = Eigen::MatrixXd::Zero(uv.rows(), 3);
			// u * A + v * B + w * C + (1 - u - v - w) * D
			switch (local_fid)
			{
			case 0:
				uvw.col(0) = uv.col(1);
				uvw.col(1) = uv.col(0);
				break;
			case 1:
				uvw.col(0) = uv.col(0);
				uvw.col(2) = uv.col(1);
				break;
			case 2:
				uvw.leftCols(2) = uv;
				uvw.col(2) = 1 - uv.col(0).array() - uv.col(1).array();
				break;
			case 3:
				uvw.col(1) = uv.col(1);
				uvw.col(2) = uv.col(0);
				break;
			default:
				log_and_throw_error("build_collision_proxy(): unknown local_fid={}", local_fid);
			}
			return uvw;
		}

		Eigen::MatrixXd u_to_uv(const Eigen::MatrixXd &uv, const int local_fid) {
			assert(uv.cols() == 1);  // Ensure uv is one-dimensional

			Eigen::MatrixXd uv2 = Eigen::MatrixXd::Zero(uv.rows(), 2);
			switch (local_fid) {
			case 0:  // Bottom edge
				uv2.col(0) = uv.col(0);  // u varies along x-axis
				uv2.col(1) = Eigen::VectorXd::Zero(uv.rows());  // v is zero along y-axis
				break;
			case 1:  // Right edge
				uv2.col(0) = Eigen::VectorXd::Ones(uv.rows());  // u is one along x-axis
				uv2.col(1) = uv.col(0);  // v varies along y-axis
				break;
			case 2:  // Top edge
				uv2.col(0) = 1 - uv.col(0).array();  // u varies along x-axis (inverted)
				uv2.col(1) = Eigen::VectorXd::Ones(uv.rows());  // v is one along y-axis
				break;
			case 3:  // Left edge
				uv2.col(0) = Eigen::VectorXd::Zero(uv.rows());  // u is zero along x-axis
				uv2.col(1) = 1 - uv.col(0).array();  // v varies along y-axis (inverted)
				break;
			default:
				log_and_throw_error("Invalid local_fid: {}", local_fid);
			}

			return uv2;
		}

		//---------- Added for hex ----------//
		//TODO: Need to define fave order
		Eigen::MatrixXd uv_to_uvw_hex(const Eigen::MatrixXd &uv, const int local_fid)
		{
			assert(uv.cols() == 2);

			Eigen::MatrixXd uvw = Eigen::MatrixXd::Zero(uv.rows(), 3);
			switch (local_fid)
			{
			case 0:
				// Bottom face (Assuming z constant at minimum)
				uvw.col(0) = uv.col(0);
				uvw.col(1) = uv.col(1);
				uvw.col(2) = Eigen::VectorXd::Zero(uv.rows());
				break;
			case 1:
				// Top face (Assuming z constant at maximum)
				uvw.col(0) = uv.col(0);
				uvw.col(1) = uv.col(1);
				uvw.col(2) = Eigen::VectorXd::Ones(uv.rows());
				break;
			case 2:
				// Front face (Assuming y constant at minimum)
				uvw.col(0) = uv.col(0);
				uvw.col(1) = Eigen::VectorXd::Zero(uv.rows());
				uvw.col(2) = uv.col(1);
				break;
			case 3:
				// Back face (Assuming y constant at maximum)
				uvw.col(0) = uv.col(0);
				uvw.col(1) = Eigen::VectorXd::Ones(uv.rows());
				uvw.col(2) = uv.col(1);
				break;
			case 4:
				// Left face (Assuming x constant at minimum)
				uvw.col(0) = Eigen::VectorXd::Zero(uv.rows());
				uvw.col(1) = uv.col(0);
				uvw.col(2) = uv.col(1);
				break;
			case 5:
				// Right face (Assuming x constant at maximum)
				uvw.col(0) = Eigen::VectorXd::Ones(uv.rows());
				uvw.col(1) = uv.col(0);
				uvw.col(2) = uv.col(1);
				break;
			default:
				log_and_throw_error("build_collision_proxy(): unknown local_fid={}", local_fid);
			}
			return uvw;
		}

		Eigen::MatrixXd extract_face_vertices(
			const basis::ElementBases &element, const int local_fid)
		{
			Eigen::MatrixXd UV(3, 2);
			// u * A + b * B + (1-u-v) * C
			UV.row(0) << 1, 0;
			UV.row(1) << 0, 1;
			UV.row(2) << 0, 0;
			const Eigen::MatrixXd UVW = uv_to_uvw(UV, local_fid);

			Eigen::MatrixXd V;
			element.eval_geom_mapping(UVW, V);

			return V;
		}

		// Eigen::MatrixXd extract_edge_vertices(
		// 	const basis::ElementBases &element, const int local_eid)
		// {
		// 	Eigen::MatrixXd UV(2, 1);
		// 	// Linear interpolation between two edge points:
		// 	// At u = 0 (start of the edge)
		// 	// At u = 1 (end of the edge)
		// 	UV.row(0) << 0;  // u = 0
		// 	UV.row(1) << 1;  // u = 1

		// 	const Eigen::MatrixXd UVW = uv_to_uvw_edge(UV, local_eid);

		// 	Eigen::MatrixXd V;
		// 	element.eval_geom_mapping(UVW, V);

		// 	// Ensures that V has two columns, representing the start and end vertices of the edge
		// 	assert(V.cols() == 2 && "The geom_mapping should return two points for an edge");

		// 	return V;
		// }
	} // namespace

	void build_collision_proxy(
		const std::vector<basis::ElementBases> &bases,
		const std::vector<basis::ElementBases> &geom_bases,
		const std::vector<LocalBoundary> &total_local_boundary,
		const int n_bases,
		const int dim,
		const double max_edge_length,
		Eigen::MatrixXd &proxy_vertices,
		Eigen::MatrixXi &proxy_faces,
		std::vector<Eigen::Triplet<double>> &displacement_map_entries,
		const CollisionProxyTessellation tessellation)
	{
		// for each boundary element (f):
		//     tessilate f with triangles of max edge length (fₜ)
		//     for each node (x) of fₜ with global index (i):
		//         for each basis (ϕⱼ) in f's parent element:
		//             set Vᵢ = x where Vᵢ is the i-th proxy vertex
		//             set W(i, j) = ϕⱼ(g⁻¹(x)) where g is the geometry mapping of f
		// caveats:
		// • if x is provided in parametric coordinates, we can skip evaluating g⁻¹
		//   - Vᵢ = g(x) instead
		// • the tessellations of all faces need to be stitched together
		//   - this means duplicate weights should be removed

		std::vector<double> proxy_vertices_list;
		std::vector<int> proxy_faces_list;
		std::vector<Eigen::Triplet<double>> displacement_map_entries_tmp;

		Eigen::MatrixXd UV;
		Eigen::MatrixXi F_local;
		if (tessellation == CollisionProxyTessellation::REGULAR)
		{
			// TODO: use max_edge_length to determine the tessellation
			int n_segments = std::ceil(1.0 / max_edge_length);
			regular_grid_triangle_barycentric_coordinates(/*n=*/n_segments, UV, F_local);
			// regular_grid_triangle_barycentric_coordinates(/*n=*/10, UV, F_local);
		}
		std::cout << "UV:" << UV << std::endl;
		std::cout << "F_local:" << F_local << std::endl;

		for (const LocalBoundary &local_boundary : total_local_boundary)
		{
			if (local_boundary.type() != BoundaryType::TRI)
				log_and_throw_error("build_collision_proxy() is only implemented for tetrahedra!");

			const basis::ElementBases elm = bases[local_boundary.element_id()];
			const basis::ElementBases g = geom_bases[local_boundary.element_id()];
			for (int fi = 0; fi < local_boundary.size(); fi++)
			{
				const int local_fid = local_boundary.local_primitive_id(fi);

				if (tessellation == CollisionProxyTessellation::IRREGULAR)
				{
					// Use the shape of f to determine the tessellation
					const Eigen::MatrixXd node_positions = extract_face_vertices(g, local_fid);
					irregular_triangle_barycentric_coordinates(
						node_positions.row(0), node_positions.row(1), node_positions.row(2),
						max_edge_length, UV, F_local);
				}

				// Convert UV to appropirate UVW based on the local face id
				Eigen::MatrixXd UVW = uv_to_uvw(UV, local_fid);

				Eigen::MatrixXd V_local;
				g.eval_geom_mapping(UVW, V_local);
				assert(V_local.rows() == UV.rows());
				// std::cout << "	fi: " << fi << ", V_local:" << V_local << std::endl;

				const int offset = proxy_vertices_list.size() / dim;
				for (const double x : V_local.reshaped<Eigen::RowMajor>())
					proxy_vertices_list.push_back(x);
				for (const int i : F_local.reshaped<Eigen::RowMajor>())
					proxy_faces_list.push_back(i + offset);

				for (const basis::Basis &basis : elm.bases)
				{
					assert(basis.global().size() == 1);
					const int basis_id = basis.global()[0].index;

					const Eigen::MatrixXd basis_values = basis(UVW);

					for (int i = 0; i < basis_values.size(); i++)
					{
						displacement_map_entries_tmp.emplace_back(
							offset + i, basis_id, basis_values(i));
					}
				}
			}
		}

		// // Print proxy_vertices_list
		// std::cout << "proxy_vertices_list:" << std::endl;
		// for (const auto& value : proxy_vertices_list) {
		// 	std::cout << value << std::endl;
		// }

		// // Print proxy_edges_list
		// std::cout << "proxy_face_list:" << std::endl;
		// for (const auto& edge : proxy_faces_list) {
		// 	std::cout << edge << std::endl;
		// }

		// // Print displacement_map_entries_tmp
		// std::cout << "displacement_map_entries_tmp:" << std::endl;
		// for (const auto& triplet : displacement_map_entries_tmp) {
		// 	std::cout << "(" << triplet.row() << ", " << triplet.col() << ") -> " << triplet.value() << std::endl;
		// }

		// stitch collision proxy together
		stitch_mesh(
			Eigen::Map<RowMajorMatrixX<double>>(proxy_vertices_list.data(), proxy_vertices_list.size() / dim, dim),
			Eigen::Map<RowMajorMatrixX<int>>(proxy_faces_list.data(), proxy_faces_list.size() / dim, dim),
			displacement_map_entries_tmp,
			proxy_vertices, proxy_faces, displacement_map_entries);
		
		// std::cout << "After stitching, proxy_vertices:" << std::endl;
		// std::cout << proxy_vertices << std::endl;
		// std::cout << "After stitching, proxy_faces:" << std::endl;
		// std::cout << proxy_faces << std::endl;
		exit(EXIT_SUCCESS); // Used for debug
	}

	//--------------- Added for 2D iga ---------------//
	//--- Save upsampled mesh function for debug
	void saveAsOBJ(const std::string& filename, const Eigen::MatrixXd& proxy_vertices, const Eigen::MatrixXi& proxy_edges) {
		std::ofstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error opening file for writing: " << filename << std::endl;
			return;
		}

		// Write vertices to file
		for (int i = 0; i < proxy_vertices.rows(); ++i) {
			file << "v " << proxy_vertices(i, 0) << " " << proxy_vertices(i, 1) << " 0" << std::endl;
		}

		// Write edges as lines to file
		for (int i = 0; i < proxy_edges.rows(); ++i) {
			// .obj files are 1-indexed, so add 1 to each index
			file << "l " << (proxy_edges(i, 0) + 1) << " " << (proxy_edges(i, 1) + 1) << std::endl;
		}

		file.close();
	}

	bool areDisplacementMapsEqual(const std::vector<Eigen::Triplet<double>>& map1, const std::vector<Eigen::Triplet<double>>& map2) {
		if (map1.size() != map2.size()) {
			std::cout << "Different number of entries." << std::endl;
			return false;
		}

		for (size_t i = 0; i < map1.size(); ++i) {
			const Eigen::Triplet<double>& triplet1 = map1[i];
			const Eigen::Triplet<double>& triplet2 = map2[i];

			if (triplet1.row() != triplet2.row() ||
				triplet1.col() != triplet2.col() ||
				triplet1.value() != triplet2.value()) {
				std::cout << "Mismatch found at index " << i << ": "
						<< "(" << triplet1.row() << ", " << triplet1.col() << ", " << triplet1.value() << ") != "
						<< "(" << triplet2.row() << ", " << triplet2.col() << ", " << triplet2.value() << ")" << std::endl;
				return false;
			}
		}

		return true;
	}

	void build_collision_proxy_quad(
		const std::vector<basis::ElementBases> &bases,
		const std::vector<basis::ElementBases> &geom_bases,
		const std::vector<LocalBoundary> &total_local_boundary,
		const int n_bases,
		const int dim,
		const double max_edge_length,
		Eigen::MatrixXd &proxy_vertices,
		Eigen::MatrixXi &proxy_edges,
		std::vector<Eigen::Triplet<double>> &displacement_map_entries,
		const CollisionProxyTessellation tessellation) 
	{
		std::vector<double> proxy_vertices_list;
		std::vector<int> proxy_edges_list;
		std::vector<Eigen::Triplet<double>> displacement_map_entries_tmp;

		//--- Define per-boundary element discretization
		Eigen::MatrixXd U;
		Eigen::MatrixXi E_local;
		if (tessellation == CollisionProxyTessellation::REGULAR)
		{
			// Define UV as linear interpolation between 0 and 1 with steps based on max_edge_length
			int n_segments = std::ceil(1.0 / max_edge_length);
			U.resize(n_segments + 1, 1);  // This defines the number of vertices
			for (int i = 0; i <= n_segments; ++i) {
				U(i, 0) = double(i) / n_segments;
			}

			E_local.resize(n_segments, 2); // Correctly define the dimensions for edges
			for (int i = 0; i < n_segments; ++i) {
				E_local(i, 0) = i;     // Start vertex of edge i
				E_local(i, 1) = i + 1; // End vertex of edge i
			}
		}
		// std::cout << "U:" << U << std::endl;
		// std::cout << "E_local has " << E_local.rows() << " rows and " << E_local.cols() << " columns." << std::endl;
		// std::cout << "E_local:" << E_local << std::endl;

		//--- Itetrate boundary elements
		for (const LocalBoundary &local_boundary : total_local_boundary)
		{
			if (local_boundary.type() != BoundaryType::QUAD_LINE)
				log_and_throw_error("build_collision_proxy() is only implemented for lines in a quad mesh context!");

			const basis::ElementBases elm = bases[local_boundary.element_id()];
			const basis::ElementBases g = geom_bases[local_boundary.element_id()];
			for (int ei = 0; ei < local_boundary.size(); ei++)
			{
				const int local_eid = local_boundary.local_primitive_id(ei);
				// std::cout << "	ei: " << ei << ", local_eid: " << local_eid << std::endl;

				Eigen::MatrixXd UV = u_to_uv(U, local_eid);

				Eigen::MatrixXd V_local;
				g.eval_geom_mapping(UV, V_local);
				assert(V_local.rows() == U.rows());
				// std::cout << "	ei: " << ei << ", V_local:" << V_local << std::endl;

				const int offset = proxy_vertices_list.size() / dim;
				for (const double x : V_local.reshaped<Eigen::RowMajor>())
					proxy_vertices_list.push_back(x);
				for (const int i : E_local.reshaped<Eigen::RowMajor>())
					proxy_edges_list.push_back(i + offset);

				for (const basis::Basis &basis : elm.bases)
				{
					assert(basis.global().size() == 1);
					const int basis_id = basis.global()[0].index;

					const Eigen::MatrixXd basis_values = basis(UV);

					for (int i = 0; i < basis_values.size(); i++)
					{
						displacement_map_entries_tmp.emplace_back(
							offset + i, basis_id, basis_values(i));
					}
				}
			}
			// std::cout << std::endl;
		}

		// // Print proxy_vertices_list
		// std::cout << "proxy_vertices_list:" << std::endl;
		// for (const auto& value : proxy_vertices_list) {
		// 	std::cout << value << std::endl;
		// }
		// // Print proxy_edges_list
		// std::cout << "proxy_edges_list:" << std::endl;
		// for (const auto& edge : proxy_edges_list) {
		// 	std::cout << edge << std::endl;
		// }
		// // Print displacement_map_entries_tmp
		// std::cout << "displacement_map_entries_tmp:" << std::endl;
		// for (const auto& triplet : displacement_map_entries_tmp) {
		// 	std::cout << "(" << triplet.row() << ", " << triplet.col() << ") -> " << triplet.value() << std::endl;
		// }
		
		//--- stitch collision proxy together
		stitch_line_mesh(
			Eigen::Map<RowMajorMatrixX<double>>(proxy_vertices_list.data(), proxy_vertices_list.size() / dim, dim),
			Eigen::Map<RowMajorMatrixX<int>>(proxy_edges_list.data(), proxy_edges_list.size() / 2, 2),
			displacement_map_entries_tmp,
			proxy_vertices, proxy_edges, displacement_map_entries);

		// std::cout << proxy_vertices << std::endl;
		// std::cout << proxy_edges << std::endl;
		saveAsOBJ("/Users/liujunyu/Desktop/Research/UVic_NYU/IGA_IPC/code/polyfem/experiments/output/collision_mesh.obj", proxy_vertices, proxy_edges);
		// bool isEqual = areDisplacementMapsEqual(displacement_map_entries_tmp, displacement_map_entries);
		// std::cout << "The displacement maps are " << (isEqual ? "equal." : "not equal.") << std::endl;
		// std::cout << "displacement_map_entries:" << std::endl;
		// for (const auto& triplet : displacement_map_entries) {
		// 	std::cout << "(" << triplet.row() << ", " << triplet.col() << ") -> " << triplet.value() << std::endl;
		// }
		exit(EXIT_SUCCESS); // Used for debug
	}

	//--------------- Added for 3D iga ---------------//
	void build_collision_proxy_hex(
		const std::vector<basis::ElementBases> &bases,
		const std::vector<basis::ElementBases> &geom_bases,
		const std::vector<LocalBoundary> &total_local_boundary,
		const int n_bases,
		const int dim,
		const double max_edge_length,
		Eigen::MatrixXd &proxy_vertices,
		Eigen::MatrixXi &proxy_faces,
		std::vector<Eigen::Triplet<double>> &displacement_map_entries,
		const CollisionProxyTessellation tessellation)
	{
		// std::cout << "Check in build_collision_proxy_hex" << std::endl;
		
		std::vector<double> proxy_vertices_list;
		std::vector<int> proxy_faces_list;
		std::vector<Eigen::Triplet<double>> displacement_map_entries_tmp;

		Eigen::MatrixXd UV;
		Eigen::MatrixXi F_local;
		if (tessellation == CollisionProxyTessellation::REGULAR)
		{
			// Adjust tessellation to use quadrilaterals instead of triangles
			regular_grid_quadrilateral_barycentric_coordinates(/*n=*/10, UV, F_local);
		}

		for (const LocalBoundary &local_boundary : total_local_boundary)
		{
			if (local_boundary.type() != BoundaryType::QUAD)
				log_and_throw_error("build_collision_proxy() is currently implemented for hexahedral elements!");

			const basis::ElementBases elm = bases[local_boundary.element_id()];
			const basis::ElementBases g = geom_bases[local_boundary.element_id()];
			for (int fi = 0; fi < local_boundary.size(); fi++)
			{
				const int local_fid = local_boundary.local_primitive_id(fi);
				//TODO: Check how does the local_primitive_id for hex elements are defined

				// if (tessellation == CollisionProxyTessellation::IRREGULAR)
				// {
				// 	const Eigen::MatrixXd node_positions = extract_face_vertices(g, local_fid);
				// 	irregular_quadrilateral_barycentric_coordinates(
				// 		node_positions.row(0), node_positions.row(1), node_positions.row(2), node_positions.row(3),
				// 		max_edge_length, UV, F_local);
				// }

				Eigen::MatrixXd UVW = uv_to_uvw_hex(UV, local_fid);

				Eigen::MatrixXd V_local;
				g.eval_geom_mapping(UVW, V_local);
				assert(V_local.rows() == UV.rows());

				const int offset = proxy_vertices_list.size() / dim;
				for (const double x : V_local.reshaped<Eigen::RowMajor>())
					proxy_vertices_list.push_back(x);
				for (const int i : F_local.reshaped<Eigen::RowMajor>())
					proxy_faces_list.push_back(i + offset);

				for (const basis::Basis &basis : elm.bases)
				{
					assert(basis.global().size() == 1);
					const int basis_id = basis.global()[0].index;

					// const Eigen::MatrixXd basis_values = basis.eval(UVW);
					Eigen::MatrixXd basis_values;
					basis.eval_basis(UVW, basis_values);

					for (int i = 0; i < basis_values.size(); i++)
					{
						displacement_map_entries_tmp.emplace_back(
							offset + i, basis_id, basis_values(i));
					}
				}
			}
		}
		
		// Final stitching of the mesh
		stitch_mesh(
			Eigen::Map<Eigen::MatrixXd>(proxy_vertices_list.data(), proxy_vertices_list.size() / dim, dim),
			Eigen::Map<Eigen::MatrixXi>(proxy_faces_list.data(), proxy_faces_list.size() / 3, 3),
			displacement_map_entries_tmp,
			proxy_vertices, proxy_faces, displacement_map_entries);
	}

	// ========================================================================

	void build_collision_proxy_displacement_map(
		const std::vector<basis::ElementBases> &bases,
		const std::vector<basis::ElementBases> &geom_bases,
		const std::vector<mesh::LocalBoundary> &total_local_boundary,
		const int n_bases,
		const int dim,
		const Eigen::MatrixXd &proxy_vertices,
		std::vector<Eigen::Triplet<double>> &displacement_map_entries)
	{
		// for each vᵢ in proxy_vertices:
		//     find closest element (t)
		//     compute v̂ᵢ = g⁻¹(v) where g is the geometry mapping of t
		//     for each basis (ϕⱼ) in t:
		//         set W(i, j) = ϕⱼ(v̂ᵢ)

		// NOTE: this is only implemented for P1 geometry
		for (const basis::ElementBases &element_bases : geom_bases)
		{
			for (const basis::Basis &basis : element_bases.bases)
			{
				if (basis.order() != 1)
					log_and_throw_error("build_collision_proxy_displacement_map() is only implemented for P1 geometry!");
			}
		}

		// --------------------------------------------------------------------
		// build a BVH to find which element each proxy vertex belongs to
		std::vector<std::array<Eigen::Vector3d, 2>> boxes(
			geom_bases.size(), {{Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()}});
		for (int i = 0; i < geom_bases.size(); i++)
		{
			const Eigen::MatrixXd nodes = geom_bases[i].nodes();
			boxes[i][0].head(dim) = nodes.colwise().minCoeff();
			boxes[i][1].head(dim) = nodes.colwise().maxCoeff();
		}

		SimpleBVH::BVH bvh;
		bvh.init(boxes);

		// --------------------------------------------------------------------
		// build wide BVH to find closest point on the surface

		// assert(dim == 3);
		// fcpw::Scene<3> scene;
		// {
		// 	std::vector<Eigen::Vector3d> vertices;
		// 	vertices.reserve(geom_bases.size() * 4);
		// 	std::vector<std::array<int, 3>> faces;
		// 	faces.reserve(geom_bases.size() * 4);
		// 	for (const basis::ElementBases &elm : geom_bases)
		// 	{
		// 		const int offset = vertices.size();
		// 		for (int i = 0; i < 4; i++)
		// 			vertices.push_back(elm.bases[i].global()[0].node);
		// 		faces.push_back({{offset + 0, offset + 1, offset + 2}});
		// 		faces.push_back({{offset + 0, offset + 1, offset + 3}});
		// 		faces.push_back({{offset + 0, offset + 2, offset + 3}});
		// 		faces.push_back({{offset + 1, offset + 2, offset + 3}});
		// 	}

		// 	scene.setObjectTypes({{fcpw::PrimitiveType::Triangle}});
		// 	scene.setObjectVertexCount(vertices.size(), /*object_id=*/0);
		// 	scene.setObjectTriangleCount(faces.size(), /*object_id=*/0);

		// 	for (int i = 0; i < vertices.size(); i++)
		// 	{
		// 		scene.setObjectVertex(vertices[i].cast<float>(), i, /*object_id=*/0);
		// 	}

		// 	// specify the triangle indices
		// 	for (int i = 0; i < faces.size(); i++)
		// 	{
		// 		scene.setObjectTriangle(faces[i].data(), i, 0);
		// 	}

		// 	scene.build(fcpw::AggregateType::Bvh_SurfaceArea, true); // the second boolean argument enables vectorization
		// }

		// --------------------------------------------------------------------

		for (int i = 0; i < proxy_vertices.rows(); i++)
		{
			Eigen::Vector3d v = Eigen::Vector3d::Zero();
			v.head(dim) = proxy_vertices.row(i);

			std::vector<unsigned int> candidates;
			bvh.intersect_box(v, v, candidates);

			// find which element the proxy vertex belongs to
			int closest_element_id = -1;
			VectorNd vhat;
			for (const unsigned int element_id : candidates)
			{
				const Eigen::MatrixXd nodes = geom_bases[element_id].nodes();
				if (dim == 2)
				{
					assert(nodes.rows() == 3);
					Eigen::RowVector3d bc;
					igl::barycentric_coordinates(
						proxy_vertices.row(i), nodes.row(0), nodes.row(1), nodes.row(2), bc);
					vhat = bc.head<2>();
				}
				else
				{
					assert(dim == 3 && nodes.rows() == 4);
					Eigen::RowVector4d bc;
					igl::barycentric_coordinates(
						proxy_vertices.row(i), nodes.row(0), nodes.row(1), nodes.row(2), nodes.row(3), bc);
					vhat = bc.head<3>();
				}
				if (vhat.minCoeff() >= 0 && vhat.maxCoeff() <= 1 && vhat.sum() <= 1)
				{
					closest_element_id = element_id;
					break;
				}
			}

			if (closest_element_id < 0)
			{
				// perform a closest point query
				log_and_throw_error("build_collision_proxy_displacement_map(): closest point query not implemented!");
				// fcpw::Interaction<3> cpq_interaction;
				// scene.findClosestPoint(v.cast<float>(), cpq_interaction);

				// closest_element_id = cpq_interaction.primitiveIndex / 4;

				// const Eigen::MatrixXd nodes = geom_bases[closest_element_id].nodes();

				// assert(dim == 3 && nodes.rows() == 4);
				// Eigen::RowVector4d bc;
				// igl::barycentric_coordinates(
				// 	proxy_vertices.row(i), nodes.row(0), nodes.row(1), nodes.row(2), nodes.row(3), bc);
				// vhat = bc.head<3>();
			}

			// compute the displacement map entries
			for (const basis::Basis &basis : bases[closest_element_id].bases)
			{
				assert(basis.global().size() == 1);
				const int j = basis.global()[0].index;
				displacement_map_entries.emplace_back(i, j, basis(vhat.transpose())(0));
			}
		}
	}

	// ========================================================================

	void load_collision_proxy(
		const std::string &mesh_filename,
		const std::string &weights_filename,
		const Eigen::VectorXi &in_node_to_node,
		const json &transformation,
		Eigen::MatrixXd &vertices,
		Eigen::VectorXi &codim_vertices,
		Eigen::MatrixXi &edges,
		Eigen::MatrixXi &faces,
		std::vector<Eigen::Triplet<double>> &displacement_map_entries)
	{
		load_collision_proxy_mesh(mesh_filename, transformation, vertices, codim_vertices, edges, faces);
		std::cout << "Finish loading collision proxy mesh" << std::endl;
		load_collision_proxy_displacement_map(weights_filename, in_node_to_node, vertices.rows(), displacement_map_entries);
		std::cout << "Finish loading collision proxy displacement map" << std::endl;
	}

	void load_collision_proxy_mesh(
		const std::string &mesh_filename,
		const json &transformation,
		Eigen::MatrixXd &vertices,
		Eigen::VectorXi &codim_vertices,
		Eigen::MatrixXi &edges,
		Eigen::MatrixXi &faces)
	{
		Eigen::MatrixXi codim_edges;
		read_surface_mesh(mesh_filename, vertices, codim_vertices, codim_edges, faces);

		if (faces.size())
			igl::edges(faces, edges);

		utils::append_rows(edges, codim_edges);

		// Transform the collision proxy
		std::array<RowVectorNd, 2> bbox;
		bbox[0] = vertices.colwise().minCoeff();
		bbox[1] = vertices.colwise().maxCoeff();

		MatrixNd A;
		VectorNd b;
		// TODO: pass correct unit scale
		construct_affine_transformation(
			/*unit_scale=*/1, transformation,
			(bbox[1] - bbox[0]).cwiseAbs().transpose(),
			A, b);
		vertices = (vertices * A.transpose()).rowwise() + b.transpose();
	}

	void load_collision_proxy_displacement_map(
		const std::string &weights_filename,
		const Eigen::VectorXi &in_node_to_node,
		const size_t num_proxy_vertices,
		std::vector<Eigen::Triplet<double>> &displacement_map_entries)
	{
#ifndef NDEBUG
		const size_t num_fe_nodes = in_node_to_node.size();
#endif
		h5pp::File file(weights_filename, h5pp::FileAccess::READONLY);
		const std::array<long, 2> shape = file.readAttribute<std::array<long, 2>>("weight_triplets", "shape");
		assert(shape[0] == num_proxy_vertices && shape[1] == num_fe_nodes);
		Eigen::VectorXd values = file.readDataset<Eigen::VectorXd>("weight_triplets/values");
		Eigen::VectorXi rows = file.readDataset<Eigen::VectorXi>("weight_triplets/rows");
		Eigen::VectorXi cols = file.readDataset<Eigen::VectorXi>("weight_triplets/cols");
		std::cout << "Load displacement map: Finish loading values, rows, cols." << std::endl;
		assert(rows.maxCoeff() < num_proxy_vertices);
		assert(cols.maxCoeff() < num_fe_nodes);

		// TODO: use these to build the in_node_to_node map
		// const Eigen::VectorXi in_ordered_vertices = file.exist("ordered_vertices") ? H5Easy::load<Eigen::VectorXi>(file, "ordered_vertices") : mesh->in_ordered_vertices();
		// const Eigen::MatrixXi in_ordered_edges = file.exist("ordered_edges") ? H5Easy::load<Eigen::MatrixXi>(file, "ordered_edges") : mesh->in_ordered_edges();
		// const Eigen::MatrixXi in_ordered_faces = file.exist("ordered_faces") ? H5Easy::load<Eigen::MatrixXi>(file, "ordered_faces") : mesh->in_ordered_faces();
		// std::cout << "Load displacement map: Finishn loading vertices, edges, and faces." << std::endl;

		displacement_map_entries.clear();
		displacement_map_entries.reserve(values.size());

		assert(in_node_to_node.size() == num_fe_nodes);
		std::cout << "Load displacement map: All assertion passed." << std::endl;
		std::cout << "Value size: " << values.size() << std::endl;
		std::cout << "Row size: " << rows.size() << std::endl;
		std::cout << "Col size: " << cols.size() << std::endl;
		std::cout << "in_node_to_node size: " << in_node_to_node.size() << std::endl;
		for (int i = 0; i < values.size(); i++)
		{
			// Rearrange the columns based on the FEM mesh node order
			displacement_map_entries.emplace_back(rows[i], in_node_to_node[cols[i]], values[i]);
		}
	}
} // namespace polyfem::mesh