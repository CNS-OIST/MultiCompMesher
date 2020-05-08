/*
Multi-Component Mesh Generator
Copyright (C) 2020 Okinawa Institute of Science and Technology, Japan.

Developer: Weiliang Chen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_to_labeling_function_wrapper.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_constant_domain_field_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Timer.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/make_mesh_3.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "mesh_repair.h"
#include "utility.h"

namespace po = boost::program_options;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::FT FT;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Submesh_domain;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;
typedef CGAL::Mesh_constant_domain_field_3<Mesh_domain::R, Mesh_domain::Index> Sizing_field;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

template <typename FT, typename P>
class FT_to_point_function_wrapper: public CGAL::cpp98::unary_function<P, FT> {
    std::function<FT(P)> function;

  public:
    typedef P Point;
    FT_to_point_function_wrapper(std::function<FT(P)> f)
        : function(f) {}
    FT operator()(P p) const {
        return function(p);
    }
};
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function> Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;

FT func(Point p, Submesh_domain* d) {
    if (d->is_in_domain_object()(p)) {
        return -1.0;
    } else {
        return 1.0;
    }
}

int main(int argc, char* argv[]) {
#ifdef CGAL_CONCURRENT_MESH_3
    std::cout << "Running in concurrency mode.\n";
#else
    std::cout << "Running in sequential mode.\n";
#endif

    try {
        po::options_description desc("Options");
        desc.add_options()("help,h", "Help message")("boundary-file,ib",
                                                     po::value<std::string>(),
                                                     "Boundary file")("component-file,ic",
                                                                      po::value<std::string>(),
                                                                      "Component file")(
            "output,o", po::value<std::string>(), "Output mesh file")(
            "size-field-file", po::value<std::string>(), "Special size field setting file.")(
            "fc-angle", po::value<double>()->default_value(25.0), "Facet criteria - Angle")(
            "fc-size", po::value<double>()->default_value(25.0), "Facet criteria - Size")(
            "fc-distance", po::value<double>()->default_value(5.0), "Facet criteria - Distance")(
            "cc-ratio",
            po::value<double>()->default_value(2.0),
            "Cell criteria - Cell radius edge ratio")("cc-size",
                                                      po::value<double>()->default_value(25.0),
                                                      "Cell criteria - Size")(
            "odt",
            po::bool_switch()->default_value(false),
            "Enable ODT mesh optimization")("odt-time",
                                            po::value<double>()->default_value(10.0),
                                            "Time limit for ODT mesh optimization, in second")(
            "lloyd",
            po::bool_switch()->default_value(false),
            "Enable Lloyd mesh optimization")("lloyd-time",
                                              po::value<double>()->default_value(10.0),
                                              "Time limit for Lloyd mesh optimization, in second")(
            "perturb",
            po::bool_switch()->default_value(false),
            "Enable mesh sliver perturbation")("perturb-time",
                                               po::value<double>()->default_value(10.0),
                                               "Time limit for sliver perturbation, in second")(
            "perturb-bound",
            po::value<double>()->default_value(0.0),
            "Targeted lower bound on dihedral angles of mesh cells for sliver perturbation, in "
            "degree")("exude",
                      po::bool_switch()->default_value(false),
                      "Enable mesh sliver exudation")("exude-time",
                                                      po::value<double>()->default_value(10.0),
                                                      "Time limit for sliver exudation, in second")(
            "exude-bound",
            po::value<double>()->default_value(0.0),
            "Targeted lower bound on dihedral angles of mesh cells for sliver exudation, in "
            "degree")("manifold",
                      po::value<uint>()->default_value(0),
                      "Mainfold restriction of the outout mesh. (0) No restriction (1) Manifold "
                      "(2) Manifold with boundaries");

        po::positional_options_description p;
        p.add("boundary-file", 1);
        p.add("component-file", 1);
        p.add("output", 1);

        po::variables_map vm;

        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << "Usage: ./MultiCompMesher boundary-file component-file output [options]\n";
            std::cout << desc;
            return 0;
        }

        if (vm.count("boundary-file")) {
            std::cout << "Boundary file is: " << vm["boundary-file"].as<std::string>() << "\n";
        } else {
            std::cerr << "Boundary file is required.\n";
            std::cerr << "Usage: ./MultiCompMesher boundary-file component-file output [options]\n";
            std::cerr << desc;
            return EXIT_FAILURE;
        }
        if (vm.count("component-file")) {
            std::cout << "Component file is: " << vm["component-file"].as<std::string>() << "\n";
        } else {
            std::cerr << "Component file is required.\n";
            std::cerr << "Usage: ./MultiCompMesher boundary-file component-file output [options]\n";
            std::cerr << desc;
            return EXIT_FAILURE;
        }

        std::string output_file;
        if (vm.count("output")) {
            if (ends_with(vm["output"].as<std::string>(), ".mesh")) {
                output_file = vm["output"].as<std::string>();
            } else {
                output_file = vm["output"].as<std::string>() + std::string(".mesh");
            }
            std::cout << "Mesh will write to: " << output_file << "\n\n";
        } else {
            std::cerr << "Output file name is required.\n";
            std::cerr << "Usage: ./MultiCompMesher boundary-file component-file output [options]\n";
            std::cerr << desc;
            return EXIT_FAILURE;
        }

        std::vector<std::string> boundary_files;
        std::vector<std::string> domain_signs;

        bool result = getFileContent(vm["boundary-file"].as<std::string>(), boundary_files);

        if (result) {
            std::cout << "boundary mesh files:\n";
            for (std::string& line: boundary_files) {
                std::cout << line << std::endl;
            }
            std::cout << "\n";
        } else {
            std::cerr << "Unable to read boundary file.\n";
            return EXIT_FAILURE;
        }

        result = getFileContent(vm["component-file"].as<std::string>(), domain_signs);

        if (result) {
            std::cout << "component domains:\n";
            for (std::string& line: domain_signs) {
                std::cout << line << std::endl;
            }
            std::cout << "\n";
        } else {
            std::cerr << "Unable to read component file.\n";
            return EXIT_FAILURE;
        }

        if (vm["manifold"].as<uint>() > 2) {
            std::cerr << "Unknown manifold option. --manifold should be 0, 1, or 2.\n";
            return EXIT_FAILURE;
        }

        const std::size_t nb_patches = boundary_files.size();

        std::map<std::pair<int, int>, double> fc_size_field_map;
        std::map<std::pair<int, int>, double> fc_length_field_map;
        std::map<int, double> cc_size_field_map;

        if (vm.count("size-field-file")) {
            std::vector<std::string> size_field_settings;
            bool result = getFileContent(vm["size-field-file"].as<std::string>(),
                                         size_field_settings);

            if (result) {
                for (std::string& line: size_field_settings) {
                    std::vector<std::string> fields;
                    boost::split(fields, line, boost::is_any_of(" "));
                    if (fields[0].compare("patch") == 0 && fields.size() == 5) {
                        int c1 = std::stoi(fields[1]);
                        int c2 = std::stoi(fields[2]);
                        if (c1 < 0 or c2 < 0 or c1 == c2) {
                            std::cerr << "Unable to define the patch with component id " << c1
                                      << " and" << c2 << "\n";
                            std::cerr
                                << "The two component ids should not be negative or the same.\n";
                            return EXIT_FAILURE;
                        }
                        double patch_size_field = std::stod(fields[3]);
                        double patch_length_field = std::stod(fields[4]);
                        auto domain_pair = c1 < c2 ? std::pair<int, int>(c1, c2)
                                                   : std::pair<int, int>(c2, c1);
                        fc_size_field_map[domain_pair] = patch_size_field;
                        fc_length_field_map[domain_pair] = patch_length_field;
                        std::cout << "Setting size field for patch between component "
                                  << domain_pair.first << " and " << domain_pair.second << ":\n";
                        std::cout << "Facet size: " << patch_size_field
                                  << " Facet length: " << patch_length_field << "\n\n";
                    } else if (fields[0].compare("component") == 0 && fields.size() == 3) {
                        int comp = std::stoi(fields[1]);
                        if (comp <= 0) {
                            std::cerr << "Unable to identify the component with id " << comp
                                      << "\n";
                            std::cerr << "The component id should be positive. Component 0 denotes "
                                         "the outer space so no sizing field is needed.\n";
                            return EXIT_FAILURE;
                        }
                        double component_size_field = std::stod(fields[2]);
                        cc_size_field_map[comp] = component_size_field;
                        std::cout << "Setting size field for component " << comp << ": "
                                  << component_size_field << "\n\n";
                    } else {
                        std::cerr << "Unknown size field setting: " << line << ".\n\n";
                        std::cerr << "patch setting example: patch 0 1 0.1 0.01\n";
                        std::cerr << "The above line sets the facet size of patch between "
                                     "component 0 and 1 to 0.1, and facet distance to 0.01\n";
                        std::cerr << "component setting example: component 1 0.1\n";
                        std::cerr << "The above line sets the cell size of component 0 to 0.1\n";
                        std::cerr << "Note that component 0 denotes the outer space, so the actual "
                                     "component id starts from 1, so patch 0 1 is on the surface "
                                     "boundary.\n\n";
                        return EXIT_FAILURE;
                    }
                }
            } else {
                std::cerr << "Unable to read size field setting file "
                          << vm["size-field-file"].as<std::string>() << ".\n";
                return EXIT_FAILURE;
            }
        }

        CGAL::Timer t;
        t.start();

        // Create domain
        std::cout << "Create domain..." << std::endl;
        std::vector<Submesh_domain*> domain_ptrs;
        Function_vector v;
        std::vector<Polyhedron> patches(nb_patches);
        CGAL::Bbox_3 bounding_box;
        bool restart_needed = false;

        for (std::size_t i = 0; i < nb_patches; ++i) {
            std::ifstream input(boundary_files[i]);
            if (!(input >> patches[i])) {
                std::cerr << "Error reading " << boundary_files[i]
                          << " as a polyhedron, try to repair the mesh.\n";
                repair(boundary_files[i]);
                restart_needed = true;
            }
            if (restart_needed) {
                continue;
            }
            Submesh_domain* d = new Submesh_domain(patches[i]);
            domain_ptrs.push_back(d);
            bounding_box += d->bbox();
            std::function<FT(Point)> func1 = [=](Point p) { return func(p, d); };
            Function f(func1);
            v.push_back(f);
        }
        if (restart_needed) {
            std::cerr << "Some boundary meshes have been repaired, please update the boundary file "
                         "list in "
                      << vm["boundary-file"].as<std::string>() << ".\n";
            return EXIT_FAILURE;
        }

        const std::size_t n_components = domain_signs.size();
        std::vector<std::string> vps;

        for (std::size_t i = 0; i < n_components; ++i) {
            int n_pos_signs = std::count(domain_signs[i].begin(), domain_signs[i].end(), '+');
            int n_neg_signs = std::count(domain_signs[i].begin(), domain_signs[i].end(), '-');
            if (n_pos_signs + n_neg_signs != nb_patches) {
                std::cerr << "Component signs should contain only + and -.\n";
                std::cerr << "The number of signs of each component should be the same as the "
                             "number of boundary meshes.\n";
                return EXIT_FAILURE;
            }
            vps.push_back(domain_signs[i]);
        }

        std::cout << "Meshing..." << std::endl;
        namespace param = CGAL::parameters;
        Mesh_domain domain(param::function = Function_wrapper(v, vps),
                           param::bounding_object = bounding_box,
                           param::relative_error_bound = 1e-6);

        Sizing_field fc_size_field(vm["fc-size"].as<double>());
        for (auto& c: fc_size_field_map) {
            fc_size_field.set_size(c.second, 2, domain.index_from_surface_patch_index(c.first));
        }

        Sizing_field fc_length_field(vm["fc-distance"].as<double>());
        for (auto& c: fc_length_field_map) {
            fc_length_field.set_size(c.second, 2, domain.index_from_surface_patch_index(c.first));
        }

        Sizing_field cc_size_field(vm["cc-size"].as<double>());
        for (auto& c: cc_size_field_map) {
            cc_size_field.set_size(c.second, 3, domain.index_from_subdomain_index(c.first));
        }

        // Set mesh criteria
        Facet_criteria facet_criteria(vm["fc-angle"].as<double>(), fc_size_field, fc_length_field);
        Cell_criteria cell_criteria(vm["cc-ratio"].as<double>(), cc_size_field);
        Mesh_criteria criteria(facet_criteria, cell_criteria);
        // Mesh generation

        CGAL::parameters::internal::Manifold_options mo;
        switch (vm["manifold"].as<uint>()) {
        case 0:
            mo = CGAL::parameters::non_manifold();
            break;

        case 1:
            mo = CGAL::parameters::manifold();
            break;

        case 2:
            mo = CGAL::parameters::manifold_with_boundary();
            break;

        default:
            mo = CGAL::parameters::non_manifold();
        }

        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(
            domain,
            criteria,
            vm["odt"].as<bool>() ? odt(time_limit = vm["odt-time"].as<double>()) : no_odt(),
            vm["lloyd"].as<bool>() ? lloyd(time_limit = vm["lloyd-time"].as<double>()) : no_lloyd(),
            vm["perturb"].as<bool>() ? perturb(time_limit = vm["perturb-time"].as<double>(),
                                               sliver_bound = vm["perturb-bound"].as<double>())
                                     : no_perturb(),
            vm["exude"].as<bool>() ? exude(time_limit = vm["exude-time"].as<double>(),
                                           sliver_bound = vm["exude-bound"].as<double>())
                                   : no_exude(),
            mo);

        // Output
        std::ofstream medit_file(output_file);
        CGAL::output_to_medit(medit_file, c3t3);

        std::cout << "Mesh has been written to " << output_file << std::endl;

        for (std::size_t i = 0; i < nb_patches; ++i) {
            delete domain_ptrs[i];
        }
        std::cout << "Time cost: " << t.time() << " sec." << std::endl;

        return EXIT_SUCCESS;
    } catch (std::exception& e) {
        std::cout << e.what() << "\n";
        return EXIT_FAILURE;
    }
}
