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

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/IO/OFF.h>
#include <iostream>
#include <string>
#include <vector>

#include "utility.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> PolyhedronWithID;

int repair(std::string& filename) {
    std::ifstream input(filename);

    std::vector<K::Point_3> points;
    std::vector<std::vector<std::size_t>> polygons;

    if (!input || !CGAL::IO::read_OFF(input, points, polygons) || points.empty()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        return EXIT_FAILURE;
    }
    CGAL::Polygon_mesh_processing::repair_polygon_soup(points, polygons);
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

    PolyhedronWithID mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

    // Number the faces because 'orient_to_bound_a_volume' needs a face <--> index map
    int index = 0;
    for (PolyhedronWithID::Face_iterator fb = mesh.facets_begin(), fe = mesh.facets_end(); fb != fe;
         ++fb)
        fb->id() = index++;

    if (CGAL::is_closed(mesh))
        CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(mesh);

    string_replace(filename, ".off", "_repaired.off");
    std::ofstream out(filename);
    out.precision(17);
    if (!(out << mesh)) {
        std::cerr << "Unable to write back repaired mesh data to " << filename << ".\n";
        return EXIT_FAILURE;
    }
    std::cout << "Repaired mesh data has been written to " << filename << ".\n";
    return EXIT_SUCCESS;
}
