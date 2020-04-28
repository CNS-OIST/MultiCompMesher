# MultiCompMesher

This program utilizes and extends the [CGAL 3D mesh generation functionality](https://doc.cgal.org/latest/Mesh_3/index.html) to
generate element-labelled multi-component tetrahedral mesh.


This program reads a list of watertight boundary surface meshes,
and a sign file descripting the spatial relationship between each 
desired component and the boundary meshes, then generates a tetrahedral mesh with each component labelled separately according to the boundaries
and the component signs.

# Prerequisite
* git
* a C++ compiler
* [CMake](https://cmake.org/) 
* [CGAL](https://www.cgal.org/) 
* [Boost](https://www.boost.org/) 

# Compilation
```
git clone https://github.com/CNS-OIST/MultiCompMesher.git
cd MultiCompMesher
mkdir build
cd build
cmake ..
make
```

# Usage
Below commands assume you are still in `MultiCompMesher/build`
    
* Basic usage
    ```
    ./MultiCompMesher BOUNDARY-FILE COMPONENT-FILE OUTPUT
    ```
    * BOUNDARY-FILE: a plain text file where each line stores
    the path to a boundary surface mesh file (in `.off` format)
    * COMPONENT-FILE: a plain text file where each line labels
    the relationship between each component and the boundaries
    using `+` (exclusive) or `-` (inclusive) signs.  
    For example,  assume that there are three boundaries, b1, b2 
    and b3. A component labelled `-++` means the tetrahedrons of the 
    component should be in b1, and should not be in b2 and b3. 
    A component labelled `+--` means the tetrahedrons of the component 
    should not be in b1, but should be in both b2 and b3.
    * OUTPUT: Path to the output mesh file. The output mesh is in .mesh
    format, which can be opened and converted to other formats in 
    [Gmsh](http://gmsh.info/).

* Advance usage  
    Other parameters can be set to control the meshing process 
    and optimize the mesh quality. You can list them using
    ```
    ./MultiCompMesher -h
    ```
    The usage of these parameters can be found in the 
    [CGAL 3D mesh generation manual](https://doc.cgal.org/latest/Mesh_3/index.html).

    If there is problem with one of more boundary meshes, such as orintation issue, 
    this program will try to repair the mesh and save
    the repaired version to another file. In this case, please change the
    corresponding file path in your `BOUNDARY-FILE` and rerun the program.

# Example

This example generates a mesh modelling a dendritic spine, with three
labelled components, the spine cytosol (excluding other two components), 
the Endoplasmic Reticulum (ER), and the Postsynaptic Density (PSD) region.

The [boundary file](example/boundaries.txt) stores the path to three
boundary surfaces.  
1. [spine boundary](example/Spine.off)
2. [ER boundary](example/ER.off)
3. [PSD boundary](example/PSD.off)

The [component file](example/components.txt) defines three components.  
1. The spine cytosol excluding ER and PSD: `-++`
2. The ER: `--+`
3. The PSD: `-+-`

The command to generate the mesh is
```
# assume still in MultiCompMesher/build

./MultiCompMesher ../example/boundaries.txt ../example/components.txt ../example/output --fc-size 0.01 --fc-distance 0.001 --cc-size 0.05 --odt --lloyd --perturb --exude
```
The mesh is written to [output.mesh](example/output.mesh), then visualized
in Gmsh. Note that each component are labelled and colorred individually. To use the mesh in [STEPS](http://steps.sourceforge.net), the user need to
export it to the Abaqus inp format or the Gmsh 2.0 ASCii format in Gmsh.
![Mesh visualization in Gmsh](example/mesh_view.png)
