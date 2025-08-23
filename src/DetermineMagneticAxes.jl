module DetermineMagneticAxes 

export determine_mag_axes

using LinearAlgebra
using Printf
using jlCCE
using readCIF

"""
    determine_mag_axes(system::SpinSystem,coordinating_atoms::Vector{String},cutoff::Float64)

input: SpinSystem (exported from jlCCE.jl), list of atom symbols (coordinating atoms of ligand), cutoff radius 

Determine the magnetic axes based on the geometry of the spin center and the nearest surrounding ligand atoms

"""
function determine_mag_axes(system::SpinSystem,coordinating_atoms::Vector{String},cutoff::Float64)
    println("")
    println("Determination of the magnetic axes")
    println("==================================\n")
    
    # identify spin center
    atomic_number_metal = symbol_to_atomic_number(system.spin_center)
    
    # call function get_spin_center: determine lattice and coordinates of the spin center
    lattice,coord_electron_spin = 
        get_spin_center(system.coord_file,atomic_number_metal,system.spin_center_index)

    # call function get_coordinates_ligand_atoms: determine coordinates of specific ligand atoms (oxygen) 
    # in the unit cell
    coords_ligand_atoms_unit_cell = []
    for (idx,coord_atom) in enumerate(coordinating_atoms)
        atomic_number_atom = symbol_to_atomic_number(coord_atom)
        coords_ligand_atoms_unit_cell = vcat(coords_ligand_atoms_unit_cell,
            get_coordinates_nuclear_spins(system.coord_file,atomic_number_atom))
    end

    # call get_list_ligand_atoms: determine coordinates of the ligand atoms 
    coordinates_distance_ligand_atoms,n_ligand_atoms = 
        get_bath_list(system.r_min,cutoff,lattice,coords_ligand_atoms_unit_cell,coord_electron_spin)
    
    println("Number of considered ligand atoms: ", n_ligand_atoms)
    print("\n")
    println("Coordinates of ligand atoms:")

    for i in 1:n_ligand_atoms
        @printf "Ligand atom %1i: [ %10.6f   %10.6f   %10.6f ] Å \n" i coordinates_distance_ligand_atoms[i][1] coordinates_distance_ligand_atoms[i][2] coordinates_distance_ligand_atoms[i][3]
    end

    if n_ligand_atoms != 4
        println("Does currently only work for 4 (near) coplanar atoms ")
        error("Not general enough routine")
    end

    if ! test_coplanar(coordinates_distance_ligand_atoms[1],coordinates_distance_ligand_atoms[2],
                       coordinates_distance_ligand_atoms[3],coordinates_distance_ligand_atoms[4], tol=0.1)
        println("Does currently only work for planar arrangements")
        # not sure, we can make a sensible descision in too many other cases anyway ...
        # we will need the computed g tensor for that
        error("Non coplanar case")
    end

    # determine the magnetic axes from the distance coordinates of specific ligand atoms
    # middle point between two ligand atoms
    # (in principle, we need two arbitrary vectors spanning the plane; the approach will only 
    #  give sensible magnetic axes if the perpendicular part of g is degenerate)
    m_12 = (coordinates_distance_ligand_atoms[1] + coordinates_distance_ligand_atoms[2])/2
    m_23 = (coordinates_distance_ligand_atoms[2] + coordinates_distance_ligand_atoms[3])/2
    
    # normalize
    x_axis = m_12/norm(m_12)
    y_axis = m_23/norm(m_23)
    
    # ensure orthogonality - Gram-Schmidt process
    # first orthogonal vector - keep the x axis vector 
    u1 = x_axis
    # second orthogonal vector - orthogonalze y axis vector 
    u2 = y_axis - dot(y_axis,u1)*u1 
    y_axis = u2/norm(u2)

    # third orthogonal vector - define z axis as cross product of x, y axis
    z_axis = cross(x_axis,y_axis)/norm(cross(x_axis,y_axis))

    # define the magnetic axes matrix R_m 
    R_m =  [x_axis y_axis z_axis] 

    # should not happen, as the cross product should guarantee a right-handed system
    if det(R_m) < 0.
        println("Something is wicked ... det < 0, check code!")
        error("Problem with magnetic axis determination")
    end

    println("Magnetic axes: ")
    @printf " x  [%10.6f %10.6f %10.6f] Å\n" R_m[1,1] R_m[2,1] R_m[3,1]
    @printf " y  [%10.6f %10.6f %10.6f] Å\n" R_m[1,2] R_m[2,2] R_m[3,2]
    @printf " z  [%10.6f %10.6f %10.6f] Å\n\n" R_m[1,3] R_m[2,3] R_m[3,3]

    return R_m
end


# Test coplanarity of 4 points in 3D
# Each point should be a 3-element vector (e.g. [x,y,z])
function test_coplanar(A::AbstractVector, B::AbstractVector, 
                       C::AbstractVector, D::AbstractVector; 
                       tol::Float64=1e-8)
    v1 = B .- A
    v2 = C .- A
    v3 = D .- A
    volume6 = dot(v1, cross(v2, v3))   # 6 * volume of tetrahedron
    return abs(volume6) < tol
end


end 





