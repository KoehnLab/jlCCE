module DetermineMagneticAxes 

export determine_mag_axes

using LinearAlgebra
using Printf
using jlCCE
using readCIF

"""
    determine_mag_axes(system::SpinSystem)

input: SpinSystem (exported from jlCCE.jl) 

Determine the magnetic axes based on the geometry of the spin center and the nearest surrounding ligand atoms

"""
function determine_mag_axes(system::SpinSystem)
    #println("Determination of the magnetic axes")
    #println("==================================\n")
    
    # identify spin center
    if system.spin_center == "V"
        atomic_number_metal = 23
    elseif system.spin_center == "Cu"
        atomic_number_metal = 29
    elseif system.spin_center == "Pd"
        atomic_number_metal = 46
    else
        print("Error currently only V, Cu and Pd \n")
        exit()
    end

    # identify nuclear spin bath 
    if system.nuc_spin_bath == "H"
        atomic_number_nuclei = 1
    elseif system.nuc_spin_bath == "O"
        atomic_number_nuclei = 8
    elseif system.nuc_spin_bath == "N"
        atomic_number_nuclei = 7
    else 
        print("Error only proton/oxygen bath \n")
        exit()
    end
    #print("Considered nuclear spins: ",system.nuc_spin_bath,"\n")

    # call function get_spin_center: determine lattice and coordinates of the spin center
    lattice,coord_electron_spin = 
        get_spin_center(system.coord_file,atomic_number_metal,system.spin_center_index)

    #println("\n")

    # call function get_coordinates_ligand_atoms: determine coordinates of specific ligand atoms (oxygen) 
    # in the unit cell
    coords_ligand_atoms_unit_cell = 
        get_coordinates_nuclear_spins(system.coord_file,atomic_number_nuclei)

    # call get_list_ligand_atoms: determine coordinates of the ligand atoms 
    coordinates_distance_ligand_atoms,n_ligand_atoms = 
        get_bath_list(system.r_min,system.r_max_ligand_atoms,lattice,coords_ligand_atoms_unit_cell,coord_electron_spin)
    
    println("Number of considered ligand atoms: ", n_ligand_atoms)

    #print("\n")
    #println("Coordinates of ligand atoms:")

    #for i in 1:n_ligand_atoms
    #    println("O", i,": ", coordinates_distance_ligand_atoms[i]," Å")
    #end

    #print("\n")
    #println("Magnetic axes will be determined based on these coordinates. \n")


    # determine the magnetic axes from the distance coordinates of specific ligand atoms
    # middle point between two ligand atoms
    m_12 = (coordinates_distance_ligand_atoms[1] + coordinates_distance_ligand_atoms[2])/2
    m_23 = (coordinates_distance_ligand_atoms[2] + coordinates_distance_ligand_atoms[3])/2
    
    # define x,y axis 
    x_axis = m_12/norm(m_12)
    y_axis = m_23/norm(m_23)
    
    # proof of orthogonality - Gram-Schmidt process
    # first orthogonal vector - keep the x axis vector 
    u1 = x_axis
    # second orthogonal vector - orthagonalze y axis vector 
    u2 = y_axis - dot(y_axis,u1)*u1 
    y_axis = u2/norm(u2)

    # third orthogonal vector - define z axis as cross product of x, y axis
    z_axis = cross(x_axis,y_axis)/norm(cross(x_axis,y_axis))

    # define the magnetic axes matrix R_m 
    R_m =  [x_axis y_axis z_axis] 
    #println("check R_m: ", R_m)

    #if det(R_m) > 0.99
     #   println("Magnetic axes: ")
     #   @printf " x  [%10.6f %10.6f %10.6f] Å\n" R_m[1,1] R_m[2,1] R_m[3,1]
     #   @printf " y  [%10.6f %10.6f %10.6f] Å\n" R_m[1,2] R_m[2,2] R_m[3,2]
     #   @printf " z  [%10.6f %10.6f %10.6f] Å\n\n" R_m[1,3] R_m[2,3] R_m[3,3]
    #else
	#println("Determinat of the matrix of the magnetic axes is not equal to 1.")
        #R_m = [y x z]
        #if det(R_m) > 0.99
	   #println("Magnetic axes: ",)
	   #@printf " x  [%10.6f %10.6f %10.6f] Å\n" R_m[2,1] R_m[1,1] R_m[3,1]
           #@printf " y  [%10.6f %10.6f %10.6f] Å\n" R_m[2,2] R_m[1,2] R_m[3,2]
           #@printf " z  [%10.6f %10.6f %10.6f] Å\n\n" R_m[2,3] R_m[1,3] R_m[3,3]
    #end

    return R_m
end


end 





