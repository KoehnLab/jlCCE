module DetermineMagneticAxes 

export System,det_mag_axes

using LinearAlgebra
using Printf

using readCIF

mutable struct System
    # file with coordinates (anything that the file read can handle)
    coord_file::String
    # name of the atom defining the spin center ("Cu", "V", ...)
    spin_center::String
    # index within unit cell to select spin center, if name is not unique
    spin_center_index::Int
    # minimum interaction radius (usually 0.)
    r_min::Float64
    # maximum interaction radius
    r_max::Float64
end

System(coord_file,spin_center,spin_center_index) = System(
    coord_file,spin_center,spin_center_index,0,10)

function det_mag_axes(system::System)
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

     
    lattice,coord_electron_spin,coords_nuclear_spins_unit_cell,coords_oxygen_unit_cell = 
        get_coordinates(system.coord_file,atomic_number_metal,system.spin_center_index,atomic_number_nuclei,system.det_mag_axes)

    # call function set_supercell_list: determine the cell list
    cell_list = set_supercell_list(system.r_max,lattice)
   
    n_cell = size(cell_list,1)
   
    distance_coordinates_el_spin_oxygen = 
        get_bath_list(system.r_min,system.r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin,coords_oxygen_unit_cell,system.det_mag_axes)
 
    n_oxygen = size(distance_coordinates_el_spin_oxygen,1)

    #println("Number of restricted oxygen: ",n_oxygen)
    
    #print("\n")
    #println("Coordinates of oxygen:")

    #for i in 1:n_oxygen
    #    println("O", i,": ", distance_coordinates_el_spin_oxygen[i]," Å")
    #end

    #print("\n")
    #println("Magnetic axes will be determined based on these coordinates. \n")
    

    # determine the magnetic axes from the distance coordinates of oxygen
    m_12 = (distance_coordinates_el_spin_oxygen[1] + distance_coordinates_el_spin_oxygen[2])/2
    x = m_12/norm(m_12)

    m_23 = (distance_coordinates_el_spin_oxygen[2] + distance_coordinates_el_spin_oxygen[3])/2
    y = m_23/norm(m_23)
    
    #println("x 1:", x)
    #println("y 1: ",y)

    # proof of orthogonality - Gram-Schmidt process
    if abs(dot(x,y)) < 1e-10
        z = cross(x,y)/norm(cross(x,y))

	#println("z 1: ",z)
    else 
        # keep x
        u1 = x 
        # orthogonalize y to x
        u2 = y - dot(y,u1)*u1 
        y = u2/norm(u2)
        z = cross(x,y)/norm(cross(x,y))
        #println("x2: ", u1)
        #println("y2: ", u2)
        #println("z2: ", z)
    end

    R_m =  [x y z] 
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

	  #else 
	  #  R_m = [x z y]
	  #if det(R_m) > 0.99
          #else
          #    if 
          #    else
          #    end
          #end
        #end
    #end

    return R_m
end


end 





