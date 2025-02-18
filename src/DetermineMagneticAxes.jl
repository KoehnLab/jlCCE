module DetermineMagneticAxes 

export SpinSystem,det_mag_axes

using LinearAlgebra
using Printf

using readCIF

mutable struct SpinSystem
    # file with coordinates (anything that the file read can handle)
    coord_file::String
    # name of the atom defining the spin center ("Cu", "V", ...)
    spin_center::String
    # index within unit cell to select spin center, if name is not unique
    spin_center_index::Int
    # nuclei defining spin bath ("H", etc.)
    nuc_spin_bath::String
    # minimum interaction radius (usually 0.)
    r_min::Float64
    # maximum interaction radius
    r_max::Float64
end

SpinSystem(coord_file,spin_center,spin_center_index) = SpinSystem(
    coord_file,spin_center,spin_center_index,"H",0,10)

function det_mag_axes(system::SpinSystem)
    println("Determination of the magnetic axes starts")
    println("========================================\n")
    
    # check cif file
    print("Coordinates will be read from: ",system.coord_file,"\n")

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
    print("\n")
    print("Current metal of the spin center: ",system.spin_center,"\n")

     # identify nuclear spin bath 
     if system.nuc_spin_bath == "H"
        atomic_number_nuclei = 1
     else 
        print("Error only proton bath \n")
        exit()
     end
     
     lattice,coord_electron_spin,coords_nuclear_spins_unit_cell,coords_oxygen_unit_cell = 
        get_coordinates(system.coord_file,atomic_number_metal,system.spin_center_index,atomic_number_nuclei)

     println("\n")
     println("lattice vectors:")
     @printf " a = [%20.6f %20.6f  %20.6f]\n" lattice[1][1] lattice[1][2] lattice[1][3]
     @printf " b = [%20.6f %20.6f  %20.6f]\n" lattice[2][1] lattice[2][2] lattice[2][3]
     @printf " c = [%20.6f %20.6f  %20.6f]\n\n" lattice[3][1] lattice[3][2] lattice[3][3]

     println("Selected nucleus ",system.spin_center," no. ",system.spin_center_index," in unit cell")
     println("coordinates:")
     @printf " x  %20.6f Å\n" coord_electron_spin[1]
     @printf " y  %20.6f Å\n" coord_electron_spin[2]
     @printf " z  %20.6f Å\n\n" coord_electron_spin[3]

     # call function set_supercell_list: determine the cell list
     cell_list = set_supercell_list(system.r_max,lattice)
   
     n_cell = size(cell_list,1)
     println("r_max = ",system.r_max)
     println("Number of replicated unit cells for spin bath: ",n_cell,"\n")
   
     distance_coordinates_el_nucs,n_nuc,distance_coordinates_el_spin_oxygen = 
        get_bath_list(system.r_min,system.r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin,coords_oxygen_unit_cell)
 
     n_oxygen = size(distance_coordinates_el_spin_oxygen,1)

     println("Number of restricted oxygen: \n ")
     print(n_oxygen)

     print("\n")
     println("Restricted oxygen coordinates:")
     print(distance_coordinates_el_spin_oxygen) 
     print("\n")

     
     

     return n_oxygen,distance_coordinates_el_spin_oxygen 
end


end 





