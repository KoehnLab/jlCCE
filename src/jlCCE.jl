module jlCCE

export SpinSystem, cce

using AtomsIO
using Unitful

# import the jl file to read cif files 
#   --> export: coordinates_nuclear_spins, distance_coordinates_el_nucs
import .readCIF

# physical constants
const hbar = 1.054571817e-34 # J s -- CODATA 2022 (rounded value)
const hbar_cgs = 1.054571817e-27 # erg s (erg = g cm^2/s^2) 

# struct to hold parameters for the run
mutable struct SpinSystem
    # file with coordinates (anything that the file read can handle)
    coord_file::String
    # name of the atom defining the spin center ("Cu", "V", ...)
    spin_center::String
    # index within unit cell to select spin center, if name is not unique
    spin_center_index::Int
    # g factor at the center (x,y,z)
    g_factor::Vector{Float64}
    # magnetic axes of spin center (column vectors for x, y, z direction)
    magnetic_axes::Matrix{Float64}
    # nuclei defining spin bath ("H", etc.)
    nuc_spin_bath::String
    # corresponding nuclear g factor
    gn_spin_bath::Float64
    # magnetic field
    B0::Vector{Float64}
    # minimum interaction radius (usually 0.)
    r_min::Float64
    # maximum interaction radius
    r_max::Float64
    # min and maxiumum of time interval
    t_min::Float64
    t_max::Float64
    # number of time steps in interval
    n_time_step::Int
end

# convenient constructor with defaults for all but the first 3 one
SpinSystem(coord_file,spin_center,spin_center_index) = SpinSystem(
    coord_file,spin_center,spin_center_index,
    [2.0,2.0,2.0],[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],"H",5.58569468,
    [0.,0.,1.],0.0,50.0,0.0,1e-3,25)

"""
    cce(system::SpinSystem)

    explanatory text goes here
"""
function cce(system::SpinSystem)
    # just some output for demo purposes (delete/modify later)
    print("Coordinates will be read from: ",system.coord_file,"\n")
    # load coordinates
    #crystal = load_system(system.coord_file)
    #lattice = ustrip.(crystal.bounding_box)
    #print("lattice vectors:\n")
    #print(lattice[1],"\n")
    #print(lattice[2],"\n")
    #print(lattice[3],"\n")

    # identify spin center
    if system.spin_center == "V"
        atomic_number_metal = 23
    elseif system.spin_center == "Cu"
        atomic_number_metal = 29
    else 
        print("Error currently only V and Cu")
        exit()
    end 

    # identify nuc spin bath 
    if system.nuc_spin_bath == "H"
        atomic_number_nuclei = 1
    else 
        print("Error only proton bath")
        exit()
    end

    # get list of spin bath nuclei
    # call function get_coordinates: determine lattice of the spin system, the coordinates of the 
        # electron spin center (x,y,z) and coordinates of the nuclear spins of the unit cell     
    lattice,coord_electron_spin,coords_nuclear_spins_unit_cell = 
        get_coordinates(system.coord_file,atomic_number_metal,atomic_number_nuclei)

    # call function set_supercell_list: determine the cell list
    cell_list = set_supercell_list(system.r_max,lattice)
    
    # call function get_bath_list: determine the distance coordinates between the electron spin center
        # and the nuclear spins 
    distance_coordinates_el_nucs = 
        get_bath_list(system.r_min,system.r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin)
    




        
    # precompute A values for electron nucleus pairs (function call)

    # initialize Intensity to 1. for all times
    time = collect(range(system.t_min,system.t_max,system.n_time_step))
    intensity = ones(system.n_time_step)
    # double loop m>n over nuclear pairs
        # call a function to:
        # compute b, c, omega for each pair, no need to store on array
        # compute v for all t values of simulation (for this m,n)
        # probabaly best: compute Intensity = Intensity .* exp(v)

    # return intensity and time
    return time, intensity
end

end