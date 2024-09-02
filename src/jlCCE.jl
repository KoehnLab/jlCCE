module jlCCE

#push!(LOAD_PATH,".")

export SpinSystem, cce

using AtomsIO
using Unitful
using LinearAlgebra

# import the jl file to read cif files 
#   --> export: distance_coordinates_el_nucs
using readCIF

# physical constants
const hbar = 1.054571817e-27 # erg s -- J s -- CODATA 2022 (rounded value)
# Bohr magneton 
const mu_b = 9.274010066e-21 # erg G^-1 J T^-1 -- CODATA 2022 (rounded value)
# nuclear magneton
const mu_n = 5.050783739e-24 # erg G^-1  J T^-1 -- CODATA 2022 (rounded value)

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
    # minimum and maxiumum of time interval
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

input: Spinsystem - cif file of the spin system, name of the metal spin center,
        index of the metal atom, g factor of the central spin, magnetic axes of 
        the spin center, nuclei defining spin center, corresponding nuclear g 
        factor, magnetic field, maximun und manimum radius around the spin center


    explanatory text goes here
"""
function cce(system::SpinSystem)
    # check cif file
    print("Coordinates will be read from: ",system.coord_file,"\n")

    # identify spin center
    if system.spin_center == "V"
        atomic_number_metal = 23
    elseif system.spin_center == "Cu"
        atomic_number_metal = 29
    else 
        print("Error currently only V and Cu")
        exit()
    end 

    # identify nuclear spin bath 
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

    print("lattice vectors:\n")
    print(lattice[1],"\n")
    print(lattice[2],"\n")
    print(lattice[3],"\n")

    # call function set_supercell_list: determine the cell list
    cell_list = set_supercell_list(system.r_max,lattice)

    print("cell list: \n")
    print(cell_list,"\n")

    # call function get_bath_list: determine the distance coordinates between the electron spin center
        # and the nuclear spins 
    distance_coordinates_el_nucs,n_nuc = 
        get_bath_list(system.r_min,system.r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin)
    # rescale distance coords to cm 
    distance_coordinates_el_nucs = distance_coordinates_el_nucs .* 1e-8

    print("Number of bath nuclei: ",n_nuc,"\n")
    print("Distance coordinates between the electron spin center and the nuclear spins: \n")
    print(distance_coordinates_el_nucs,"\n") 
    

    # calculation of the gryomagnetic ratio 
    gamma_electron = (system.g_factor[1] * mu_b) / hbar
    gamma_n = (system.gn_spin_bath * mu_n) / hbar
    print("gamma el: ",gamma_electron,"\n")
    print("gamma n: ",gamma_n,"\n")

    # precompute A values for electron nucleus pairs
    A_n = zeros(n_nuc)
    for i in 1:size(distance_coordinates_el_nucs)[1]
        r_i_x_B0 = cross(distance_coordinates_el_nucs[i], system.B0)
        r_i_dot_B0 = dot(distance_coordinates_el_nucs[i], system.B0)
        theta_i = atan(norm(r_i_x_B0), r_i_dot_B0)
        r_i_norm = norm(distance_coordinates_el_nucs[i])
        print("Distance between el spin and nuc spins: ", r_i_norm,"\n") 
        #print(gamma_n,gamma_electron,hbar,theta_i,r_i_norm)
        A_n[i] = -gamma_n * gamma_electron * hbar * (1 - 3 * cos(theta_i)^2) / r_i_norm^3
   end
   print("A_n:",A_n,"\n")

    # initialize Intensity to 1. for all times
    time = collect(range(system.t_min,system.t_max,system.n_time_step))
    intensity = ones(system.n_time_step)
    # double loop m>n over nuclear pairs
        # call a function to:
        # compute b, c, omega for each pair, no need to store on array
        # compute v for all t values of simulation (for this m,n)
        # probabaly best: compute Intensity = Intensity .* exp(v)

    #coordinates_nuclear_spins = Vector{}[] 
    #for i in eachindex(distance_coordinates_el_nucs)
     #   coords_nuclear_spins = distance_coordinates_el_nucs[i] + coord_electron_spin
     #   push!(coordinates_nuclear_spins,coord_electron_spin)
    #end 
    #print("coord nuc spins: ",coordinates_nuclear_spins,"\n")
    
    # loop over all pairs of bath nuclei
    for n in 1:n_nuc-1
        for m in n+1:n_nuc
            r_nm = (distance_coordinates_el_nucs[m] - distance_coordinates_el_nucs[n]) 
            #print("coord nuc m: ",coordinates_nuclear_spins[m],"\n")
            #print("coord nuc n: ",coordinates_nuclear_spins[n],"\n")
            r_nm_x_B0 = cross(r_nm, system.B0)
            r_nm_dot_B0 = dot(r_nm, system.B0)
            theta_nm = atan(norm(r_nm_x_B0), r_nm_dot_B0)
            b_nm = -0.25 * gamma_n^2 * hbar * (1 - 3 * cos(theta_nm)^2) / norm(r_nm)^3
            c_nm = (A_n[n] - A_n[m]) / (4 * b_nm)
            w_nm = 2 * b_nm * sqrt(1 + c_nm^2)
            #print("r_nm: ",r_nm,"\n")
            #print("b_nm: ",b_nm,"\n")
            #print("c_nm: ",c_nm,"\n")
            #print("w_nm: ",w_nm,"\n")
            for j in 1:size(time)[1]
                v_nm = -((c_nm^2) / (1 + c_nm^2)^2) * (cos(w_nm * time[j]) - 1)^2
                #print("v_nm: ",v_nm,"\n")
                    #sim[j] = sum(v_nm) 
                intensity[j] = intensity[j] * exp(v_nm)
            end 
        end
    end        



    # return intensity and time
    return time, intensity
end

end