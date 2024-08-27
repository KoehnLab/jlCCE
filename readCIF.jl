module readCIF

export get_coordinates, set_supercell_list, get_bath_list

using AtomsIO
using LinearAlgebra
using Unitful

"""
    get_coordinates(cif_file::String,atomic_number_metal::Int,atomic_number_nuclei::Int)

extract the relevant information from the CIF

input: cif_file - name of the file
       atomic_number_metal - atomic number describing the metal center
       atomic_number_nuclei - atomic number describing the nuclei of the spin bath

returns: lattice - lattice vectors of the crystal
        coord_metal_atom - coordinates of the spin center (metal atom)
        coords_nuclear_spins_unit_cell - coordinates of the nuclei of the spin bath

"""
function get_coordinates(cif_file::String,atomic_number_metal::Int,atomic_number_nuclei::Int)
    # extract relevant information from CIF
    system = load_system(cif_file)

    # coordinates in unit cell and translation vectors (lattice)
    coords_atoms_unit_cell = ustrip.(position(system))
    lattice = ustrip.(system.bounding_box)
    
    # unused:
    # atoms  = atomic_symbol(system)
    # atomic numbers for identifying spin centers
    atomic_n = atomic_number(system)
    
    # get the coordinates of the spin center on the metal atom
    idx_metal_atom = findfirst(x -> x==atomic_number_metal,atomic_n)
    coord_electron_spin = coords_atoms_unit_cell[idx_metal_atom]

    # get the coordinates of all bath spins
    idx_nuclei = findall(x -> x==atomic_number_nuclei,atomic_n)
    coords_nuclear_spins_unit_cell = coords_atoms_unit_cell[idx_nuclei]

    return lattice,coord_electron_spin,coords_nuclear_spins_unit_cell
end

"""
    set_supercell_list(r_max,lattice)

get the size of the supercell such that r_max lies fully in it
simple version without considering the placement of the origin in
the central unit cell

input: r_max   - radius around spin center
       lattice - a vector of three vectors giving the unit cell

returns: cell_list - a list of all relevant translations of the unit cell
       to build the supercell
"""
function set_supercell_list(r_max::Float64,lattice)
    
    # compute how man times the lattice in x/y/z directions goes into r_max
    # use the ceil function to get the upper bound of that value as integer
    nx = ceil(Int,r_max/sqrt(dot(lattice[1],lattice[1])))
    ny = ceil(Int,r_max/sqrt(dot(lattice[2],lattice[2])))
    nz = ceil(Int,r_max/sqrt(dot(lattice[3],lattice[3])))

    # set up a list of all required neighboring unit cells
    cell_list = Vector{Int}[]
    for ii = -nx:nx
        for jj = -ny:ny
            for kk = -nz:nz
                push!(cell_list,[ii,jj,kk])
            end
        end
    end
    return cell_list
end

#"""
#    build_supercell(cell_list,lattice,coord_unit_cell)

#build the complete supercell (of spin lattice)
#note: we actually do not need this function, except for visualization
#"""
# function build_supercell(cell_list::Matrix{Int},lattice,coord_unit_cell)

    # initialize the vector containing the coordinates
#    supercell = Vector{}[]

    # loop over cell list
#    for Tidx in eachindex(cell_list)
	# construct the shift vector for this cell
#        shift = cell_list[Tidx][1]*lattice[1] + cell_list[Tidx][2]*lattice[2] + cell_list[Tidx][3]*lattice[3]

	# push the shifted coordinates of the unit cell to the supercell vector
#	for Aidx in eachindex(coord_unit_cell)
#	    push!(supercell, coord_unit_cell[Aidx] + shift)
#        end
    
#    end
#    return supercell
# end

"""
    get_bath_list(r_min,r_max,lattice,coords_spins_unit_cell,coord_spin_center)

get a list of all bath nuclei with their distances from the spin_center
and another list with their (absolute) coordinates to later compute their
relative positions

"""
function get_bath_list(r_min::Float64,r_max::Float64,lattice,coords_spins_unit_cell,coord_spin_center)

    # call set_supercell_list here (as we only need it here)
    cell_list = set_supercell_list(r_max,lattice)

    # initialize outputs lists for spin bath
    coordinates_nuclear_spins = Vector{}[]
    distance_coordinates_el_nucs = Vector{}[]
    distance_coordinates_nuc_nuc = Vector{}[]

    # loop over cell list
    for Tidx in eachindex(cell_list)
        # construct the shift vector for this cell
        shift = cell_list[Tidx][1]*lattice[1] + cell_list[Tidx][2]*lattice[2] + cell_list[Tidx][3]*lattice[3]
    
        # calculate the (shifted) coordinates of all bath nuclei and their distances to the spin center
        for Aidx in eachindex(coords_spins_unit_cell)
            # coordinates of the shifted nuclear spins
            coords_nuclear_spins = coords_spins_unit_cell[Aidx] + shift 
            # substract the coordinates of the spin center to obtain the distance coordinates of the bath nuclei from the spin center
            distance_coords_el_nucs = coords_nuclear_spins - coord_electron_spin
            # distances of the bath nuclei to the spin center
            distance_el_nucs = norm(distance_coords_el_nucs) 
         
            # restricted 
            if distance_el_nucs > r_min && distance_el_nucs <= r_max
                push!(coordinates_nuclear_spins,coords_nuclear_spins) 
                push!(distance_coordinates_el_nucs, distance_coords_el_nucs)
            end
        end
        
    end
    return coordinates_nuclear_spins,distance_coordinates_el_nucs

end
   
# quick testing - remove or comment out before using as module
cif_file = "/home/suchaneck/masterarbeit/cif_files/vopc.cif"
atomic_number_metal = 23
atomic_number_nuclei = 1
r_min = 0.0
r_max = 10.0

lattice,coord_electron_spin,coords_nuclear_spins_unit_cell = get_coordinates(cif_file,atomic_number_metal,atomic_number_nuclei)

cell_list = set_supercell_list(r_max,lattice)

print("cell_list:\n")
print(cell_list)
print("\n")

#supercell = build_supercell(cell_list,lattice,coords_nuclear_spins_unit_cell)

#print("supercell:\n")
#print(supercell)
#print("\n")

coordinates_nuclear_spins,distance_coordinates_el_nucs = get_bath_list(r_min,r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin)

print("coordinates of the nuclear spins of the bath: \n")
print(coordinates_nuclear_spins)

print("distance coordinates between the spin center and nuclear bath spins: \n")
print(distance_coordinates_el_nucs)

end
