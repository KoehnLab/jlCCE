module readCIF

export get_spin_center, get_coordinates_nuclear_spins, get_coordinates, set_supercell_list, get_bath_list, get_list_ligand_atoms

using AtomsIO
using LinearAlgebra
using Unitful
#using AtomsBase


"""
    get_spin_center(cif_file::String,atomic_number_metal::Int,idx_metal::Int)

extract the relevant information from the CIF 

input: cif_file - name of the file
       atomic_number_metal - atomic number describing the metal center
       idx_metal  -  counter to select among several metal centers in the unit cell

returns: lattice - lattice vectors of the crystal
        coord_metal_atom - coordinates of the spin center (metal atom)

""" 
function get_spin_center(cif_file::String,atomic_number_metal::Int,idx_metal::Int)
    # extract relevant information from CIF
    system = load_system(cif_file)

    # coordinates in unit cell and translation vectors (lattice)
    coords_atoms_unit_cell = ustrip.(position(system))
    lattice = ustrip.(system.bounding_box)

    # atomic numbers for identifying spin centers
    atomic_n = atomic_number(system)

    # get the coordinates of the spin center on the metal atom
    idx_list = findall(x -> x==atomic_number_metal,atomic_n)
    if idx_metal > size(idx_list)[1]
        println("requested index: ", idx_metal)
        println("number of such centers in unit cell: ", size(idx_list))
        error("requested index for selecting spin center is too large")
    end
    idx_metal_atom = idx_list[idx_metal]
    coord_electron_spin = coords_atoms_unit_cell[idx_metal_atom]

    return lattice,coord_electron_spin
end


"""
    get_coordinates_nuclear_spins(cif_file::String,atomic_number_ligand_atoms::Int)

extract the coordinates of nuclear spins in the spin bath

input: cif_file - name of the file
       atomic_number_nuclei - atomic number describing the nuclei of the spin bath
                                                                    
returns: coords_nuclear_spins_unit_cell - coordinates of nuclear spins in the bath 

"""
function get_coordinates_nuclear_spins(cif_file::String,atomic_number_nuclei::Int)
    # extract relevant information from CIF
    system = load_system(cif_file)

    # coordinates in unit cell and translation vectors (lattice)
    coords_atoms_unit_cell = ustrip.(position(system))

    # atomic numbers for identifying spin centers
    atomic_n = atomic_number(system)

    # get the coordinates of all bath spins
    idx_nuclei = findall(x -> x==atomic_number_nuclei,atomic_n)
    coords_nuclear_spins_unit_cell = coords_atoms_unit_cell[idx_nuclei]

    return coords_nuclear_spins_unit_cell
end


"""
    get_coordinates(cif_file::String,atomic_number_metal::Int,idx_metal::Int,atomic_number_nuclei::Int)

extract the relevant information from the CIF

input: cif_file - name of the file
       atomic_number_metal - atomic number describing the metal center
       idx_metal  -  counter to select among several metal centers in the unit cell
       atomic_number_nuclei - atomic number describing the nuclei of the spin bath

returns: lattice - lattice vectors of the crystal
        coord_metal_atom - coordinates of the spin center (metal atom)
        coords_nuclear_spins_unit_cell - coordinates of the nuclei of the spin bath

"""
function get_coordinates(cif_file::String,atomic_number_metal::Int,idx_metal::Int,atomic_number_nuclei::Int)
    # extract relevant information from CIF
    system = load_system(cif_file)
    #print("cif file:",cif_file,"\n")
    #print("System: ",system,"\n")

    # coordinates in unit cell and translation vectors (lattice)
    coords_atoms_unit_cell = ustrip.(position(system))
    lattice = ustrip.(system.bounding_box)
    
    # unused:
    # atoms  = atomic_symbol(system)
    # atomic numbers for identifying spin centers
    atomic_n = atomic_number(system)
    #print("atomic_n:",atomic_n,"\n")
    
    # get the coordinates of the spin center on the metal atom
    idx_list = findall(x -> x==atomic_number_metal,atomic_n)
    if idx_metal > size(idx_list)[1]
        println("requested index: ", idx_metal)
        println("number of such centers in unit cell: ", size(idx_list))
        error("requested index for selecting spin center is too large")
    end
    idx_metal_atom = idx_list[idx_metal]
    coord_electron_spin = coords_atoms_unit_cell[idx_metal_atom]

    # get the coordinates of all bath spins
    idx_nuclei = findall(x -> x==atomic_number_nuclei,atomic_n)
    coords_nuclear_spins_unit_cell = coords_atoms_unit_cell[idx_nuclei]

    return lattice,coord_electron_spin,coords_nuclear_spins_unit_cell
end


#"""
#    get_coordinates_ligand_atoms(cif_file::String,atomic_number_ligand_atoms::Int)

#extract the coordinates of specific ligand atoms

#input: cif_file - name of the file
#       atomic_number_ligand_atoms - atomic number describing the atom of the ligand
                                                                    

#returns: coords_ligand_atoms_unit_cell - coordinates of the specific ligand atoms 

#"""
#function get_coordinates_ligand_atoms(cif_file::String,atomic_number_ligand_atom::Int)
    # extract relevant information from CIF
#    system = load_system(cif_file)

    # coordinates in unit cell and translation vectors (lattice)
#    coords_atoms_unit_cell = ustrip.(position(system))

    # atomic numbers for identifying spin centers
#    atomic_n = atomic_number(system)

    # get the coordinates of specific ligand atoms within the unit cell --> to determine the magnetic axes
#    idx_ligand_atoms = findall(x -> x==atomic_number_ligand_atom,atomic_n)
#    coords_ligand_atoms_unit_cell = coords_atoms_unit_cell[idx_ligand_atoms]

#    return coords_ligand_atoms_unit_cell
#end


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
function get_bath_list(r_min::Float64,r_max::Float64,lattice,coords_nuclear_spins_unit_cell,coord_spin_center)
    # call set_supercell_list here (as we only need it here)
    cell_list = set_supercell_list(r_max,lattice)

    # initialize outputs lists for spin bath
    coordinates_nuclear_spins = Vector{}[]

    n_nuc = 0

    # loop over cell list
    for Tidx in eachindex(cell_list)
        # construct the shift vector for this cell
        shift = cell_list[Tidx][1]*lattice[1] + cell_list[Tidx][2]*lattice[2] + cell_list[Tidx][3]*lattice[3]
    
        # calculate the (shifted) coordinates of all bath nuclei and their distances to the spin center
        for Aidx in eachindex(coords_nuclear_spins_unit_cell)
            # coordinates of the shifted nuclear spins
            coords_shifted_nuclear_spins = coords_nuclear_spins_unit_cell[Aidx] + shift 
            # substract the coordinates of the spin center to obtain the distance coordinates of the bath nuclei from the spin center
            coords_distance_nuclear_spins = coords_shifted_nuclear_spins - coord_spin_center
            # distances of the bath nuclei to the spin center
            distance_el_nucs = norm(coords_distance_nuclear_spins) 
         
            # restricted 
            if distance_el_nucs > r_min && distance_el_nucs <= r_max 
                push!(coordinates_nuclear_spins, coords_distance_nuclear_spins)
                n_nuc = n_nuc+1 
            end
        end
	end 

    return coordinates_nuclear_spins,n_nuc        
end
   

"""
    get_list_ligand_atoms(r_max,lattice,coord_spin_center,r_max_ligand_atoms)

get a list of specific ligand atoms with their distances from the spin center 
and another list with their (absolue) coordinates to later compute their
relative positions

"""
function get_list_ligand_atoms(r_max_ligand_atoms::Float64,lattice,coords_ligand_atoms_unit_cell,coord_spin_center)
    # call set_supercell_list here (as we only need it here)
    cell_list = set_supercell_list(r_max_ligand_atoms,lattice)
 
    # initialize outputs lists for the ligand atoms
    distance_coordinates_ligand_atoms = Vector{}[]

    for Tidx in eachindex(cell_list)
    # construct the shift vector for this cell
    shift = cell_list[Tidx][1]*lattice[1] + cell_list[Tidx][2]*lattice[2] + cell_list[Tidx][3]*lattice[3]
    # calculate the (shifted) coordinates of all oxygen and their distances to the electron spin center for Aidx in eachindex(coords_oxygen_unit_cell)
        for Aidx in eachindex(coords_ligand_atoms_unit_cell)
            # coordinates of oxygen atoms in the crystal structure
            coords_ligand_atoms = coords_ligand_atoms_unit_cell[Aidx] + shift  
            distance_coords_ligand_atoms = coords_ligand_atoms - coord_spin_center
            distance_el_ligand_atoms = norm(distance_coords_ligand_atoms)

            # restricted oxygens (nearest oxygen around the electron spin) 
            if distance_el_ligand_atoms <= r_max_ligand_atoms
            push!(distance_coordinates_ligand_atoms,distance_coords_ligand_atoms) 
            end 
        end
    end

    return distance_coordinates_ligand_atoms
end


# quick testing - remove or comment out before using as module
#cif_file = "/home/suchaneck/masterarbeit/cif_files/vopc.cif"
#atomic_number_metal = 23
#atomic_number_nuclei = 1
#r_min = 0.0
#r_max = 10.0

#lattice,coord_electron_spin,coords_nuclear_spins_unit_cell = get_coordinates(cif_file,atomic_number_metal,atomic_number_nuclei)

#cell_list = set_supercell_list(r_max,lattice)

#print("cell_list:\n")
#print(cell_list)
#print("\n")

#supercell = build_supercell(cell_list,lattice,coords_nuclear_spins_unit_cell)

#print("supercell:\n")
#print(supercell)
#print("\n")

#coordinates_nuclear_spins,distance_coordinates_el_nucs = get_bath_list(r_min,r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin)

#print("coordinates of the nuclear spins of the bath: \n")
#print(coordinates_nuclear_spins)

#print("distance coordinates between the spin center and nuclear bath spins: \n")
#print(distance_coordinates_el_nucs)

end
