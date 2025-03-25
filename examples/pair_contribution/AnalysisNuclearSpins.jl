push!(LOAD_PATH,"../../src")

using jlCCE
using readCIF
using DetermineMagneticAxes
using RotationMatrices

using LinearAlgebra
using Tables, CSV, DataFrames
using BenchmarkTools
using Printf
using Distances, StatsBase


# determine the coordinates of the nuclear spins wrt the spin center - using readCIF

# input parameter
coord_file = "../../../cif_files/pddbm2.pdb"
spin_center = "Pd"
atomic_number_metal = 46
spin_center_index = 1
atomic_number_nuclei = 1
r_max = 35.
r_min = 0.

# call readCIF functions to determine the coordinates of the nuclear spins wrt the spin center
lattice,coord_electron_spin = 
    get_spin_center(coord_file,atomic_number_metal,spin_center_index)
      
println("\n")
println("lattice vectors:")
@printf " a = [%20.6f %20.6f  %20.6f]\n" lattice[1][1] lattice[1][2] lattice[1][3]
@printf " b = [%20.6f %20.6f  %20.6f]\n" lattice[2][1] lattice[2][2] lattice[2][3]
@printf " c = [%20.6f %20.6f  %20.6f]\n\n" lattice[3][1] lattice[3][2] lattice[3][3]

println("Selected nucleus ",spin_center," no. ",spin_center_index," in unit cell")
println("coordinates:")
@printf " x  %20.6f Å\n" coord_electron_spin[1]
@printf " y  %20.6f Å\n" coord_electron_spin[2]
@printf " z  %20.6f Å\n\n" coord_electron_spin[3]

coords_nuclear_spins_unit_cell = 
    get_coordinates_nuclear_spins(coord_file,atomic_number_nuclei)

# call function set_supercell_list: determine the cell list
cell_list = set_supercell_list(r_max,lattice)

n_cell = size(cell_list,1)
println("r_max = ",r_max)
println("Number of replicated unit cells for spin bath: ",n_cell,"\n")

distance_coordinates_el_nucs,n_nuc= 
    get_bath_list(r_min,r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin)

#print("\n")
#println("Coordinates of the nuclear spins in Å: ", distance_coordinates_el_nucs)

# Konvertiere die Liste von Vektoren in eine Nx3-Matrix - anschließend als csv speichern
coordinates_nuclear_spins = reduce(vcat, permutedims.(distance_coordinates_el_nucs))
#print("coords: ", coordinates_nuclear_spins)
df = DataFrame(coordinates_nuclear_spins, [:x, :y, :z])  
CSV.write("coordinates_nuclear_spins_cudbm2.csv", df)

# Häufigkeit der Abstände (gerundet auf 2 Dezimalstellen)
distance_counts = countmap(round.(coordinates_nuclear_spins, digits=2))

# Abstände zwischen allen Kernspins berechnen
pairwise_distances = []

for i in 1:size(coordinates_nuclear_spins, 1)
    for j in i+1:size(coordinates_nuclear_spins, 1)  # sicherstellen, dass j > i
        push!(pairwise_distances, euclidean(coordinates_nuclear_spins[i, :], coordinates_nuclear_spins[j, :]))
    end
end

# Häufigkeit der Abstände zwischen Kernspins bestimmen
pairwise_distance_counts = countmap(round.(pairwise_distances, digits=2))

# Ergebnisse ausgeben
println("Häufigkeiten der Abstände zum Elektronenspin: ", distance_counts)
println("Häufigkeiten der Abstände zwischen Kernspins: ", pairwise_distance_counts)

# Speichern der Häufigkeiten als CSV mit zwei Spalten
df_electron = DataFrame(Distance=collect(keys(distance_counts)), Count=collect(values(distance_counts)))
df_pairwise = DataFrame(Distance=collect(keys(pairwise_distance_counts)), Count=collect(values(pairwise_distance_counts)))

CSV.write("distance_counts_to_electron.csv", df_electron)
CSV.write("pairwise_distance_counts.csv", df_pairwise)

# (Optional) Ergebnisse speichern
#CSV.write("distance_counts_to_electron.csv", DataFrame(Distance=keys(distance_counts), Count=values(distance_counts)))
#CSV.write("pairwise_distance_counts.csv", DataFrame(Distance=keys(pairwise_distance_counts), Count=values(pairwise_distance_counts)))

#print("n nuc: ",n_nuc)