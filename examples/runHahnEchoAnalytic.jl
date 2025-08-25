push!(LOAD_PATH,"../src")

using jlCCE
using jlCCEtools
using DetermineMagneticAxes
using RotationMatrices
using Tables, CSV
using BenchmarkTools

# determine magnetic axes from geometry
spinsystem = SpinSystem("../../cif_files/pddbm2.pdb","Pd",1,"O",2.0)
spinsystem.r_min = 0.
R_m = determine_mag_axes(spinsystem)
#println("R_m: ",R_m)
#spinsystem.magnetic_axes = R_m


# cce settings
spinsystem = SpinSystem("../../cif_files/pddbm2.pdb","Pd",1)

# define the magnetic field
B0 = [0.,0.,1.]
theta = 0.
phi = 0.
rot_mat = rotate_solid(deg2rad(theta),deg2rad(phi))
spinsystem.B0 = R_m * (rot_mat * B0)
#println("B0: ",spinsystem.B0)

# define the spinsystem for thr run of the simulation
spinsystem.magnetic_axes = R_m
spinsystem.s_el = 0.5
#spinsystem.g_factor = [1.99,1.99,2.00]
#spinsystem.g_factor = [2.0379 ,2.0379,2.150]
#spinsystem.g_factor = [2.0424 ,2.0554,2.225]
spinsystem.g_factor = [2.051,2.051,2.258] # cudbm2 
spinsystem.r_max = 35.0					 
spinsystem.r_min = 0.
spinsystem.r_max_bath = 10.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 50
spinsystem.use_exp = false
spinsystem.simulation_type="highfield_analytic"

# run
# run CCE - determine the Hahn echo intensity
times,R,r_12,intensity = cce(spinsystem)
print("\n")
print("Simulated intensity: ",intensity,"\n")
print("\n")
# determine the coherence time Tm
Tm = 2*get_decay_time(times,intensity)*10^6
println("Simulated coherence time: ",Tm,"\n")

CSV.write("distances_cudbm2.csv", Tables.table([R r_12])) 

#CSV.write("HEsignal_cudbm2.csv", Tables.table([times intensity]))


