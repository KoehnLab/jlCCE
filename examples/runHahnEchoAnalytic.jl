push!(LOAD_PATH,"../src")

using jlCCE
using jlCCEtools
using DetermineMagneticAxes
using RotationMatrices
using Tables, CSV
using BenchmarkTools

# set system - use the simple constructor
#spinsystem = SpinSystem("cudbm2.pdb","Pd",1)  # <--- note: put this file into the folder

# determine magnetic axes from geometry
spinsystem = SpinSystem("cudbm2.pdb","Pd",1,"O",2.0)
R_m = determine_mag_axes(spinsystem)
#println("R_m: ",R_m)
spinsystem.magnetic_axes = R_m

# define the magnetic field
B0 = [0., 0., 1.]
theta = 0.
phi = 0.
rot_mat = rotate_solid(deg2rad(theta),deg2rad(phi))
spinsystem.B0 = R_m * (rot_mat * B0)

# define the spinsystem for thr run of the simulation
spinsystem = SpinSystem("cudbm2.pdb","Pd",1)
spinsystem.s_el = 0.5
spinsystem.g_factor = [2.051,2.051,2.258] 
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
times,intensity = cce(spinsystem)
print("\n")
print("Simulated intensity: ",intensity,"\n")
print("\n")
#println("Simulation time: ", times)
# determine the coherence time Tm
Tm = 2*get_decay_time(times,intensity)*10^6
println("Simulated coherence time: ",Tm,"\n")



#CSV.write("echo.csv", Tables.table([times intensity intensityNE intensityEX iCCE1 iCCE2]))


