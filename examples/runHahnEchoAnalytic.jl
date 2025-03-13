push!(LOAD_PATH,"../src")

using jlCCE
using Tables, CSV
using BenchmarkTools

# set system - use the simple constructor
spinsystem = SpinSystem("../cudbm2.pdb","Pd",1)  # <--- note: put this file into the folder

# determine magnetic axes from geometry
R_m = determine_mag_axes(spinsystem)
println("R_m: ",R_m)
spinsystem.magnetic_axes = R_m

# define the magnetic field
B0 = [0., 0., 1.]
theta = 0.
phi = 0.
rot_mat = rotate_solid(deg2rad(theta[i]),deg2rad(phi[j]))
spinsystem.B0 = R_m * (rot_mat * B0)


# modify the values (SpinSystem creates a mutable object):
spinsystem.s_el = 0.5
spinsystem.g_factor = [2.051,2.051,2.258] 
spinsystem.r_max = 35.0					 
spinsystem.r_min = 0.
r_max = 35
spinsystem.r_max_bath = 10.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 50
spinsystem.use_exp = false
spinsystem.simulation_type="highfield_analytic"

# run
# run CCE - determine the Hahn echo intensity 
times,intensity = cce(spinsystem)
# determine Tm
Tm = 2*get_decay_time(time,intensity)*10^6
#print("\n")
#print("Simulated intensity: ",intensity,"\n")
#print("\n")


#CSV.write("echo.csv", Tables.table([times intensity intensityNE intensityEX iCCE1 iCCE2]))


