push!(LOAD_PATH,"../src")

using jlCCE
using jlCCEtools
using DetermineMagneticAxes
using RotationMatrices
using Tables, CSV
using BenchmarkTools
using Printf

# determine magnetic axes from geometry
spinsystem = SpinSystem("cudbm2.pdb","Pd",1)

R_m = determine_mag_axes(spinsystem,["O"],3.0)

# define the magnetic field
B0 = [1., 0., 0.]
theta = 0.
phi = 0.
rot_mat = rotate_solid(deg2rad(theta),deg2rad(phi))
spinsystem.B0 = R_m * (rot_mat * B0)
#println("B0: ",spinsystem.B0)

# define the spinsystem for thr run of the simulation
spinsystem.magnetic_axes = R_m
spinsystem.s_el = 0.5
spinsystem.g_factor = [2.0837,2.0200,2.0200]
spinsystem.r_max = 35.0					 
spinsystem.r_min = 0.
spinsystem.r_max_bath = 15.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 11
spinsystem.use_exp = false
spinsystem.simulation_type="highfield_analytic"
spinsystem.report_pair_contrib = false
spinsystem.pair_log_file = "pair_log.txt"

# run
# run CCE - determine the Hahn echo intensity
times,intensity = cce(spinsystem)
print("\n")
print("      tau/µs      2*tau/µs        s(2*tau)\n")
for (tau,sHE) in zip(times,intensity)
    @printf " %12.6f  %12.6f   %12.6e\n" tau*1e6 2*tau*1e6 sHE
end
print("\n")
# determine the coherence time Tm
Tm = 2*get_decay_time(times,intensity)*10^6
@printf "Simulated coherence time: %12.6f µs\n\n" Tm



#CSV.write("echo.csv", Tables.table([times intensity intensityNE intensityEX iCCE1 iCCE2]))


