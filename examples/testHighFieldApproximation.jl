push!(LOAD_PATH,"../src")

using jlCCE
#      vvvvvvvvvv  you might need to install these in addition
using Tables, CSV
using DetermineMagneticAxes
using RotationMatrices
#using BenchmarkTools
using Printf
using LinearAlgebra

BLAS.set_num_threads(Threads.nthreads())

@printf "Threads set for BLAS:  %3i\n" BLAS.get_num_threads()
@printf "Threads set for Julia: %3i\n" Threads.nthreads()

# set system - use the simple constructor
spinsystem = SpinSystem("cudbm2.pdb","Pd",1)  # <--- note: Pd is the host crystal; it's actually Cu = Spin 1/2

# determine magnetic axes from geometry
R_m = determine_mag_axes(spinsystem,["O"],3.0)

# define the magnetic field -- 100 Tesla (really high field)
B0 = [0.,0.,100.]*10000 # 10 kGaus = 1 Tesla
# orientation of B0 relative to xy plane of magnetic axis system: in plane
theta = 90.
phi = 0.
rot_mat = rotate_solid(deg2rad(theta),deg2rad(phi))
spinsystem.B0 = R_m * (rot_mat * B0)  

# modify the values (SpinSystem creates a mutable object):
spinsystem.g_factor = [2.0200,2.0200,2.0837]   # anisotropic
spinsystem.magnetic_axes = R_m

# defaults, no need to define here
#spinsystem.nuc_spin_bath = "H"
#spinsystem.gn_spin_bath = 5.58569468
#spinsystem.use_exp = false

spinsystem.r_min = 0.
spinsystem.r_max = 30.
spinsystem.r_max_bath = 10.


spinsystem.t_min = 0.
spinsystem.t_max = 20e-6 # 20 µs (delay time)
spinsystem.n_time_step = 21

println("Running simulation for 100 T")

# do analytic high-field approx
@time times,intensity = cce(spinsystem)
    
# and now the "exact" solutions
spinsystem.simulation_type="exact"
    
@time times,intensityEX,iCCE0,iCCE1,iCCE2 = cce(spinsystem)

@printf "    2*tau/µs    High-Field         Exact\n"
for (tau,intHFA,intEX) in zip(times,intensity,intensityEX)
   @printf " %12.6f  %12.6e  %12.6e\n" 2*tau*1e6 intHFA intEX
end
println("")

CSV.write("cu_in_pd-dbm2-echo-100T.csv", Tables.table([times intensity intensityEX iCCE0 iCCE1 iCCE2]))
 
println("All done ...")
