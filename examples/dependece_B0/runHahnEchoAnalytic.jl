push!(LOAD_PATH,"../../src")

using jlCCE
using jlCCEtools
using DetermineMagneticAxes
using Tables, CSV, DataFrames
using BenchmarkTools
using rotation_B0

# set system - use the simple constructor
system = System("../cudbm2.pdb","Pd",1)
spinsystem = SpinSystem("test","Pd",1)  

# determine magnetic axes from geometry
system.det_mag_axes = true
R_m = det_mag_axes(system)
#R_m = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
spinsystem.magnetic_axes = R_m

# determine the rotated magnetic field
B0 = [0., 0., 1.]
theta = 0:5:90 
phi = 40:5:45 #0:5:90 


# modify the values (SpinSystem creates a mutable object):
spinsystem.s_el = 0.5
spinsystem.g_factor = [2.051,2.051,2.258] 
#spinsystem.g_factor = [2.051,2.051,2.051] # isotropic for the test system 
#spinsystem.B0 = [-0.748365, 0.567226, -0.343808]
spinsystem.r_max = 35.0					 
spinsystem.r_min = 0.
r_max = 35
spinsystem.r_max_bath = 10.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 50
spinsystem.use_exp = false
spinsystem.det_mag_axes = false
spinsystem.simulation_type="highfield_analytic"

# run
tab_angle_theta = Vector{Float64}()
tab_angle_phi = Vector{Float64}()
tab_Tm = Vector{Float64}()

df = DataFrame(Theta=Float64[], Phi=Float64[], Tm=Float64[])

for i in 1:size(theta)[1] 
    for j in 1:size(phi)[1]
	# determine theta for the loop
        angle_theta = deg2rad(theta[i])
        push!(tab_angle_theta, rad2deg(angle_theta))
	# determine phi for the loop 
        angle_phi = deg2rad(phi[j])
        push!(tab_angle_phi, rad2deg(angle_phi))
	
	# define theta, phi for the test system
	#angle_theta = deg2rad(0)
	#angle_phi = deg2rad(0)

	# determine rotated B0  
        B0_rot = rotation(B0,angle_theta,angle_phi,R_m)
        spinsystem.B0 = B0_rot

	# simulate intensity
        print("r_max: ",r_max,"\n")   
	time,intensity = cce(spinsystem,angle_theta,angle_phi)
        
	# determine Tm
	Tm = 2*get_decay_time(time,intensity)*10^6
	push!(tab_Tm, Tm)

	push!(df, (angle_theta, angle_phi, Tm))

	println("Simulated intensity: ",intensity,"\n") 
	println("Determined coherent time: ", Tm, "\n") 

	#CSV.write("Hahn_echo_intensity_diffB0.csv", Tables.table([angle_theta angle_phi times intensity]))
    end
end 


#CSV.write("echo.csv", Tables.table([tab_angle_theta tab_angle_phi tab_Tm]))
CSV.write("echo.csv", df)



