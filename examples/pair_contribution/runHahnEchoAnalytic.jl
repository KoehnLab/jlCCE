### PAIR CONTRIBUTION ####

push!(LOAD_PATH,"../../src")

using jlCCE
using jlCCEtools
using DetermineMagneticAxes
using RotationMatrices
using Tables, CSV, DataFrames
using BenchmarkTools
using Printf

#spinsystem = SpinSystem("../cudbm2.pdb","Pd",1)  

# determine magnetic axes from geometry
spinsystem = SpinSystem("../../../cif_files/pddbm2.pdb","Pd",1,"O",2.0)
R_m = determine_mag_axes(spinsystem)
println("R_m: ",R_m)

# cce settings 
spinsystem = SpinSystem("../../../cif_files/pddbm2.pdb","Pd",1)

# determine the rotated magnetic field
B0 = [0., 0., 1.]
theta = 0.
phi = 0.
rot_mat = rotate_solid(deg2rad(theta),deg2rad(phi))
spinsystem.B0 = R_m * (rot_mat * B0) 
#println("B0: ",spinsystem.B0)

# modify the values (SpinSystem creates a mutable object):
spinsystem.magnetic_axes = R_m
spinsystem.s_el = 0.5
spinsystem.g_factor = [2.051,2.051,2.258] # anisotropic 
#spinsystem.g_factor = [2.12,2.12,2.12] # g_iso
#spinsystem.r_max = 5.:0.2:35.
r_max = 3.8:0.2:35.
spinsystem.r_min = 0.
spinsystem.r_max_bath = 10.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 50
spinsystem.use_exp = false
spinsystem.simulation_type="highfield_analytic"

# run
tab_angle_theta = Vector{Float64}()
tab_angle_phi = Vector{Float64}()
tab_Tm = Vector{Float64}()

df = DataFrame(theta=Float64[], phi=Float64[], Tm=Float64[])

for i in 1:size(r_max)[1] 
    #for j in 1:size(phi)[1]
        spinsystem.r_max = r_max[i]

        rot_mat = rotate_solid(deg2rad(theta[i]),deg2rad(phi[j]))
		spinsystem.B0 = (rot_mat * B0)
		push!(tab_angle_theta, theta[i])
        push!(tab_angle_phi, phi[j])
	
		# simulate intensity
        #print("r_max: ",r_max,"\n")   
		time,intensity = cce(spinsystem)


		# determine Tm
		Tm = 2*get_decay_time(time,intensity)*10^6
		push!(tab_Tm, Tm)

		push!(df, (theta[i], phi[j], Tm))

        println("Maximum distance between spin center and nuclear spins: ", r_max[i],"Å \n")

		println("Simulated intensity: ",intensity,"\n") 
		println("Determined coherence time: ", Tm, "\n") 
    #end
end 


#CSV.write("echo.csv", Tables.table([tab_angle_theta tab_angle_phi tab_Tm]))
CSV.write("echo_cudbm2_anisotropic_360_degree_B0_cartesian_coords.csv", df)



