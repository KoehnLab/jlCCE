# run ParameterSpaceScreening

push!(LOAD_PATH,"../../src")

using LinearAlgebra
using Tables, CSV, DataFrames
using Printf
using ParameterSpaceScreening 
using RotationMatrices

# set parameters - R,R_12,alpha,ge,gn,time,B0
#R = collect(range(1.,35.,150))
#R_12 = collect(range(1.,10.,50))
alpha = deg2rad(45)
#alpha_list = collect(range(0.,pi,5))

# scan of alpha for fixed distances R,R_12
R_12 = 35.
R = 35.
#alpha_list = 0.:2:360.

# time-dependent simulation
#time_list = (101.:1.:150.) .* 1e-6

g_e = 2.002;
g_n = 5.586;
time = 10e-6
B0_initial = [0., 0., 1.];
theta = 0.:1:180.
phi = 0.:1:360.

#tab_alpha = Vector{Float64}()
#tab_time = Vector{Float64}()
tab_theta = Vector{Float64}() 
tab_phi = Vector{Float64}()
tab_R = Vector{Float64}()
tab_R_12 = Vector{Float64}()
tab_modulation_depth = Vector{Float64}()
tab_nuclear_zq_frequency = Vector{Float64}()
tab_pair_contribution = Vector{Float64}()
tab_intensity = Vector{Float64}()


# scan of alpha
#for i in 1:length(alpha_list)[1]
    #alpha = deg2rad(alpha_list[i])
    #push!(tab_alpha,alpha)

# time-dependent simulation
#for t in 1:length(time_list)[1]
#    time = time_list[t]
#    println("time: ",time)
#    push!(tab_time,time)
 
df = DataFrame(theta=Float64[], phi=Float64[], R=Float64[],R_12=Float64[],modulation_depth=Float64[], nuclear_zq_frequency=Float64[], pair_contribution=Float64[], intensity=Float64[])

# scan of the magnetic field direction
for m in 1:length(theta)[1]
    for n in 1:length(phi)[1]
        println("Theta: ",theta[m])
        println("Phi: ",phi[n])
    
        rot_mat = rotate_solid(deg2rad(theta[m]),deg2rad(phi[n]))
		B0 = (rot_mat * B0_initial)

        #df = DataFrame(theta=Float64[], phi=Float64[], R=Float64[],R_12=Float64[],modulation_depth=Float64[], nuclear_zq_frequency=Float64[], pair_contribution=Float64[], intensity=Float64[])

        for i = 1:length(R)[1]
            for j = 1:length(R_12)[1]
                vec_r1,vec_r2,vec_r12 = get_distance_vectors(R[i],R_12[j],alpha)
                modulation_depth,nuclear_zq_frequency,pair_contribution,intensity = pss(vec_r1,vec_r2,vec_r12,g_e,g_n,time,B0)
            
                push!(tab_R,R[i])
                push!(tab_R_12,R_12[j])
                push!(tab_modulation_depth, modulation_depth)
                push!(tab_nuclear_zq_frequency, nuclear_zq_frequency)
                push!(tab_pair_contribution, pair_contribution)
                push!(tab_intensity, intensity)

                push!(df, (theta[m], phi[n], R[i], R_12[j], modulation_depth, nuclear_zq_frequency, pair_contribution, intensity))
            end
        end

        #theta_str = string(theta[m]) 
        #phi_str = string(phi[n])
        #file_name = "parameter_space_screening_theta_$(theta_str)_phi_$(phi_str)_R_$(R)_R12_$(R_12).csv" 
        #CSV.write(file_name, df)
        #CSV.write("parameter_space_screening_alpha_$alpha.csv", df)
    end
end

CSV.write("parameter_space_screening_magneticfieldscan_R_$(R)_R12_$(R_12).csv", df)