# run ParameterSpaceScreening: scan of alpha for fixed distances r12,R

push!(LOAD_PATH,"../../src")

using LinearAlgebra
using Tables, CSV, DataFrames
using Printf
using ParameterSpaceScreening 
using RotationMatrices

# set parameters - R,R_12,alpha,ge,gn,time,B0
R = 20.
R_12 = 15.
alpha_list = 0.:2:360.
 
g_e = 2.002;
g_n = 5.586;
time = 10e-6
B0 = [0., 0., 1.];

df = DataFrame(R=Float64[],R_12=Float64[], alpha=Float64[], modulation_depth=Float64[], nuclear_zq_frequency=Float64[], pair_contribution=Float64[], intensity=Float64[])

# scan of alpha
for k in 1:length(alpha_list)[1]
    alpha = deg2rad(alpha_list[k])

    vec_r1,vec_r2,vec_r12 = get_distance_vectors(R,R_12,alpha)
    modulation_depth,nuclear_zq_frequency,pair_contribution,intensity = pss(vec_r1,vec_r2,vec_r12,g_e,g_n,time,B0)
            
    push!(df, (R, R_12, alpha_list[k], modulation_depth, nuclear_zq_frequency, pair_contribution, intensity))
end

CSV.write("parameter_space_screening_alphascan_R_$(R)_R12_$(R_12).csv", df)