# run ParameterSpaceScreening

push!(LOAD_PATH,"../../src")

using LinearAlgebra
using Tables, CSV, DataFrames
using Printf
using ParameterSpaceScreening 
using RotationMatrices

# set parameters - R,R_12,alpha,ge,gn,time,B0
R = collect(range(1.,35.,150))
R_12 = collect(range(1.,35.,150))
alpha = deg2rad(45) 
g_e = 2.002;
g_n = 5.586;
time = 10e-6
B0 = [0., 0., 1.];
 
df = DataFrame(R=Float64[],R_12=Float64[],modulation_depth=Float64[], nuclear_zq_frequency=Float64[], pair_contribution=Float64[], intensity=Float64[])

for i = 1:length(R)[1]
    for j = 1:length(R_12)[1]
        vec_r1,vec_r2,vec_r12 = get_distance_vectors(R[i],R_12[j],alpha)
        modulation_depth,nuclear_zq_frequency,pair_contribution,intensity = pss(vec_r1,vec_r2,vec_r12,g_e,g_n,time,B0)
            
        push!(df, (R[i], R_12[j], modulation_depth, nuclear_zq_frequency, pair_contribution, intensity))
    end
end

alpha = rad2deg(alpha)
CSV.write("parameter_space_screening_alpha_$(alpha)_degree.csv", df)


