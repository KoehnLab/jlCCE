# run ParameterSpaceScreening

push!(LOAD_PATH,"../../src")

using LinearAlgebra
using Tables, CSV, DataFrames
using ParameterSpaceScreening 

# set parameters - R,R_12,alpha,ge,gn,time,B0
R = collect(range(1.,30.,100))
R_12 = collect(range(1.,30.,100))
alpha = pi/4

g_e = 2.002;
g_n = 5.586;
time = 10e-6
B0 = [0., 0., 1.];

# summarize results in tables
tab_R = Vector{Float64}()
tab_R_12 = Vector{Float64}()
tab_modulation_depth = Vector{Float64}()
tab_nuclear_zq_frequency = Vector{Float64}()
tab_pair_contribution = Vector{Float64}()
tab_intensity = Vector{Float64}()

df = DataFrame(R=Float64[],R_12=Float64[],modulation_depth=Float64[], nuclear_zq_frequency=Float64[], pair_contribution=Float64[], intensity=Float64[])

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

        push!(df, (R[i], R_12[j], modulation_depth, nuclear_zq_frequency, pair_contribution, intensity))
    end
end

CSV.write("parameter_space_screening.csv", df)
