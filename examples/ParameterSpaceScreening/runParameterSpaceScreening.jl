# run ParameterSpaceScreening

using LinearAlgebra
using ParameterSpaceScreening 

# set parameters - R,R_12,alpha,ge,gn,time,B0
R = collect(range(1.,30.,50))
R_12 = collect(range(1.,30.,50))
alpha = pi/4

g_e = 2.002;
g_n = 5.586;
time = 10e-6
B0 = [1., 0., 0.];

# summarize results in tables
tab_R = Vector{Float64}()
tab_R_12 = Vector{Float64}()
tab_modulation_depth = Vector{Float64}()
tab_nuclear_zq_frequency = Vector{Float64}()
tab_pair_contribution = Vector{Float64}()
tab_intensity = Vector{Float64}()

df = DataFrame(modulation_depth=Float64[], nuclear_zq_frequency=Float64[], pair_contribution=Float64[], intensity=Float64[])

for i = 1:length(R)
    for j = 1:length(R_12)
        vec_r1,vec_r2,vec_r12 = get_distance_vectors(R[i],R_12[j],)
        modulation_depth,nuclear_zq_frequency,pair_contribution,intensity = pss(vec_r1,vec_r2,vec_r12,ge,gn,time,B0)
        push!(tab_R,R[i])
        push!(tab_R_12,R_12[j])
        push!(tab_modulation_depth, modulation_depth)
        push!(tab_nuclear_zq_frequency, nuclear_zq_frequency)
        push!(tab_pair_contribution, pair_contribution)
        push!(tab_intensity, intensity)
    end
end

CSV.write("parameter_space_screening.csv", df)