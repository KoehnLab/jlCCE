push!(LOAD_PATH,"../../src")

using QuadGK
using Lebedev
using Trapz
using LinearAlgebra
using Tables, CSV, DataFrames
using Statistics
using RotationMatrices

# set parameters
r = collect(range(1.,35.,150))
r12 = collect(range(1.,35.,150))
time_hahn_echo = collect(range(0,15e-6,50))
B0 = [0.,0.,1.]


"""
    get_mod_depth(alpha,R,R12,B0)

determine the modulation depth depending on the distances R,R12 scanning alpha

input: R,R12 - Jacobi coordinates of the distances between the spins
       alpha - angle between the distance vectors
       B0 - magnetic field 

returns lambda (modulation depth)

"""
function get_mod_depth(alpha,R,R12,B0)
    #r1_norm = sqrt( R^2 - R*R12*cos(alpha) + 0.5*R12*cos(alpha)^2 - 0.25*R12^2*sin(alpha)^2 ) 
    #r1_norm = sqrt( R^2 + R*R12*cos(alpha) + 0.5*R12*cos(alpha)^2 + 0.25*R12^2*sin(alpha)^2 )
    g_e = 2.002;
    g_n = 5.586; 
    mu_b = 9.274010066e-21
    mu_n = 5.050783739e-24
    hbar = 1.054571817e-27

    gamma_s = (g_e * mu_b)/ hbar
    gamma_n = (g_n * mu_n) / hbar

    vec_r = [R, 0., 0.]
    vec_r12 = [R12 * cos(alpha), R12 * sin(alpha), 0.]
    vec_r1 = vec_r + 1/2 * vec_r12
    vec_r2 = vec_r - 1/2 * vec_r12

    aacm = 1e-8
    r1 = vec_r1 * aacm
    r2 = vec_r2 * aacm
    r12 = vec_r12 * aacm

    r1_cross_B0 = cross(r1,B0)
    r1_dot_B0 = dot(r1,B0)
    theta1 = atan(norm(r1_cross_B0),r1_dot_B0)
    r1_norm = max(norm(r1),1e-8)

    r2_cross_B0 = cross(r2,B0)
    r2_dot_B0 = dot(r2,B0)
    theta2 = atan(norm(r2_cross_B0),r2_dot_B0)
    r2_norm = max(norm(r2),1e-8)

    A1 = - gamma_s*gamma_n*hbar * ((1-3*cos(theta1))/r1_norm^3)
    A2 = - gamma_s*gamma_n*hbar * ((1-3*cos(theta2))/r2_norm^3)

    r12_cross_B0 = cross(r12,B0);
    r12_dot_B0 = dot(r12,B0);
    theta_12 = atan(norm(r12_cross_B0),r12_dot_B0);
    r12_norm = norm(r12)
    b12 = - gamma_n^2 * hbar * ((1 - 3* cos(theta_12)^2)/r12_norm^3); 

    c12 = (A1 - A2)/(4*b12)

    lambda = 4 * (c12^2/((1+c12^2)^2))
    #println("modulation depth: ", lambda)
    return lambda
end


"""
    get_nzq_frequency(alpha,R,R12,B0)

determine the nuclear zero-quantum frequency depending on the distances R,R12 scanning alpha

input: R,R12 - Jacobi coordinates of the distances between the spins
       alpha - angle between the distance vectors
       B0 - magnetic field 

returns w12 (nuclear zero-quantum frequency)

"""
function get_nzq_frequency(alpha,R,R12,B0)
    g_e = 2.002;
    g_n = 5.586; 
    mu_b = 9.274010066e-21
    mu_n = 5.050783739e-24
    hbar = 1.054571817e-27

    gamma_s = (g_e * mu_b)/ hbar
    gamma_n = (g_n * mu_n) / hbar

    vec_r = [R, 0., 0.]
    vec_r12 = [R12 * cos(alpha), R12 * sin(alpha), 0.]
    vec_r1 = vec_r + 1/2 * vec_r12
    vec_r2 = vec_r - 1/2 * vec_r12

    aacm = 1e-8
    r1 = vec_r1 * aacm
    r2 = vec_r2 * aacm
    r12 = vec_r12 * aacm

    r1_cross_B0 = cross(r1,B0)
    r1_dot_B0 = dot(r1,B0)
    theta1 = atan(norm(r1_cross_B0),r1_dot_B0)
    r1_norm = max(norm(r1),1e-8)

    r2_cross_B0 = cross(r2,B0)
    r2_dot_B0 = dot(r2,B0)
    theta2 = atan(norm(r2_cross_B0),r2_dot_B0)
    r2_norm = max(norm(r2),1e-8)

    A1 = - gamma_s*gamma_n*hbar * ((1-3*cos(theta1))/r1_norm^3)
    A2 = - gamma_s*gamma_n*hbar * ((1-3*cos(theta2))/r2_norm^3)

    r12_cross_B0 = cross(r12,B0);
    r12_dot_B0 = dot(r12,B0);
    theta_12 = atan(norm(r12_cross_B0),r12_dot_B0);
    r12_norm = norm(r12)
    b12 = (- 1/4) * gamma_n^2 * hbar * ((1 - 3* cos(theta_12)^2)/r12_norm^3); 

    c12 = (A1 - A2)/(4*b12)

    w12 = 2. * b12 * sqrt(1 + c12^2)
    #println("nuclear zero-quantum frequency: ", w12)
    return w12 
end


"""
    get_pair_contribution(alpha,R,R12,B0)

determine the pair contribution using the averaged modulation depth and the averaged nuclear zero-quantum frequency

input: averaged_modulation_depth - averaged modulation depth
       averaged_nzq_frequency - averaged nuclear zero-quantum frequency
       time_hahn_echo - time set for which the pair contribution is to be computed

returns pair_contribution 

"""
function get_pair_contribution(averaged_modulation_depth,averaged_nzq_frequency,time_hahn_echo)
    pair_contribution = ones(size(time_hahn_echo)[1])
    for t in 1:size(time_hahn_echo)[1]
        pair_contribution[t] = 1 - averaged_modulation_depth * sin(1/4 * averaged_nzq_frequency * time_hahn_echo[t])^4
    end    

    return pair_contribution
end



df = DataFrame(r=Float64[],r12=Float64[], averaged_modulation_depth=Float64[], averaged_nzq_frequency=Float64[], time_hahn_echo=String[], pair_contribution=String[])

# run Gauss-Kronrod integration depending on r,r12 - averaging over alpha 
for i = 1:length(r)[1]
    for j = 1:length(r12)[1]
        R = r[i]
        R12 = r12[j] 
        #println("distances: ",R," ",R12)
        
        int, err = quadgk(alpha -> get_mod_depth(alpha,R,R12,B0), 0, 2*pi, rtol=1e-8)
       #println("integration value: ",int)
       #println("error: ",err)

        averaged_modulation_depth = int/(2*pi)
       #println("averaging modulation depth value: ",averaged_modulation_depth)

        int, err = quadgk(alpha -> get_nzq_frequency(alpha,R,R12,B0), 0, 2*pi, rtol=1e-8)
        #println("integration value: ",int)
        #println("error: ",err)

        averaged_nzq_frequency = int/(2*pi)
        #println("averaging modulation depth value: ",averaged_modulation_depth)

        averaged_pair_contribution = get_pair_contribution(averaged_modulation_depth,averaged_nzq_frequency,time_hahn_echo)

        push!(df, (R, R12, averaged_modulation_depth, averaged_nzq_frequency, join(time_hahn_echo,","), join(averaged_pair_contribution, ",")))

    end
end

CSV.write("gauss_kronrod_integration_alpha.csv", df)
