using QuadGK
using Trapz
using LinearAlgebra
using Tables, CSV, DataFrames
using Statistics

# Integration hier ausprobieren! mit den package von Andreas

# read out csv file - parameter space screening: alpha scan 
df = CSV.read("C:\\Users\\sucha\\iCloudDrive\\Promotion\\jlCCE\\examples\\ParameterSpaceScreening\\parameter_space_screening_alphascan_R_15.0_R12_2.0.csv", DataFrame)

# convert string of the intensity values to vector of the type Float64
#df.intensity = [parse.(Float64, split(x, ",")) for x in df.intensity]

# averaging intensity and coherence time
#intensity_average = mean(df.intensity)

# set R, R12, alpha 
R12 = 10.
R = 15.
println("distances: ", R, " ", R12)

alpha = df.alpha
#println("alpha: ", alpha)

# function of the modulation depth depending on alpha
function get_mod_depth(alpha,R,R12)
    #r1_norm = sqrt( R^2 - R*R12*cos(alpha) + 0.5*R12*cos(alpha)^2 - 0.25*R12^2*sin(alpha)^2 ) 
    #r1_norm = sqrt( R^2 + R*R12*cos(alpha) + 0.5*R12*cos(alpha)^2 + 0.25*R12^2*sin(alpha)^2 )
    g_e = 2.002;
    g_n = 5.586; 
    mu_b = 9.274010066e-21
    mu_n = 5.050783739e-24
    hbar = 1.054571817e-27
    B0 = [0.,0.,1.]

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
    r1_norm = norm(r1)

    r2_cross_B0 = cross(r2,B0)
    r2_dot_B0 = dot(r2,B0)
    theta2 = atan(norm(r2_cross_B0),r2_dot_B0)
    r2_norm = norm(r2)

    A1 = - gamma_s*gamma_n*hbar * ((1-3*cos(theta1))/r1_norm^3)
    A2 = - gamma_s*gamma_n*hbar * ((1-3*cos(theta2))/r2_norm^3)

    r12_cross_B0 = cross(r12,B0);
    r12_dot_B0 = dot(r12,B0);
    theta_12 = atan(norm(r12_cross_B0),r12_dot_B0);
    r12_norm = norm(r12)
    b12 = (- 1/4) * gamma_n^2 * hbar * ((1 - 3* cos(theta_12)^2)/r12_norm^3); 

    c12 = (A1 - A2)/(4*b12)

    lambda = c12^2/((1+c12^2)^2) 
    #println("modulation depth: ", lambda)
    return lambda
end

#lambda = get_mod_depth.(alpha,R,R12)
#print("lamdba: ",lambda)

# one-dimensional numerical integration ("quadrature") using adaptive Gauss–Kronrod quadrature 
int, err = quadgk(alpha -> get_mod_depth(alpha,R,R12), 0, 2*pi, rtol=8)
println("integration value: ",int)
println("error: ",err)

alpha_average = int/(2*pi)
println("averaging alpha value: ",alpha_average, " rad, ", rad2deg(alpha_average), " deg")

# trapezoidal intgration 
#integral_trapez = trapz(alpha, df.modulation_depth)  # Näherungsweise Integration
#println("Näherungsintegral mit Trapezregel: ", integral_trapez)

#alpha_average = integral_trapez/(2*pi)
#println("averaging alpha value: ",alpha_average)