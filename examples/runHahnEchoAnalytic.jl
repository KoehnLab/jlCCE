push!(LOAD_PATH,"../src")

using jlCCE
using Tables, CSV
using BenchmarkTools

# set system - use the simple constructor
spinsystem = SpinSystem("test_rotated","Pd",1)  # <--- note: put this file into the folder

# modify the values (SpinSystem creates a mutable object):
spinsystem.g_factor = [2.051,2.051,2.258]  # isotropic here

# defaults, no need to define here
spinsystem.magnetic_axes = [1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 1.0]
#spinsystem.magnetic_axes = [-0.0930991  -0.656122   0.748421; 0.423093   -0.708907  -0.567275; 0.901291    0.258756   0.343604]
spinsystem.B0 = [0., 0., 1.] # 10 kGaus = 1 Tesla
spinsystem.r_min = 0.
spinsystem.r_max = 35.
r_max = 10
spinsystem.r_max_bath = 10.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 50

# rotation parameters
spinsystem.theta = 90.
spinsystem.phi = 0. 

# run
#intensity = zeros(20,size(r_max)[1])
#for i in 1:size(r_max)[1]
    times,intensity = cce(spinsystem)
    #spinsystem.use_exp = false
    times,intensityNE = cce(spinsystem)
    spinsystem.simulation_type="highfield_analytic"
    #times,intensityEX,iCCE1,iCCE2 = cce(spinsystem)
    #print("\n")
    #print("Simulated intensity: ",intensity,"\n")
    #print("\n")
#end


#CSV.write("echo.csv", Tables.table([times intensity intensityNE intensityEX iCCE1 iCCE2]))



#print("\n")
#print("Our great result:\n")
#print("time of the Hahn echo: ",time_hahn_echo,"\n")
#print("Simulated intensity: ",intensity,"\n")
#print("\n")

