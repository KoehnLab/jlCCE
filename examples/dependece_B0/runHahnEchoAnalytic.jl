push!(LOAD_PATH,"../../src")

using jlCCE
using Tables, CSV
using BenchmarkTools

# set system - use the simple constructor
spinsystem = SpinSystem("../cudbm2.pdb","Pd",1)  # <--- note: put this file into the folder

# modify the values (SpinSystem creates a mutable object):
spinsystem.s_el = 0.5
spinsystem.g_factor = [2.051,2.051,2.258] 
spinsystem.magnetic_axes = [-0.0930991  -0.656122   0.748421; 0.423093   -0.708907  -0.567275; 0.901291    0.258756   0.343604]
spinsystem.B0 = [0.7484, -0.5673, 0.3436]
spinsystem.r_min = 0.
r_max = 35
spinsystem.r_max_bath = 10.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 50

# run
    spinsystem.r_max = 35.0 
    print("r_max: ",r_max,"\n")   
    times,intensity = cce(spinsystem)
    spinsystem.use_exp = false
    #times,intensityNE = cce(spinsystem)
    spinsystem.simulation_type="highfield_analytic"
    #times,intensityEX,iCCE1,iCCE2 = cce(spinsystem)
    #print("\n")
    print("Simulated intensity: ",intensity,"\n")
    #print("\n")


#CSV.write("echo.csv", Tables.table([times intensity intensityNE intensityEX iCCE1 iCCE2]))



#print("\n")
#print("Our great result:\n")
#print("time of the Hahn echo: ",time_hahn_echo,"\n")
#print("Simulated intensity: ",intensity,"\n")
#print("\n")

