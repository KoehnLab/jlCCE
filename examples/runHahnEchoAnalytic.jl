push!(LOAD_PATH,"../src")

using jlCCE
using Tables, CSV
using BenchmarkTools

# set system - use the simple constructor
spinsystem = SpinSystem("cudbm2.pdb","Pd",1)  # <--- note: put this file into the folder
#spinsystem = SpinSystem("pddbm2.cif","Pd",1)  # <--- note: put this file into the folder

# modify the values (SpinSystem creates a mutable object):
spinsystem.g_factor = [2.051,2.051,2.258]  # isotropic here
# pddbm2: g = 2.25
# vmnt3: 1.9846,1.9846,1.9846



# defaults, no need to define here
#voacac2.magnetic_axes = [1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 1.0]
spinsystem.magnetic_axes = [-0.0930991  -0.656122   0.748421; 0.423093   -0.708907  -0.567275; 0.901291    0.258756   0.343604]
#voacac2.nuc_spin_bath = "H"
#voacac2.gn_spin_bath = 5.58569468
spinsystem.B0 = [2.9597, -2.4458, 1.3158]*10000 # 10 kGaus = 1 Tesla
spinsystem.r_min = 0.
#r_max = range(4.5,35.5,length=150) # AA
#r_max = range(4.5,5.5,length=150) # AA
r_max = 10
spinsystem.r_max_bath = 10.
spinsystem.t_min = 0.
spinsystem.t_max = 15e-6 # s
spinsystem.n_time_step = 20


# run
#intensity = zeros(20,size(r_max)[1])
#for i in 1:size(r_max)[1]
    spinsystem.r_max = 10.0 #r_max[i]
    print("r_max: ",r_max,"\n")
    times,intensity = cce(spinsystem)
    spinsystem.use_exp = false
    times,intensityNE = cce(spinsystem)
    spinsystem.simulation_type="exact"
    times,intensityEX,iCCE1,iCCE2 = cce(spinsystem)
    #print("\n")
    #print("Simulated intensity: ",intensity,"\n")
    #print("\n")
#end


CSV.write("echo.csv", Tables.table([times intensity intensityNE intensityEX iCCE1 iCCE2]))



#print("\n")
#print("Our great result:\n")
#print("time of the Hahn echo: ",time_hahn_echo,"\n")
#print("Simulated intensity: ",intensity,"\n")
#print("\n")

