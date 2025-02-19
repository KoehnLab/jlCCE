push!(LOAD_PATH,"../../src")

using DetermineMagneticAxes
using LinearAlgebra

# set system - use the simple constructor
spinsystem = SpinSystem("../cudbm2.pdb","Pd",1) 

# run
R_m = det_mag_axes(spinsystem)










# example: cudbm2

# coordinates of Cu and O
#Cu = ([2.979859,    1.445000,    5.902247]);
#O1 = ([3.726962,    2.999820,    6.841884]);
#O2 = ([3.992780,    1.791800,    4.268505]);
#O3 = ([2.232756,   -0.109820,    4.962609]);
#O4 = ([1.966938,    1.098200,    7.535989]);

# determine x,y,z contribution of the magnetic axes

#         z
#         |
#    O2   |   O1
#        Cu ------ x
#    O3  /    O4
#       /
#      y

# middle point between O1 & O4, O3 & O4
#M_o1o4 = (O1 + O4)/2
#M_o3o4 = (O3 + O4)/2

# determnine x and y direction of the magnetic axes
#x = (M_o1o4 - Cu)/norm(M_o1o4 - Cu)
#y = (M_o3o4 - Cu)/norm(M_o3o4 - Cu)

# determine z direction of the magnetic axes
#z =  (cross(x,y))#/norm(cross(x,y))
#norm_z = norm(z)
#println("Norm of z: ",norm_z)

# matrix of the magnetic axes
#R_m = [x y z]

#println("Magnetix axes: ")
#@printf " x = [%2.4f %2.4f  %2.4f]\n" x[1] x[2] x[3]
#@printf " y = [%2.4f %2.4f  %2.4f]\n" y[1] y[2] y[3]
#@printf " z = [%2.4f %2.4f  %2.4f]\n\n" z[1] z[2] z[3]

# überprüfen, ob die Achsen korrekt gesetzt wurden
#det_R_m = det(R_m)
#println(det_R_m)

