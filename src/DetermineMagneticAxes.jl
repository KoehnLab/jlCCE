# determine magnetic axes from the geometry 
# test
# Bsp mit Cudbm2 

using LinearAlgebra

Cu = ([2.979859,    1.445000,    5.902247]);
O1 = ([3.726962,    2.999820,    6.841884]);
O2 = ([3.992780,    1.791800,    4.268505]);
O3 = ([2.232756,   -0.109820,    4.962609]);
O4 = ([1.966938,    1.098200,    7.535989]);


# determine z direction of the magnetic axes
#         z
#         |
#    O    |   O
#        Cu ------ x 
#    O    |   O
#         y 

# distance between Cu and O1, O4 --> spannen die Molekülebene auf
r_cuo1 = Cu - O1
r_cuo4 = Cu - O4 

# z ist orthogonal zu den beiden Vektoren
z =  (cross(r_cuo1,r_cuo4))/norm(cross(r_cuo1,r_cuo4))

# middle point between O1 & O4, O3 & O4
M_o1o4 = (O1 + O4)/2
M_o3o4 = (O3 + O4)/2

# determnine x and y
x = (M_o1o4 - Cu)/norm(M_o1o4 - Cu)
y = (M_o3o4 - Cu)/norm(M_o3o4 - Cu)

println(x)
println(y)
println(z)

# magnetix axes matrix
R_m = [x y z]
println(R_m)

# überprüfen, ob die Achsen korrekt gesetzt wurden
det_R_m = det(R_m)
println(det_R_m)



