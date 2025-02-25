module rotation_B0

export rotation 

push!(LOAD_PATH,"../../src")

#using DetermineMagneticAxes
using LinearAlgebra

# determine the magentic axes from geometry
#system = System("../cudbm2.pdb","Pd",1)
#R_m = det_mag_axes(system)


function rotation(B0::Vector{Float64},theta::Float64,phi::Float64,R_m) 
     print("\n")
     println("Start Rotation of the magnetic field")
     println("==================================== \n")
     println("given B0: ", B0)
     println("Angle around the y axis: ", rad2deg(theta))
     println("Angle around the z axis: ", rad2deg(phi))

     # define rotation matrixes -- cartesian coordinate system 
     # rotation around y axis
     rot_mat_y = [cos(theta) 0 sin(theta); 0 1 0; +sin(theta) 0 cos(theta)]
     # rotation around z axis
     rot_mat_z = [cos(phi) +sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]

     # rotate B0
     B0_rot = transpose(R_m) *(rot_mat_z * (rot_mat_y * B0))
     println("Rotated magnetic field: ", B0_rot, "\n")

     return B0_rot
end


end



# rotation around a specific axis -- magnetic axes
#X = R_m[1,:]
#Y = R_m[2,:]
#Z = R_m[3,:]

# rotation around X
#rot_mat_X = [cos(psi)+X[1]^2*(1-cos(psi)) X[1]*X[2]*(1-cos(psi))-X[3]*sin(psi) X[1]*X[3]*(1-cos(psi))+X[2]*sin(psi);
#	     X[2]*X[1]*(1-cos(psi))+X[3]*sin(psi) cos(psi)+X[2]^2*(1-cos(psi)) X[2]*X[3]*(1-cos(psi))-X[1]*sin(psi);
#	     X[3]*X[1]*(1-cos(psi)).X[2]*sin(psi) X[3]'X[2]*(1-cos(psi))+X[1]*sin(psi) cos(psi)*X[3]^2*(1-cos(psi))] 

# rotation around Y
#rot_mat_Y = [cos(psi)+Y[1]^2*(1-cos(psi)) Y[1]*Y[2]*(1-cos(psi))-Y[3]*sin(psi) Y[1]*Y[3]*(1-cos(psi))+Y[2]*sin(psi);
#           Y[2]*X[1]*(1-cos(psi))+X[3]*sin(psi) cos(psi)+X[2]^2*(1-cos(psi)) X[2]*X[3]*(1-cos(psi))-X[1]*sin(psi);
#             X[3]*X[1]*(1-cos(psi)).X[2]*sin(psi) X[3]'X[2]*(1-cos(psi))+X[1]*sin(psi) cos(psi)*X[3]^2*(1-cos(psi))]

# rotate B0
#B0_rot = rot_mat_x * (rot_mat_y * (rot_mat_z * B0))
#println("new B0: ", B0_rot) 




