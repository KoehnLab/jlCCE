module RotationMatrices

export rotate_x,rotate_y,rotate_z,rotate_solid 

push!(LOAD_PATH,"../../src")

using LinearAlgebra

"""
    rotate_x(angle::Float64)

returns rotation matrix for rotation (counter clockwise) around x axis by angle (mathematically positive) 

"""
function rotate_x(angle::Float64)
     rot_mat = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)]
     return rot_mat
end


"""
    rotate_y(angle::Float64)

returns rotation matrix for rotation (counter clockwise) around y axis by angle (mathematically positive) 

"""
function rotate_y(angle::Float64)
     rot_mat = [cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle)]
     return rot_mat
end


"""
    rotate_z(angle::Float64)

returns rotation matrix for rotation (counter clockwise) around z axis by angle (mathematically positive)

"""
function rotate_z(angle::Float64)
     rot_mat = [cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1]
     return rot_mat
end

"""
    rotate_solid(theta::Float64,phi::Float64)

returns rotation matrix for rotation (counter clockwise) around y axis by theta and z axis by phi

"""
function rotate_solid(theta::Float64,phi::Float64)
     rot_mat = rotate_z(phi) * rotate_y(theta)
     return rot_mat
end

end
