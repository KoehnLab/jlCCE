module SpinBase

using LinearAlgebra

export unit_mat, sz_mat, sp_mat, sm_mat, sy_mat, sx_mat

"""
    unit_mat(dim)

returns a unit matrix of ComplexF64 type

"""
function unit_mat(dim::Integer)
    u = (1. + 0. * im) * I(dim)
    return u
end

"""
    sz_mat(S)

returns the Z component of a spin operator in the basis
|M=S> |M=S-1> ... |M=-S>
the result is of ComplexF64 type

See also [`sx_mat`](@ref), [`sy_mat`](@ref), [`sp_mat`](@ref), [`sm_mat`](@ref).

"""
function sz_mat(S::Float64)
    dim = convert(Integer,(2. * S + 1.))
    smat = zeros(ComplexF64,dim,dim)
    for ii = 1:dim
        smat[ii,ii] = S - ii * 1.0 + 1.0
    end
    return smat
end

"""
    sp_mat(S)

returns the matrix elements of the raising operator in the basis
|M=S> |M=S-1> ... |M=-S>
the result is of ComplexF64 type

See also [`sm_mat`](@ref), [`sx_mat`](@ref), [`sy_mat`](@ref), [`sz_mat`](@ref).

"""
function sp_mat(S::Float64)
    dim = convert(Integer,(2. * S + 1.))
    smat = zeros(ComplexF64,dim,dim)
    SSp1 = S * (S + 1.)
    M = S - 1.
    for ii = 1:dim-1
        smat[ii,ii+1] = sqrt(SSp1 - M * (M + 1.))
        M = M - 1.
    end
    return smat
end

"""                               
    sm_mat(S)                                         
                                                          
returns the matrix elements of the lowering operator in the basis
|M=S> |M=S-1> ... |M=-S>
the result is of ComplexF64 type

See also [`sp_mat`](@ref), [`sx_mat`](@ref), [`sy_mat`](@ref), [`sz_mat`](@ref).

"""
function sm_mat(S::Float64)
    dim = convert(Integer,(2. * S + 1.))
    smat = zeros(ComplexF64,dim,dim)
    SSp1 = S * (S + 1.)
    M = S
    for ii = 2:dim
        smat[ii,ii-1] = sqrt(SSp1 - M * (M - 1.))
        M = M - 1.
    end
    return smat
end

"""
    sx_mat(S)

returns the matrix elements of the S_x operator in the basis
|M=S> |M=S-1> ... |M=-S>
the result is of ComplexF64 type

See also [`sy_mat`](@ref), [`sz_mat`](@ref), [`sp_mat`](@ref), [`sm_mat`](@ref).

"""
function sx_mat(S::Float64)
    return 0.5*(sp_mat(S) + sm_mat(S))
end

"""
    sy_mat(S)

returns the matrix elements of the S_y operator in the basis
|M=S> |M=S-1> ... |M=-S>
the result is of ComplexF64 type

See also [`sx_mat`](@ref), [`sz_mat`](@ref), [`sp_mat`](@ref), [`sm_mat`](@ref).

"""
function sy_mat(S::Float64)
    return -0.5im * (sp_mat(S) - sm_mat(S))
end


end
