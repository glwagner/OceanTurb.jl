# Equations and LinearOperators for OceanTurb.jl

export
  LinearOperator,
  DifferentialOperator,
  Equation,
  StandardEquation,
  calc_rhs!

abstract type LinearOperator end

struct DifferentialOperator{A} <: LinearOperator 
  L::A
end

struct DiagonalOperator{A} <: LinearOperator 
  L::A
end

abstract type Equation end

"""
    StandardEquation(L, N)

Defines the PDE for a solution vector Φ such that

∂t(Φ) = L*Φ + N(Φ) ,
      = rhs
"""
struct StandardEquation{LO,NO} <: Equation
  L::LO # Linear operator
  N::NO # Nonlinear operator
end

function calc_rhs!(rhs, model)
  #equation.L(rhs, model)
  model.equation.N(rhs, model)
  nothing
end

#
# Zero-flux boundary conditions
#

function zero_flux!(f::FaceField)
  f.data[1] = 0
  f.data[f.grid.nz+1] = 0
  nothing
end

function zero_flux!(c::CellField)
  c.data[0] = c.data[1]
  c.data[end] = c.data[end-1]
  nothing
end
