# Equations and LinearOperators for OceanTurb.jl

export
  LinearOperator,
  DifferentialOperator
  Equation,
  StandardEquation,
  calc_rhs!(rhs, model)

abstract type LinearOperator end

struct DifferentialOperator{A} <: LinearOperator 
  L::A
end

struct DiagonalOperator{A} <: LinearOperator 
  L::A
end

struct NullOperator <: LinearOperator end

abstract type Equation end

"""
    Equation(L, N)

Defines the PDE for a solution vector Φ such that

∂t(Φ) = L*Φ + N(Φ).
      = rhs
"""
struct StandardEquation{LO} <: Equation
  L::LO # Linear operator
  N::Function # Nonlinear operator
end

function calc_rhs!(rhs, model)
  #equation.L(rhs, model)
  equation.N(rhs, model)
  nothing
end
