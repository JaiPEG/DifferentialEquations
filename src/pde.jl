include("function.jl")

"""
Neumann boundary conditions on a 1-dimensional domain.
"""
struct Neumann{T}
	a::T
	b::T
end

"""
Dirichlet boundary conditions on a 1-dimensional domain.
"""
struct Dirichlet{T}
	a::T
	b::T
end

"""
Impose Neumann boundary conditions. Modify only the endpoints.
"""
function neumann!(bc::Neumann{R}, f::Hat{R, a, b, n}) where {R, a, b, n}
	if n >= 3
		h = (b - a) / (n - 1)
		f.coeffs[1] = f.coeffs[2] - bc.a * h
		f.coeffs[n] = f.coeffs[n - 1] + bc.a * h
	end
end

"""
Impose Dirichlet boundary conditions. Modify only the endpoints.
"""
function dirichlet!(bc::Dirichlet{R}, f::Hat{R, a, b, n}) where {R, a, b, n}
	if n >= 2
		f.coeffs[1] = bc.a
		f.coeffs[n] = bc.b
	end
end

"""
Evolve a function by one time-step under the diffusion equation with Neumann
boundary conditions using the given autonomous ode solver.

ode :: ((Hat{R, a, b, n} -> Hat{R, a, b, n}), Hat{R, a, b, n}, R) -> R
	Autonomous ODE solver
	E.g.: rk2step
bc :: Neumann
	Neumann boundary conditions
f :: Hat{R, a, b, n}
	Solution at a fixed time t
h :: R
	Time-step

-> :: Hat{R, a, b, n}
	Approximation to the solution at fixed time t + h

See also: diffusionD
"""
function diffusionN(ode, bc::Neumann, f::Hat{R, a, b, n}, h::R)::Hat{R, a, b, n} where {R, a, b, n}
	# state vector is just f since only the first time-derivative is involved
	g = ode(deriv2, f, h)
	# Note: deriv2 keeps boundary slopes constant
	# So we would get automatic Neumann BC for those initial slopes
	neumann!(bc, g)
	g
end

"""
Evolve a function by one time-step under the diffusion equation with Dirichlet
boundary conditions using the given autonomous ode solver.

ode :: ((Hat{R, a, b, n} -> Hat{R, a, b, n}), Hat{R, a, b, n}, R) -> R
	Autonomous ODE solver
	E.g.: rk2
bc :: Dirichlet
	Dirichlet boundary conditions
f :: Hat{R, a, b, n}
	Solution at a fixed time t
h :: R
	Time-step

-> :: Hat{R, a, b, n}
	Approximation to the solution at fixed time t + h

See also: diffusionN
"""
function diffusionD(ode, bc::Dirichlet, f::Hat{R, a, b, n}, h::R)::Hat{R, a, b, n} where {R, a, b, n}
	# state vector is just f since only the first time-derivative is involved
	g = ode(deriv2, f, h)
	dirichlet!(bc, g)
	g
end

"""
Evolve a state vector by one time-step under the wave equation using the
given autonomous ode solver.

ode :: ((V -> V), V, R) -> R
	Autonomous ODE solver
	E.g.: rk2
s :: Vector{V}
	State vector.
	First element is the time derivative of the solution at time t.
	Second element is the space derivative of the solution at time t.
h :: R
	Time-step
-> :: Vector{V}
	State vector at time t + h
"""
function wave(ode, s::Vector{V}, h::R)::Vector{V} where {R, V}
	rhs(s) = [deriv(s[2]), deriv(s[1])]
	ode(rhs, s, h)
end

"""
Evolve a state vector by one time-step under the wave equation using the
given autonomous ode solver and apply fixed Dirichlet boundary conditions.
Specifically, set the endpoints of the time deriviative of the solution
to zero.

ode :: ((Hat{R, a, b, n} -> Hat{R, a, b, n}), Hat{R, a, b, n}, R) -> R
	Autonomous ODE solver
	E.g.: rk2
s :: Vector{Hat{R, a, b, n}}
	State vector.
	First element is the time derivative of the solution at time t.
	Second element is the space derivative of the solution at time t.
h :: R
	Time-step
-> :: Vector{Hat{R, a, b, n}}
	State vector at time t + h
"""
function waveD0(ode, s::Vector{Hat{R, a, b, n}}, h::R)::Vector{Hat{R, a, b, n}} where {R, a, b, n}
	newS = wave(ode, s, h)
	# Set time derivative of solution to zero at the endpoints.
	# Essentially apply a fixed Dirichlet boundary condition on the
	# original solution.
	# Idea by Benjamin de Bruyne.
	newS[1].coeffs[1] = 0.0
	newS[1].coeffs[end] = 0.0
	newS
end
