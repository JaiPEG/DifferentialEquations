include("util.jl")

# Define data structures and algorithms for representing and manipulating
# functions.

#
# Definitions
#

# Define F(K, V, a, b, n) to be the vector space of functions from n uniformly
# sampled points between a and b to some K-vector space V.

# Define the hat basis to be the basis of F(R, R, a, b, n), where R is some
# approximation to the field of real numbers (e.g.: Float64, Rational{Int64}),
# consisting of continuous piecewise linear functions where the ith such
# function is 1 at the ith sample point and 0 at the other sample points.

"""
Represent a function in the vector space spanned by the hat basis functions.
Store just the coefficients of the representation in the hat basis.

R::DataType, an ordered field
a::R
b::R
n::Int
"""
struct Hat{R, a, b, n}
	coeffs::Vector{R}
end

function Base.:(==)(f::Hat{R, a, b, n}, g::Hat{R, a, b, n})::Bool where {R, a, b, n}
	f.coeffs == g.coeffs
end

#
# Vector space
#

"""
Return the identity of the vector space Hat{R, a, b, n}.
"""
function zeroHat(a::R, b::R, n::Int)::Hat{R, a, b, n} where {R}
	Hat{R, a, b, n}(zeros(n))
end

"""
Return the identity of the vector space Hat{R, a, b, n} using the
interval size and number of samples from the given element.
"""
function zeroHat(f::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(zeros(n))
end

function Base.:+(f::Hat{R, a, b, n}, g::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(f.coeffs + g.coeffs)
end

function Base.:-(f::Hat{R, a, b, n}, g::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(f.coeffs - g.coeffs)
end

function Base.:-(f::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(-f.coeffs)
end

function Base.:*(c::R, f::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(c * f.coeffs)
end

#
# Projection
#

"""
Sample a function [a, b] -> R at n points uniformly spaced between a and b.
"""
function sampleUniform(a::R, b::R, n::Int, f)::Vector{R} where {R}
	h = (b - a) / (n - 1)
	[f(a + h * i) for i in 0:(n - 1)]
end

"""
Project a function [a, b] -> R into the Hat{R} subspace. This equates to
sampling the function at n points uniformly spaced between a and b.
"""
function projectHat(a::R, b::R, n::Int, f)::Hat{R, a, b, n} where {R}
	Hat{R, a, b, n}(sampleUniform(a, b, n, f))
end

#
# Evaluation
#

"""
Evaluate the ith hat basis function at a point x between a and b where i goes
from 0 to N - 1.

Currently throws an exception if x is not in the interval [a, b].

Really, this shouldn't be a publically visible function.
"""
function hat(a::R, b::R, n::Int, i::Int, x::R)::R where {R}
	if !(0 <= i < n)
		throw(DomainError(i, "Hat basis index out of bounds."))
	# TODO: consider removing
	elseif !(a <= x <= b)
		throw(DomainError(x, "Hat basis function evaluation out of bounds."))
	else
		h = (b - a) / R(n - 1)
		xip = a + h * (i - 1)
		xi  = a + h * i
		xin = a + h * (i + 1)
		if x <= xip
			return R(0)
		elseif x <= xi
			return (x - xip) / h
		elseif x <= xin
			return -(x - xin) / h
		else
			return R(0)
		end
	end
end

"""
Evaluate a function in the hat basis representation at a point x between a and
b.
"""
function eval(f::Hat{R, a, b, n}, x::R)::R where {R, a, b, n}
	sum(f.coeffs[i + 1] * hat(a, b, n, i, x) for i in 0:(n - 1))
end

#
# Quadrature
#

"""
Compute the quadrature (numerical integration) for a function in the hat basis
representation over its whole domain.
"""
function quad(f::Hat{R, a, b, n})::R where {R, a, b, n}
	h = (b - a) / (n - 1)
	s = R(0)
	for i in 0:(n - 1)
		# weight of each basis function
		# boundary basis functions have half contribution
		w = (i == 0 || i == n - 1 ? R(1//2) : R(1))
		s += w * h * f.coeffs[i + 1]
	end
	s
end

"""
Compute the quadrature (numerical integration) for a function in the hat basis
representation over the given interval.
"""
function quad(f::Hat{R, a, b, n}, x1::R, x2::R)::R where {R, a, b, n}
	# Bound too far left
	if x1 < a
		throw(DomainError(x1, "Quadrature bound out of bounds."))
	# Bound too far right
	elseif x2 > b
		throw(DomainError(x2, "Quadrature bound out of bounds."))
	# Bounds switched
	elseif x1 > x2
		-quad(f, x2, x1)
	else
		# println("Bounds: ", x1, " ", x2)
		h = (b - a) / (n - 1)

		function trap(height, base1, base2)
			0.5*height*(base1 + base2)
		end

		# Map into sample range
		y1 = mapRange(a, b, 1.0, R(n), x1)
		y2 = mapRange(a, b, 1.0, R(n), x2)
		# println("Sample r: ", y1, " ", y2, " ")

		# Get inner sample bounds
		s1 = Int(round(y1, RoundUp))
		s2 = Int(round(y2, RoundDown))
		# println("Sample b: ", s1, " ", s2, " ")

		# Get inner sample bound coordinates
		r1 = mapRange(1.0, R(n), a, b, R(s1))
		r2 = mapRange(1.0, R(n), a, b, R(s2))

		# Bounds within same sample bin
		if s1 > s2
			f1 = f.coeffs[s2]
			f2 = f.coeffs[s2 + 1]
			dx1 = x1 - r2
			dx2 = x2 - x1
			dx3 = r1 - x2
			m = (f2 - f1)/h
			u = f1 + m*dx1
			v = f2 - m*dx3
			A = trap(dx2, u, v)
			# println("Intra A: ", A)
			A
		# Bounds in different sample bins (typical case)
		else
			# First partial bin
			A1 = R(0)
			if s1 > 1 && s1 <= n
				f1 = f.coeffs[s1 - 1]
				f2 = f.coeffs[s1]
				dx1 = x1 - r1 + h
				dx2 = r1 - x1
				m = (f2 - f1)/h
				u = f1 + m*dx1
				A1 = trap(dx2, u, f2)
			end
			# Full bins
			A2 = R(0)
			if s1 == s2
				# 0 bins => 0 evals
			else
				# s2 - s1 bins => s2 - s1 + 1 evals
				for i in s1:s2
					# weight of each basis function
					# boundary basis functions have half contribution
					w = (i == s1 || i == s2 ? R(1//2) : R(1))
					A2 += w * h * f.coeffs[i]
				end
			end
			# Last partial bin
			A3 = R(0)
			if s2 < n && s1 >= 1
				f3 = f.coeffs[s2]
				f4 = f.coeffs[s2 + 1]
				dx3 = x2 - r2
				dx4 = r2 + h - x2
				m = (f4 - f3)/h
				v = f4 - m*dx4
				A3 = trap(dx3, f3, v)
			end
			# println("First A:  ", A1)
			# println("Middle A: ", A2)
			# println("Last A:   ", A3)
			A1 + A2 + A3
		end
	end
end

#
# Norm
#

"""
Compute the L2 norm of a function in the hat basis representation.
"""
function normL2(f::Hat{R, a, b, n}) where {R, a, b, n}
	fabs2 = Hat{R, a, b, n}([abs2(c) for c in f.coeffs])
	sqrt(quad(fabs2))
end

#
# Derivatives
#

"""
Compute the derivative of a function in the hat basis representation and
project an extension into the hat basis representation.

CAVEAT: Technically the derivative is a discontinuous (piecewise constant)
function and it's undefined at the sample points! In order to project into the
hat basis representation on the same sample points, we must extend the domain
of the derivative as follows:
- At the sample midpoints, take the average of the limit of the derivative on
either side.
- At the sample endpoints, take the inner limit of the derivative.
"""
function deriv(f::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	h = (b - a) / (n - 1)
	dcoeffs = similar(f.coeffs)
	for i in 1:n
		if i == 1
			dcoeffs[1] = (f.coeffs[2] - f.coeffs[1]) / h
		elseif i == n
			dcoeffs[n] = (f.coeffs[n] - f.coeffs[n - 1]) / h
		else
			dcoeffs[i] = (f.coeffs[i + 1] - f.coeffs[i - 1]) / 2h
		end
	end
	Hat{R, a, b, n}(dcoeffs)
end

"""
Compute an approximation of the second derivative of a function in the hat
basis representation. (This is not taking the projection of the actual second
derivative, which is zero almost everywhere). In principal, this does two
applications of derivative but with less \"information loss\".

CAVEAT: Same as deriv
"""
function deriv2(f::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	h = (b - a) / (n - 1)
	dcoeffs = similar(f.coeffs)
	for i in 1:n
		if i == 1
			dcoeffs[1] = (f.coeffs[3] - 2*f.coeffs[2] + f.coeffs[1]) / h^2
		elseif i == n
			dcoeffs[n] = (f.coeffs[n] - 2*f.coeffs[n - 1] + f.coeffs[n - 2]) / h^2
		else
			dcoeffs[i] = (f.coeffs[i + 1] - 2*f.coeffs[i] + f.coeffs[i - 1]) / h^2
		end
	end
	Hat{R, a, b, n}(dcoeffs)
end

#
# Miscellaneous
#

function bounds(f::Hat{R, a, b, n})::Tuple{Int, Int} where {R, a, b, n}
	(a, b)
end

function samples(f::Hat{R, a, b, n})::Int where {R, a, b, n}
	n
end

function concat(f::Hat{R, a, b, m}, g::Hat{R, b, c, n})::Hat{R, a, c, m + n - 1} where {R, a, b, c, m, n}
	# TODO: floating point error
	if f.coeffs[m] != g.coeffs[1]
		throw("Last and first sample do not match.")
	end
	Hat{R, a, c, m + n - 1}(vcat(f.coeffs, g.coeffs[2:end]))
end
