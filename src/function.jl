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

function Base.:+(f::Hat{R, a, b, n}, g::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(f.coeffs + g.coeffs)
end

function Base.:-(f::Hat{R, a, b, n}, g::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(f.coeffs - g.coeffs)
end

function Base.:*(c::R, g::Hat{R, a, b, n})::Hat{R, a, b, n} where {R, a, b, n}
	Hat{R, a, b, n}(c * f.coeffs)
end

#
# Projection
#

"""
Sample a function [a, b] -> R at n points uniformly spaced between a and b.
"""
function sampleUniform(a::R, b::R, n::Int, f)::Vector{R} where {R}
	h = R(1) / (n - 1)
	[f(h * i) for i in 0:(n - 1)]
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
		h = 1 / R(n - 1)
		xip = h * (i - 1)
		xi  = h * i
		xin = h * (i + 1)
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

#
# Norm
#

"""
Compute the L2 norm of a function in the hat basis representation.
"""
function normL2(f::Hat{R, a, b, n}) where {R, a, b, n}
	fabs2 = Hat{R, a, b, n}([abs2(c) for c in f.coeffs])
	sqrt(quadrature(fabs2))
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

function concat(f::Hat{R, a, b, m}, g::Hat{R, b, c, n})::Hat{R, a, c, m + n - 1} where {R, a, b, c, m, n}
	# TODO: floating point error
	if f.coeffs[m] != g.coeffs[1]
		throw("Last and first sample do not match.")
	end
	Hat{R, a, c, m + n - 1}(vcat(f.coeffs, g.coeffs[2:end]))
end
