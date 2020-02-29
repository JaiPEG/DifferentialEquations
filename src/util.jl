"""
Evaluate the unique affine map between the given intervals at the given
point.

[a1, b1] : Old interval
[a2, b2] : New interval
x : Old point
-> : New point
"""
function mapRange(a1::T, b1::T, a2::T, b2::T, x::T)::T where {T}
	(x - a1)*(b2 - a2)/(b1 - a1) + a2
end

"""
Clamp the given point to the given interval.

If x < a, a
Elseif x > b, b
Else, x
"""
function clamp(a::T, b::T, x::T)::T where {T}
	if x < a
		a
	elseif x > b
		b
	else
		x
	end
end

"""
Return [e, f(e), ..., f^{n-1}(e)].
"""
function accum(f, e::T, n::Int)::Vector{T} where {T}
	if n == 0
		T([])
	else
		ret = Vector{T}(undef, n)
		ret[1] = e
		for i in 2:n
			ret[i] = f(ret[i - 1])
		end
		ret
	end
end

"""
Equivalent to accum(f, e, n + 1)[end] but more efficient.
Note the n + 1.
"""
function apply(f, e::T, n::Int)::T where {T}
	if n == 0
		e
	else
		ret = e
		for i in 1:n
			ret = f(ret)
		end
		ret
	end
end

"""
Return the binary representation of an integer as a vector of integers from
0 to 1 where the first element of the vector is the least significant bit.
"""
function binaryLSB(n::Int)::Vector{Int}
	if n < 0
		throw(DomainError(n, "Negative input."))
	elseif n == 0
		Int[]
	else
		q, r = divrem(n, 2)
		vcat([r], binaryLSB(q))
	end
end

"""
Return the binary representation of an integer as a vector of integers from
0 to 1 where the first element of the vector is the most significant bit.
"""
function binaryMSB(n::Int)::Vector{Int}
	reverse(binaryLSB(n))
end

"""
Return the binary representation of an integer as a vector of integers from
0 to 1. Alias to binaryMSB.
"""
binary = binaryMSB

"""
Fast powering algorithm for any associative binary operation f with identity
e.
"""
function pow(f, e::T, x::T, n::Int)::T where {T}
	# fast powering algorithm
	mask = binaryLSB(n)
	squares = accum(x -> f(x, x), x, length(mask))
	ret = e
	for (bit, square) in zip(mask, squares)
		if bit == 1
			ret = f(ret, square)
		end
	end
	ret
end
