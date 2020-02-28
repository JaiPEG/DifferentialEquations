"""
Evaluate the unique affine map between the given intervals at the given
point.

[a1, b1] : Old interval
[a2, b2] : New interval
x : Old point
-> : New point
"""
function mapRange(a1, b1, a2, b2, x)
	(x - a1)*(b2 - a2)/(b1 - a1) + a2
end

"""
Clamp the given point to the given interval.

If x < a, a
Elseif x > b, b
Else, x
"""
function clamp(a, b, x)
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
