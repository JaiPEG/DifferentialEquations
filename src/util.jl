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
