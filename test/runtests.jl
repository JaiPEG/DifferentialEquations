using Test
include("../src/util.jl")
include("../src/function.jl")
include("../src/ode.jl")
include("../src/pde.jl")

#
# Test Function
#

@testset "Function" begin
	println("Testing quad")
	println("\tx^2")
	a = 0.0
	b = 1.0
	n = 1000
	f = projectHat(a, b, n, x -> x^2)
	numTests = 100
	totErr = 0.0
	for xs in zip(rand(numTests), rand(numTests))
		x1 = xs[1]
		x2 = xs[2]
		err = abs(quad(f, x1, x2) - (x2^3 - x1^3)/3.0)
		@test err < 1e-3
		totErr += err
	end
	println("\tAverage error: ", totErr / numTests)
	println("\tsin(x)")
	a = 0.0
	b = 2.0*pi
	n = 1000
	f = projectHat(a, b, n, sin)
	numTests = 100
	totErr = 0.0
	for xs in zip(rand(numTests), rand(numTests))
		x1 = xs[1]
		x2 = xs[2]
		err = abs(quad(f, x1, x2) - (-cos(x2) + cos(x1)))
		@test err < 1e-3
		totErr += err
	end
	println("\tAverage error: ", totErr / numTests)
end
