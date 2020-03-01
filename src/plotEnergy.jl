using Gadfly, Printf

include("util.jl")
include("function.jl")
include("ode.jl")
include("pde.jl")

# Discretization parameters
a = 0.0
b = 2.0*pi
ns = [10, 20, 40, 80, 160, 320, 640, 1280]

# Simulation parameters
t1 = 0.0
t2 = 100.0
c = 1.0

Ess = Vector(undef, length(ns))
Esf = Vector{Function}(undef, length(ns))

@printf("Space interval: %f, %f\n", a, b)
@printf("Time interval:  %f, %f\n", t1, t2)
for (i, n) in enumerate(ns)
	# Discretization parameters
	h = (b - a)/(n - 1)

	# Simulation parameters
	ht = 0.1*h # make time-step proportional to resolution
	nt = Int(round((t2 - t1)/ht)) + 1
	ht = (t2 - t1)/(nt - 1) # recalculate time-step

	# Initial conditions
	f0 = projectHat(a, b, n, sin)
	fx0 = deriv(f0)
	# ft0 = Hat{Float64, a, b, n}(rand(n)) # ?
	# ft0 = Hat{Float64, a, b, n}(zeros(n)) # Standing wave
	ft0 = fx0 # Interesting dynamics
	s0 = [ft0, fx0]

	# Plotting
	s = s0
	Ess[i] = Vector{Float64}(undef, nt)

	@printf("Using %d samples and taking %d time-steps.\n", n, nt)
	for k in 1:nt
		s = waveD0(rk2, s, ht)
		ft, fx = s
		Ess[i][k] = 0.5*quad(ft^2 + c^2*fx^2)
	end
	# Use Hat since it's eval implementation does linear interpolation
	Esf[i] = x -> eval(Hat{Float64, t1, t2, nt}(Ess[i]), x)
end

p = plot(Esf, t1, t2,
	Guide.xlabel("Time"), Guide.ylabel("Energy"),
	Guide.colorkey(title="Samples", labels=[string(n) for n in ns])
)
fileName = "output/plotEnergy/energy.svg"
draw(SVG(fileName, 6inch, 4inch), p)
println("Saved to " * fileName * ".")
