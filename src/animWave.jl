using Gadfly, Printf

include("util.jl")
include("function.jl")
include("ode.jl")
include("pde.jl")

# Discretization parameters
a = 0.0
b = 2.0*pi
n = 100
h = (b - a)/(n - 1)

# Simulation parameters
t1 = 0.0
t2 = 10.0
ht = 0.1*h # make time-step proportional to resolution
nt = Int(round((t2 - t1)/ht)) + 1
ht = (t2 - t1)/(nt - 1) # recalculate time-step
frameRate = 30.0
frames = Int(round((t2 - t1)*frameRate))
frameToTime(frame) = mapRange(1.0, Float64(frames), t1, t2, Float64(frame))
timeToSteps(t) = Int(round(mapRange(t1, t2, 0.0, Float64(nt), t)))
frameToSteps(frame) = timeToSteps(frameToTime(frame))

# Initial conditions
f0 = projectHat(a, b, n, sin)
fx0 = deriv(f0)
# ft0 = Hat{Float64, a, b, n}(rand(n)) # ?
# ft0 = Hat{Float64, a, b, n}(zeros(n)) # Standing wave
ft0 = fx0 # Interesting dynamics
s0 = [ft0, fx0]

# Plotting
xs = [i for i in range(a, b, length=1000)]
# Simulate by frames
@printf("Simulation will last %d steps and produce %d frames.\n",
	frameToSteps(frames), frames)
for frame in 1:frames
	# TODO: Inefficient since simulating from s0 everytime
	k = frameToSteps(frame)
	fx = apply(s -> waveD0(rk2, s, ht), s0, k)[2]
	f = x -> quad(fx, a, x)
	p = plot(x=xs, y=[f(x) for x in xs], Geom.line,
		Scale.x_continuous(minvalue=a, maxvalue=b),
		Scale.y_continuous(minvalue=-2.1, maxvalue=2.1));
	# @sprintf expects hard-coded format string so we have to choose a fixed number of digits to print
	# Four is reasonable but this will need to be increased if frames > 9999
	fileName = @sprintf("output/animWave/frame-%04d.svg", frame)
	draw(SVG(fileName, 6inch, 4inch), p)
	print("\rSaved to " * fileName * ".")
end
println("")
