"""
Given the autonomous initial value problem,
	y' = f(y), y(0) = y0
take a second-order Runge-Kutta step with step size h.

f :: V -> V
	Takes an ordinate vector and returns its derivative
y0 :: V
	Initial vector
h :: R
	Step size

-> :: V
	Approximation of f(y0 + h)
"""
function rk2(f, y0::V, h::R)::V where {V, R}
	k0 = f(y0)
	y1 = y0 + h/2*k0
	k1 = f(y1)
	y2 = y0 + h*k1
	y2
end
