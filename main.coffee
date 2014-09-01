# See Lloyd N. Trefethen, "Spectral Methods in MATLAB"
# http://people.maths.ox.ac.uk/trefethen/spectral.html
# p27.m - Solve KdV eq. $u_t + uu_x + u_{xxx} = 0$ on [-pi,pi] by
#         FFT with integrating factor v = exp(-ik^3t)*u-hat.

# Set up grid
N = 256
dt = 0.4/N.pow(2)  # 0.4 in Trefethen
x = 2*pi/N * linspace(-N/2, N/2-1, N)

# Integrating factor
k = linspace(0, N/2-1, N/2)
    .concat([0])
    .concat(linspace(-N/2+1, -1, N/2-1)) 
ik3 = j*k.pow(3)

# FFT helper functions
zeros = (x) -> nm.zeros(1, x.length)[0]  # Zero vector same length as x
r2c = (x) -> complex x, zeros(x)  # Real to complex (for fft input)
fft = (x) -> r2c(x).fft()
ifftRe = (x) -> x.ifft().x  # Real part of ifft 

# 4th order Runge-Kutta
g = -0.5*j*dt*k
E = exp(dt*ik3/2)
E2 = E.pow(2)
tx = (g, x) -> g * fft(ifftRe(x).pow 2)
computeUV = (v) ->
    a = tx g, v
    b = tx g, E*(v + a/2)
    c = tx g, E*v + b/2
    d = tx g, E2*v + E*c
    v = E2*v + (E2*a + 2*E*(b+c) + d)/6
    u = ifftRe v
    {u, v}

# Hyperbolic secant
sech = (x) -> 2 / (exp(x) + exp(-x))

# Initial condition for soliton
soliton = (A, x1) ->
    3*A.pow(2) * (sech(.5*A*(x+x1))).pow(2)
    
lineChart = new $blab.LineChart
    id: "solitons"
    xLabel: "x"
    yLabel: "u"
    xLim: [-pi, pi]
    yLim: [0, 1000]
    xTicks: 7
    yTicks: 5
    click: (x, y) -> initSoliton(x, y) 

@kdv = new $blab.BasicAnimation
    numSnapshots: 2000
    delay: 10
    strobeInterval: 10000
    snapshotFunction: ->
    strobeFunction: ->

# Plot control
count = 0 # number of solitons in plot
u0 = null # initial condition
initSoliton = (xS, yS) ->
    @kdv.stopAnimation()
    count++
    if count is 1
        u0 = zeros(x)
        @kdv.n = 0
    x1 = -xS
    A = sqrt(yS/3)
    u0 += soliton(A, x1)
    lineChart.plot(x, u0)
    if count is 2
        count = 0
        animateSolitons(u0)

# Animated soliton plot
animateSolitons = (u) ->
    v = fft u
    @kdv.snapshotFunction = ->
        {u, v} = computeUV(v)
        lineChart.plot(x, u)
    @kdv.animate()
    
# Run: two solitons
initSoliton(-1, 800)
initSoliton(0, 200)
