# WRITTEN BY HOYOUNG YOO 06-03-2020
# SUMMER 2020, COMPUTATIONAL BOOTCAMP PS 2

cd("C:\\Users\\stell\\Dropbox\\Hoyoung\\Coursework\\10_2020summer\\Computational\\02homework\\hw02")

### use packages
using Optim, Random, Distributions, Plots, Parameters, ForwardDiff, LinearAlgebra, SharedArrays

## (Question 1) :: HIMMELBLAU FUNCTION

##(a) Draw a surface plot

function Himmelblau(vec::Vector{Float64})
    x, y = vec[1], vec[2]
    val = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
    val
end

x_grid = collect(-4:0.01:4)
y_grid = collect(-4:0.01:4)
nx = length(x_grid)
ny = length(y_grid)
z_grid = zeros(nx, ny)

for i = 1:nx, j = 1:ny
    z_grid[i,j] = Himmelblau([x_grid[i], y_grid[j]])
end

Plots.surface(x_grid, y_grid, z_grid)
savefig("Question1.png")

#four local minima appear

##(b) Newton's method

#Gradient
function g(G, vec::Vector{Float64})
    x, y = vec[1], vec[2]
    G[1] = 2*(x^2 + y - 11)*2*x + 2*(x + y^2 - 7)
    G[2] = 2*(x^2 + y - 11) + 2*(x + y^2 - 7)*2*y
    G
end

#Hessian
function h(H, vec::Vector{Float64})
    x, y = vec[1], vec[2]
    H[1,1] = 4*(x^2 + y - 11) + 8*x + 2
    H[1,2] = 4*x + 4*y
    H[2,1] = 4*x + 4*y
    H[2,2] = 2 + 4*(x + y^2 - 7) + 8*y^2
    H
end

#guess
guess = [0.0, 0.0]
guess = [1.0, 1.0]
guess = [-1.0, -1.0]
opt_newton = optimize(Himmelblau, g, h, guess)

##(c) Nelder-Mead
opt_nm = optimize(Himmelblau, guess)

## (Question 2) :: ACKLEY FUNCTION

function Ackley(vec::Vector{Float64})
    x, y = vec[1], vec[2]
    val = -20*exp(-0.2*sqrt(0.5*(x^2 + y^2))) - exp(0.5*(cos(2*π*x) + cos(2*π*y))) + ℯ + 20
end

##(a) Draw a surface plot

z_grid = zeros(nx, ny)

for i = 1:nx, j = 1:ny
    z_grid[i,j] = Ackley([x_grid[i], y_grid[j]])
end

Plots.surface(x_grid, y_grid, z_grid)
savefig("Question2_a.png")
Plots.contourf(x_grid, y_grid, z_grid)
savefig("Question2_b.png")

#global minima 0.0 at (0.0, 0.0)

##(b) Minimize the function using LBFGS and Nelder-Mead

guess = [1.0, 0.0]

opt_lbfgs = optimize(Ackley, guess, LBFGS())
opt_nm = optimize(Ackley, guess)

## (Question 3) :: RASTRIGIN FUNCTION

function Rastrigin(x::Array{Float64,1})
    A = 10
    n = length(x)
    val = A*n + sum(x.^2 .- A*cos.(2*π*x))
end

##(a) n = 1
xgrid = collect(-5.12:0.01:5.12)
nx = length(xgrid)
zgrid = zeros(nx)

for i = 1:nx
    x = [xgrid[i]]
    zgrid[i] = Rastrigin(x)
end

plot(xgrid,zgrid)
savefig("Question3.png")

guess = [1.0]
opt = optimize(Rastrigin, guess, LBFGS())

##(b) n = 2
xgrid1 = collect(-5.12:0.01:5.12)
nx1 = length(xgrid1)
xgrid2 = collect(-5.12:0.01:5.12)
nx2 = length(xgrid2)

zgrid = zeros(nx1, nx2)

for i = 1:nx1
    for j = 1:nx2
        zgrid[i, j] = Rastrigin([xgrid1[i], xgrid2[j]])
    end
end

Plots.surface(xgrid1, xgrid2, zgrid)
savefig("Question3_a.png")
Plots.contourf(xgrid1, xgrid2, zgrid)
savefig("Question3_b.png")

##(c) LBFGS & Nelder-Mead
guess = [1.0, 1.0]

opt_lbfgs = optimize(Rastrigin, guess, LBFGS())
opt_nm = optimize(Rastrigin, guess)

## (Question 4) :: LINEAR APPROXIMATION

a = 0.01 ;
b = 100 ;
n = 1000 ;
x = 50 ;

function random(x::Float64)
    val = (x-a)^2 + b
end

function linearapprox(a,b,n,x)

    X = collect(range(a, length = n, stop = b))
    X_after = X[findfirst(z->z>x, X)]
    X_before = X[findfirst(z->z>x, X)-1]
    m = (random(X_after)-random(X_before))/(X_after-X_before)
    y_intercept = random(X_after) - m*X_after
    y = m*x + y_intercept
    println("linear approximation = ", y, ", true value = ", f(X_after))

end

## (Question 5) :: OPTIMAL APPROXIMATION

function lfunction(x::Float64)
    val = log(x+1)
end

# create linear approximation function
function linearapprox(a,b,n,x)

    X = collect(range(a, length = n, stop = b))
    X_after = X[findfirst(z->z>x, X)]
    X_before = X[findfirst(z->z>x, X)-1]
    m = (lfunction(X_after)-lfunction(X_before))/(X_after-X_before)
    y_intercept = lfunction(X_after) - m*X_after
    approx = m*x + y_intercept

end

# set values
a = 0 ;
b = 100 ;
n = 11 ;

# (A) plot
x_candidates_A = collect(a:0.01:b)
n_A = length(x_candidates_A)
plot_A = zeros(n_A-1, 3)

for i in 1:n_A-1
    plot_A[i,1] = linearapprox(a, b, n, x_candidates_A[i]) # approx
    plot_A[i,2] = lfunction(x_candidates_A[i]) # true
    plot_A[i,3] = abs(linearapprox(a, b, n, x_candidates_A[i]) - lfunction(x_candidates_A[i])) # diff
end

plot(x_candidates_A[1:n_A-1], plot_A[:,1], label = "linear approximation")
plot!(x_candidates_A[1:n_A-1], plot_A[:,2], label = "true value")
plot!(x_candidates_A[1:n_A-1], plot_A[:,3], label = "difference")
savefig("Question5_a.png")

# (B)
residual_sum = sum(plot_A[:,3])

# (C) Nelder-Mead :: need to be fixed
function residual(sol::Array{Float64,9})
    a = 0
    b = 100
    X = [a, x_candidates_A[sol], b]
    residual_sum = 0
    for i in 1:length(X)-1
        residual = abs(linearapprox(a, b, n, X[i]) - lfunction(X[i]))
        residual_sum = residual_sum + residual
    end
    println(residual_sum)
end

guess = [10 20 30 40 50 60 70 80 90]
opt = optimize(residual, guess)

## (Question 5) :: AGAIN

function approx_log(grid::Vector{Float64})
    n = length(grid)
    func_grid = zeros(n)
    f(x) = log(1+x)

    for i = 1:n
        func_grid[i] = f(grid[i])
    end

    grid_interp = interpolate(grid, BSpline(Linear()))
    func_interp = interpolate(func_grid, BSpline(Linear()))
    errs, func_approx, func_true = zeros(1001), zeros(1001), zeros(1001)

    for i = 0:1000
        x = i/10
        findmatch(index) = abs(grid_interp(index) - x)
        x_ind = optimize(index->findmatch(index), 1.0, n).minimizer
        func_approx[i+1] = func_interp(x_ind)
        errs[i+1] = abs(f(x) - func_interp(x_ind))
        func_true[i+1] = f(x)
    end
    errs, func_approx, func_true
end

x_grid = collect(0:0.1:100)
grid_approx = collect((0.0:10.0:100.0))
errs, func_approx, func_true = approx_log(grid_approx)
sum(errs)
Plots.plot(x_grid, errs)
Plots.plot(x_grid, [func_approx, func_true])

function eval_points(guess::Vector{Float64})

    if minimum(guess)<0 || maximum(guess)>100
        return Inf
    else
        grid = vcat(0.0, guess, 100.0)
        errs, func_approx, func_true = approx_log(grid)
        return sum(errs)
    end
end

guess_init = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]
opt = optimize(guess->eval_points(guess), guess_init; g_tol = 1e-4)

grid_approx = vcat(0.0, opt.minimizer, 100.0)
errs, func_approx, func_true = approx_log(grid_approx)

Plots.plot(x_grid, errs)
Plots.plot(x_grid, [func_approx, func_true])
sum(errs)
