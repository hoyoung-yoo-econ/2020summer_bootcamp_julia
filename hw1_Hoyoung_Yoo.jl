# WRITTEN BY HOYOUNG YOO 05-30-2020
# SUMMER 2020, COMPUTATIONAL BOOTCAMP PS 1

cd("C:\\Users\\stell\\Dropbox\\Hoyoung\\Coursework\\10_2020summer\\Computational\\02homework\\hw01")

### use packages
using Random, Distributions, Plots, Parameters, ForwardDiff

### Question 1 :: factorial
function myfactorial(n)
    fact = 1
    for m in 1:n
        fact *= m
    end
    fact
end

myfactorial(10)

### Question 2 :: polynomial

function p(x, a)
    poly = zero(x)
    for (i, v) in enumerate(a)
        poly += v*(x^(i-1))
    end
    poly
end

a = [1, 2, 3];
x = 2;
p(x,a)

### Question 3 :: approximation for pi

dist = Uniform(-1,1)

function pi(N::Int)
    n_landed_in_circle = 0
    for i = 1:N
        x = rand(dist)
        y = rand(dist)
        r2 = x*x + y*y
        if r2 < 1.0
            n_landed_in_circle += 1
        end
    end
    return n_landed_in_circle / N*4
end

pi(1000)

### Question 4 :: DGP

# set distribution
N = 50;
n = 20;

# parameters
a = 0.1 ;
b = 0.2 ;
c = 0.5 ;
d = 1.0 ;
σ = 0.1 ;

# simulation
x1 = randn(N)
x2 = randn(N)

beta = zeros(n, 4)

for i = 1:n
    w = randn(N)
    y = a*x1 + b*(x1.^2) + c*x2 .+ d + σ*w
    X = [x1 x1.^2 x2 ones(N)]
    beta[i,:] = (inv(X'*X)*(X'*y))'
end

# plot histograms
# a
a = histogram(beta[:,1], bins = 10, color=:blue, fillalpha=.5, label = "a")
plot!([0.1], seriestype="vline", label = false, lw = 2)
# b
b = histogram(beta[:,2], bins = 10, color=:blue, fillalpha=.5, label = "b")
plot!([0.2], seriestype="vline", label = false, lw = 2)
# c
c = histogram(beta[:,3], bins = 10, color=:blue, fillalpha=.5, label = "c")
plot!([0.5], seriestype="vline", label = false, lw = 2)
# d
d = histogram(beta[:,4], bins = 10, color=:blue, fillalpha=.5, label = "d")
plot!([1.0], seriestype="vline", label = false, lw = 2)
# altogether
plot(a, b, c, d, layout = (2,2))
savefig("Question4.png")

### Question 5 :: first-passage time

# random walk
function firstpassage(a, α)
a = 0
    t_max = 200;
    X = ones(t_max+1);
    X[t_max+1] = a;
    σ = 0.2;
    ϵ = randn(t_max+1);
    for i in 2:(t_max-1)
            X[i] = α*X[i-1] + σ*ϵ[i]
        if X[i] <= a && i <= 199
            return i-1
        elseif X[i] > a && i == 199
            return 200
        end
    end
end

firstpassage(0, 0.8)

# plot histogram
T = zeros(100,3);

for i in 1:100
    global T1
    T[i, 1] = firstpassage(0, 0.8)
    T[i, 2] = firstpassage(0, 1.0)
    T[i, 3] = firstpassage(0, 1.2)
end

# means
m1 = mean(T[:,1])
m2 = mean(T[:,2])
m3 = mean(T[:,3])

alpha1 = histogram(T[:,1], bins=(0:5:205), label = "\\alpha = 0.8")
plot!([m1], seriestype="vline", label = m1, lw = 2)
alpha2 = histogram(T[:,2], bins=(0:5:205), label = "\\alpha = 1.0")
plot!([m2], seriestype="vline", label = m2, lw = 2)
alpha3 = histogram(T[:,3], bins=(0:5:205), label = "\\alpha = 1.2")
plot!([m3], seriestype="vline", label = m3, lw = 2)
# altogether
plot(alpha1, alpha2, alpha3, layout = (3,1))
savefig("Question5.png")

### Question 6 :: Newton's method

#Question what does lambda mean?
#global in while in function

maxiter = 10000
tol = 0.001
x_0 = 30

D(f) = x -> ForwardDiff.derivative(f,x)
f(x) = (x-1)^3
f_prime = D(f)

function fixedpointmap(f, f_prime, x_0, tol, maxiter)

    X = zeros(maxiter) ;
    X[1] = x_0 - f(x_0)/f_prime(x_0)
    n = 1 ;
    d = tol +1 ;

    while d > tol
        X[n+1] = X[n] - f(X[n])/f_prime(X[n])
        d = abs(X[n+1] - X[n])
        n = n+1
        println("iteration = ", n, " norm = ", diff)
    end

    println("solution = ", X[n])
end

fixedpointmap(f, f_prime, x_0, tol, maxiter)
