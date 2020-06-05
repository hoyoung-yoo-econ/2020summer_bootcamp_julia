# Types with default values & keyword constructors
@with_kw struct Primitives
    θ::Float64 = 0.36 #capital share
    δ::Float64 = 0.025 #depreciation rate
    β::Float64 = 0.99 #discount rate
    K_min::Float64 = 0.01 #capital lower bound
    K_max::Float64 = 45 #capital upper bound
    nk::Int64 = 1000 #number of capital grid points
    Π::Array{Float64,2} = [0.977 0.023; 0.074 0.926] #markov process
    states::Array{Float64,2} = [1.25 0.2] #states of productivity
    nz::Int64 = length(states) #number of markov states
    k_grid::Array{Float64,1} = collect(range(K_min, length = nk, stop = K_max)) #capital grid
end

# Structure that holds model results
mutable struct Results
    val_func::Array{Float64,2} #value function
    pol_func::Array{Float64,2} #policy function
end

#  function for initializing model primitives and results
function Initialize()
    prim = Primitives() #initialize primitives
    val_func = zeros(prim.nk, prim.nz) #initial value function guess
    pol_func = zeros(prim.nk, prim.nz) #initial policy function guess
    res = Results(val_func, pol_func) #initialize results struct
    prim, res #return deliverables
end

# Bellman Operator
function Bellman(prim::Primitives, res::Results)
    @unpack val_func, pol_func = res #unpack value function
    @unpack k_grid, nz, nk, β, δ, θ, Π, states = prim #unpack model primitives
    v_next = zeros(nk, nz) #next guess of value function to fill

    for z = 1:nz #productivity state space
        for i_k = 1:nk #capital state space
            k = k_grid[i_k]
            budget = states[z]*k^θ + (1-δ)*k
            max_utility = -Inf

            for i_kp = 1:nk
                kp = k_grid[i_kp]
                c = budget - kp
                if c>0
                    val = log(c) + β*sum(Π[z,:].*val_func[i_kp,:])
                    if val > max_utility
                        max_utility = val
                        res.pol_func[i_k,z] = kp
                    end
                end
            end
            v_next[i_k,z] = max_utility
        end
    end
    v_next
end

#Solve the model
function Solve_model()
    n = 0;
    error, tol = 1000, 0.001;

    while error > tol
        n += 1
        v_next = Bellman(prim, res)
        error = maximum(abs.(v_next .- res.val_func))
        res.val_func = v_next
        println("Current error is ", error)
    end
    println("Value fucntion is converged in ", n, " iterations")
    prim, res
end
