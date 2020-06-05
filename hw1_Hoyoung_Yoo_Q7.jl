# WRITTEN BY HOYOUNG YOO 05-30-2020 REFERENCE GARRETT_CODE
# SUMMER 2020, COMPUTATIONAL BOOTCAMP PS 1 - Q #7 Capital Investment

cd("C:\\Users\\stell\\Dropbox\\Hoyoung\\Coursework\\10_2020summer\\Computational\\02homework\\hw01")

using Parameters, Plots #import the libraries we want
include("hw1_Hoyoung_Yoo_Q7_model.jl") #import the functions that solve our growth model

prim, res = Initialize()
@elapsed prim, res = Solve_model(prim::Primitives, res::Results)
@unpack val_func, pol_func = res
@unpack k_grid = prim

### Draw Plots
# value function
plot(k_grid, val_func, title = "Value Functions", label = ["Good State" "Bad State"])
savefig("Question7_value_functions.png")

# policy function
plot(k_grid, pol_func, title = "Policy Functions", label = ["Good State" "Bad State"])
savefig("Question7_policy_functions.png")
