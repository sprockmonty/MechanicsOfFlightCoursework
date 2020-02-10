using Revise
using LinearAlgebra
using ControlSystems
using Plots
using Printf
plotly()
# constant definitions

ρ_∞ = 1.112 
θ_0 = -2
V_∞ = 28
m = 318
S = 6.11
l_t = 4.63
I_yy = 432
S_t = 1.14
l_w = 0.3
C_Lαw = 5.55
C_Lαt = 4.3
C_Lδt = 0.45
b = 15
T_0 = 0
e = 0.8
g = 9.81

C_Lα = C_Lαw + S_t * C_Lαt / S
q_∞ = 0.5 * ρ_∞ * V_∞^2

# Question 1
# state matrix
A = [ 
      [
      2 * ( m * g * sind(θ_0) -T_0) / (m * V_∞)
      g * cosd(θ_0) * (1 - 2 * S * C_Lα / (π * b^2 * e) ) / V_∞
      0
      -g * cosd(θ_0)
      ]';

      [
      -2 * g * cosd(θ_0) / V_∞
      (m * g * sind(θ_0) - T_0) / (m * V_∞) - ρ_∞ * V_∞ * S * C_Lα / (2 * m)
      V_∞ - ρ_∞ * V_∞ * S_t * l_t * C_Lαt / (2 * m)
      -g * sind(θ_0)
      ]';

      [
      0
      ρ_∞ * V_∞ * (S * l_w * C_Lαw - S_t * l_t * C_Lαt) / (2 * I_yy)
      -ρ_∞ * V_∞ * S_t * l_t^2 * C_Lαt / (2 * I_yy)
      0
      ]';
     
      [
      0
      0
      1
      0
      ]'
    ] 

# input matrix
 B = [ 
        0
        -q_∞ * S_t * C_Lδt / m
        -q_∞ * S_t * l_t * C_Lδt / I_yy
        0
     ]

# sensory matrix
C = Matrix(I,4,4)

# direct term
D = 0

# eigenvalue plots
λ = eigvals(A)
scatter(λ, xlabel="Re(λ)", ylabel="Im(λ)")

# Question 2
# state space
sys = ss(A,B,C,D)
stepplot(sys, 200)
phugω = imag(λ[4]) # phugiod frequency
@printf("The phugoid frequency is %f and the time period is %f", phugω, 2 * π / phugω)


# Question 4

