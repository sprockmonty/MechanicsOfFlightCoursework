using Revise
using LinearAlgebra
using ControlSystems
using Plots
plotly()
# constant definitions
const ρ_∞ = 1.112 
const θ_0 = -2
const V_∞ = 28
const m = 318
const S = 6.11
const l_t = 4.63
const I_yy = 432
const S_t = 1.14
const l_w = 0.3
const C_Lαw = 5.55
const C_Lαt = 4.3
const C_Lδt = 0.45
const b = 15
const T_0 = 0
const e = 0.8
const g = 9.81

const C_Lα = C_Lαw + S_t * C_Lαt / S

# state matrix
const A = [ 
            [
            2 * ( m * g * sind(θ_0) -T_0) / (m * V_∞)
            g * cosd(θ_0) * (1 - 2 * S * C_Lα / (π * b^2 * e) ) / V_∞
            0
            -cosd(θ_0)
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
const B = [ 
        0
        -q_∞ * S_t * C_Lδt / m
        -q_∞ * S_t * l_t * C_Lδt / I_yy
        0
    ]

# sensory matrix
const C = Matrix(I,4,4)

# direct term
const D = 0


sys = ss(A,B,C,D)
stepplot(sys, 200)
