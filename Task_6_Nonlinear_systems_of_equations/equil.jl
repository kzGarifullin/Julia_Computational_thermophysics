using LinearAlgebra
function broydensys(f, x, J; maxiter=50, xtol=1e-6, ftol=1e-6)
    δx = float(similar(x))
    yp, yn = similar.((δx, δx))
    x = float(copy(x))
    B = float(copy(J))
    yn .= f(x)
    for i in 1:maxiter
        yp .= yn
        δx .= .- (B \ yp)
        x .+= δx
        yn .= f(x)

        norm(δx) < xtol && return x
        norm(yn) < ftol && return x
        
        g = B * δx
        B .+= (1 / dot(δx, δx)) .* (yn .- yp .- g) .* δx'
    end
    error("Превышено число итераций.")
end

function jacobianfd(f, x; y=f(x), δ=sqrt(eps())*max(norm(x), 1))
    m, n = size(y, 1), size(x, 1)
    J = zeros(m, n)
    x = float(copy(x))
    for j in 1:n
        x[j] += δ
        J[:, j] .= (f(x) .- y) ./ δ
        x[j] -= δ
    end
    return J
end


P_r(V_r::Float64, T_r::Float64) = 8*T_r/(3*V_r - 1) - 3 / V_r / V_r
mu(V_r::Float64, T_r::Float64)  = -T_r * log((3*V_r - 1)/(2*exp(-0.5))) + T_r/(3*V_r - 1) - 9 / (4*V_r)

function f(x, T_r)
    V_g, V_l = x
    return [
        P_r(V_g, T_r) - P_r(V_l, T_r),
        mu(V_g, T_r) - mu(V_l, T_r),
    ]
end
#F1(z::Array{Float64}, K::Array{Float64}) = G -> F(G::Float64 ,z::Array{Float64}, K::Array{Float64})
f(T_r) = x -> f(x, T_r)
T_range = [0.85, 0.88, 0.9, 0.93, 0.95, 0.97, 0.99] 

for k in 1:7
    T = T_range[k]
    println("T: ", T)
    x_left = [0.5, 0.8]
    ts_l = (x_left[2] - x_left[1]) / 7.0
    x_right = [1.2, 3.1]
    ts_r = (x_right[2] - x_right[1]) / 7.0
    x = [x_left[1] + ts_l * (k-1), x_right[2] - ts_r * (k-1)]
    #println(x)
    #println(jacobianfd(f, x))
    
    rfdjac = broydensys(f(T), x, jacobianfd(f(T), x), maxiter=50)
    println(rfdjac)
    
end

