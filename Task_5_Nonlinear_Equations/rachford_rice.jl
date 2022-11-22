using Plots
function bisection(f, x₁, x₂; xtol=eps(), ftol=eps())
    if x₁ > x₂; x₁, x₂ = x₂, x₁; end
    y₁, y₂ = f(x₁), f(x₂)

    sign(y₁) == sign(y₂) && error("Функция должна иметь разные знаки в концах отрезка")
    abs(y₁) < ftol && return x₁
    abs(y₂) < ftol && return x₂
    
    maxiter = ceil(Int, log2((x₂-x₁)/xtol))
    
    for i in 1:maxiter
        xnew = (x₂ + x₁) / 2
        ynew = f(xnew)
        
        if sign(y₂) == sign(ynew)
            x₂, y₂ = xnew, ynew
        elseif sign(y₁) == sign(ynew)
            x₁, y₁ = xnew, ynew
        else
            return xnew
        end
        abs(ynew) < ftol && return xnew
    end
    return (x₂ + x₁)/2
end

function ridders(f, x₁, x₂; maxiter=25, xtol=eps(), ftol=eps())
    if x₁ > x₂; x₁, x₂ = x₂, x₁; end
    y₁, y₂ = f.((x₁, x₂))
    y₁ * y₂ > 0 && error("Функция должна иметь разные знаки в концах отрезка")
    y₁ == 0 && return x₁
    y₂ == 0 && return x₂
    
    for i in 1:maxiter
        xmid = (x₁ + x₂) / 2
        ymid = f(xmid)
        xnew = xmid + (xmid - x₁) * sign(y₁) * ymid / sqrt(ymid^2 - y₁*y₂)
        ynew = f(xnew)

        ynew == 0 && return xnew
        
        if sign(ynew) == sign(y₂)
            x₂, y₂ = xnew, ynew
        elseif sign(ynew) == sign(y₁)
            x₁, y₁ = xnew, ynew
        end
        
        abs(ynew) < ftol && return xnew
        abs(x₁ - x₂) < xtol && return xnew
    end
    error("Число итераций превышено.")
end

function F(G::Float64 ,z::Array{Float64}, K::Array{Float64})
    l = length(z)
    F = 0
    for i in 1:l
        #println(i)
        F = F + z[i] * (K[i] - 1) / (G*(K[i] - 1) + 1)
    end
    #println(l)   
    return F
end

function plot_graph(foo, ts)
    #scatter(ts, ys; label="Узлы интерполяции", legend=:top, xlabel="Temperature, K", ylabel="Viscosity, μPa * s")
    xs = range(first(ts), last(ts); length=200)
    plot(xs, foo.(xs); title = "Molecule: ", label="Кубический сплайн")
end



z = [0.2463, 0.2208, 0.2208, 0.3121]
K = [40, 25, 0.6, 0.005]
F1(z::Array{Float64}, K::Array{Float64}) = G -> F(G::Float64 ,z::Array{Float64}, K::Array{Float64})
F_new = F1(z, K)



root = bisection(F_new, 1.0/(1 - maximum(K)), 1.0/(1 - minimum(K))) 
println(root)

z1 = [0.9, 0.1]
K1 = [1.5, 0.01]
F_new_1 = F1(z1, K1)
root1 = bisection(F_new_1, 1.0/(1 - maximum(K1)), 1.0/(1 - minimum(K1))) 
println(root1)
#display(plot_graph(F_new, [1.0/(1 - maximum(K)), 1.0/(1 - minimum(K))]))
#println(F(6.0, z, K))
