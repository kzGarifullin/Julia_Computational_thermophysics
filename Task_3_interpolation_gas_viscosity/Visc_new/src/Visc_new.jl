module Visc_new

using LinearAlgebra
export Substance, viscos
struct Substance{T}
    M::T
    sigma::T
    k::T
    T_min::T
    T_max::T
end
function horner(c, x)
    ans = last(c)
    for i in lastindex(c)-1:-1:1
        ans = (ans * x) + c[i]
    end
    return ans
end
polynomial(c) = x -> horner(c, x)
function spinterp(t, y)
    n = size(t, 1) - 1
    
    In = I(n)
    E = In[1:end-1, :]
    J = diagm(0 => ones(n), 1 => -ones(n-1)) 
    Z = zeros(n, n)
    h = [t[k+1] - t[k] for k in 1:n]
    H = diagm(0 => h)
    
    # 1.а Значения на левой границе
    AL = [In Z Z Z]
    vL = y[1:end-1]
    
    # 1.б Значения на правой границе
    AR = [In H H^2 H^3]
    vR = y[2:end]
    
    # 2. Непрерывность первой производной
    A1 = E * [Z J 2*H 3*H^2]
    v1 = zeros(n-1)

    # 3. Непрерывность второй производной
    A2 = E * [Z Z J 6*H]
    v2 = zeros(n-1)
    
    # 4. Not-a-knot
    nakL = [zeros(1, 3*n) 1 -1 zeros(1, n-2)]  # слева
    nakR = [zeros(1, 3*n) zeros(1, n-2) 1 -1]  # справа
    
    # Собираем систему и решаем
    A = [AL; AR; A1; A2; nakL; nakR]
    v = [vL; vR; v1; v2; 0; 0]
    coefs = A \ v
    
    # Разбираем коэффициенты
    a = coefs[1:n]
    b = coefs[n+1:2*n]
    c = coefs[2*n+1:3*n]
    d = coefs[3*n+1:4*n]
    
    S = [polynomial([a[k], b[k], c[k], d[k]]) for k in 1:n]
    
    return function (x)
        if x < first(t) || x > last(t)
            return NaN
        elseif x == first(t)
            return first(y)
        else
            k = findlast(x .> t)  # k такое, что x ∈ (tₖ₋₁, tₖ)
            return S[k](x - t[k])
        end
    end
end

function viscos(T, struc::Substance)  
    T_zv = T / struc.k
    eta = 266.93 /100 * sqrt(T * struc.M) * inter_small(T/struc.k) / struc.sigma / struc.sigma / inter_om(T/struc.k) 
    return eta
end

f_eta = [1.0014, 1.0002, 1.0000, 1.0000, 1.0001, 1.0004, 1.0014, 1.0025, 1.0034, 1.0049, 1.0058, 1.0075, 1.0079, 1.0080, 1.0080]
T_list = [0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0, 50.0, 100.0, 400.0]
omega = [2.785, 2.628, 2.492, 2.368, 2.257, 2.156, 2.065, 1.982, 1.906, 1.841, 1.780, 1.725, 1.675, 1.627, 1.587, 1.549, 1.514, 1.482, 1.452, 1.424, 1.399, 1.375, 1.353, 1.333, 1.314, 1.296, 1.279, 1.264, 1.248, 1.234, 1.221, 1.209, 1.197, 1.186, 1.175, 1.156, 1.138, 1.122, 1.107, 1.093, 1.081, 1.069, 1.058, 1.048, 1.039, 1.030, 1.022, 1.014, 1.007, 0.9999, 0.9932, 0.9870, 0.9811, 0.9755, 0.9700, 0.9649, 0.9600, 0.9553, 0.9507, 0.9464, 0.9422, 0.9382, 0.9343, 0.9305, 0.9269, 0.8963, 0.8727, 0.8538, 0.8379, 0.8242, 0.7432, 0.7005, 0.6718, 0.6504, 0.6335, 0.6194, 0.6076, 0.5973, 0.5882, 0.5320, 0.5016, 0.4811]
T_om = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90, 5.00, 6.0 ,7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 ,90.0, 100.0, 200.0, 300.0, 400.0]
inter = spinterp(T_list, f_eta)
inter_om = spinterp(T_om, omega)
inter_small = spinterp(T_list[1:12], f_eta[1:12])


end # module Visc_new
