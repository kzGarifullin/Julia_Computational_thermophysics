
function romberg(f, a, b; atol=1e-10, maxstep::Integer=100)
    maxstep = max(1, maxstep)  # хотя бы одно разбиение
    I = Matrix{Float64}(undef, maxstep+1, maxstep+1)
    I[1, 1] = (b - a) * (f(a) + f(b)) / 2
    for i in 2:maxstep+1
        let hc = (b - a) / 2^(i-1), np = 2^(i-2)
            I[i, 1] = I[i-1, 1] / 2 + hc * sum(f, (a + hc * (2i-1) for i in 1:np))
        end
        for k in i-1:-1:1
            I[k, i-k+1] = (2^i*I[k+1, i-k] - I[k, i-k]) / (2^i - 1)
        end
        abs(I[1, i] - I[2, i-1]) < atol && return I[1, i]
    end
    error("Точность не удовлетворена.")
end

foo2(x) = sin(100*x*exp(-x*x))
a, b, atol = 0, 3, 1e-6;

println(romberg(foo2, -1.0/3, 3))
