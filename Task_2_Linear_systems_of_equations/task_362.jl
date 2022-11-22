using LinearAlgebra

#! Это велосипед
#! https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.Tridiagonal
function diagM(M::Tridiagonal)
	n,m = size(M)
	a = zeros(n)
	b = zeros(n)
	c = zeros(n)
	a[1] = 0
	c[n] = 0
	for i in 1:n
                b[i] = M[i, i]
        end
	for i in 2:n
               a[i] = M[i, i-1]
        end
	for i in 1:(n-1)
		c[i] = M[i, i+1]
	end
	return(a, b, c)
end

function tridiagsolve(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, f::Vector{Float64})
	x = zeros(length(b))
	n = length(b) 
	alpha = zeros(length(b)) 
	beta = zeros(length(b))
	alpha[2] = -c[1]/b[1]
	beta[2] = f[1]/b[1]
	k = n - 1
	for i in 2:k
		alpha[i+1] = -c[i]/(b[i] + a[i]*alpha[i])
		beta[i+1] = (-a[i]*beta[i] + f[i])/(b[i] + a[i]*alpha[i]) 
	end
	x[n] = (-a[n] * beta[n] + f[n]) / (b[n] + a[n] * alpha[n])
	for i in k:-1:1
		x[i] = alpha[i+1] * x[i+1] + beta[i+1]
	end
	return(x)
end

#! Здесь хорошо: функция переиспользована.
function tridiagsolve(A::Tridiagonal, f::Vector{Float64})
	a, b, c = diagM(A)
	return(tridiagsolve(a, b, c, f))
end

L = [4.0 0.0 0.0 0.0; 1.0 2.0 0.0 0.0; 1.0 2.0 1.0 0.0; 1.0 3.0 2.0 1.0]

b = [6.0, 8.0, 1.0, 2.0, 1.0]
a = [0.0, 3.0, 2.0, 1.0, 6.0]
c = [2.0, 6.0, 3.0, 1.0, 0.0] 
f = [1.0, 2.0, 3.0, 4.0, 5.0]
A = Tridiagonal([6 2 0 0 0; 3 8 6 0 0; 0 2 1 3 0; 0 0 1 2 1; 0 0 0 6 1])
a_vi = [0.0, 2.0, 2.0, 3.0]
b_vi = [1.0, 2.0, 3.0, 4.0]
c_vi = [2.0, 3.0, 5.0, 0.0]
f_vi = [1.0, 2.0, 3.0, 4.0]
println(tridiagsolve(a,b,c,f))
println(tridiagsolve(a_vi,b_vi,c_vi,f_vi))
println(diagM(A))
println(tridiagsolve(A, f))
