using LinearAlgebra
function forwardsub(L::AbstractMatrix, b::Vector{Float64})
        x = Vector{Float64}([])
        n,m = size(L)
        if(n != m || n != length(b))
                return x
        end
        B = LowerTriangular(L)
        if(B != L)
                return ("Error! Matrix must be LowerTriangular")
        end
        if(det(L) == 0)
                return("Error! Singular matrix")
        end
        k=0
        for i in 1:length(b)
                sum = 0
                for j in 1:k
                        sum = sum + L[i, j]*x[j]
                end
                k+=1
                push!(x, (b[i] - sum)/L[i, i])

        end
        return(x)
end
function backwardsub(U::AbstractMatrix, b::Vector{Float64})

        x = zeros(length(b))
        n,m = size(U)
        if(n != m || n != length(b))
                return("No solutions!")
        end
        B = UpperTriangular(U)
        if(B != U)
                return ("Error! Matrix must be UpperTriangular")
        end
        if(det(U) == 0)
                return("Error! Singular matrix")
        end
        k = length(b) + 1
        for i in length(b):-1:1

                sum = 0
                d = i+1
                for j in d:length(b)

                        sum = sum + U[i, j]*x[j]
                end
                k-=1
                x[i] = (b[i]-sum)/U[i,i]

        end
        return(x)
end

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
function tridiagsolve(A::Tridiagonal, f::Vector{Float64})
        a, b, c = diagM(A)
        return(tridiagsolve(a, b, c, f))
end

function solution_a()    
	A = [8.0 9.0 4.0 -1.0; 0.0 4.0 1.0 0.0; 0.0 0.0 -1.0 6.0; 0.0 0.0 0.0 11.0]
	b = [9.0, 3.0, -1.0, 2.0]
	x = backwardsub(A, b)
	return (x, b - A * x)
end

function solution_b()
	A = Tridiagonal([-2 1 0 0 0; 1 -2 1 0 0; 0 1 -2 1 0; 0 0 1 -2 1; 0 0 0 1 -2])
	b = [1.0, 1.0, 1.0, 1.0, 1.0]
	x = tridiagsolve(A, b)
	return (x, b - A * x)
end

"Нестабильное LU-разложение квадратной матрицы `A`. Возвращает `L`, `U`."
function lufact(A::AbstractMatrix)
    n = size(A, 1)
    L = diagm(0 => ones(n))
    U = zeros(n, n)
    Aₖ = float(copy(A))	    
    for k in 1:n-1
   	 U[k, :] .= Aₖ[k, :]
	 L[:, k] .= Aₖ[:, k] ./ U[k, k]
	 Aₖ .-= L[:, k] * U[k, :]'
    end
    U[n, n] = Aₖ[n, n]
    return LowerTriangular(L), UpperTriangular(U)
end

function solution_c()
	A = [1 8 -3 9; 0 4 10 -2; 8 2 -5 1; 3 1 6 12]
	b = [3.0, 6.0, 1.0, 4.0] 
        #! Ок, оно отработало, но в общем случае всё-такие PLU надо использовать.
	L, U = lufact(A)
	z = forwardsub(L, b)
	x = backwardsub(U, z)
	return x, b - A * x
end


println(solution_a())
println(solution_b())
println(solution_c())

