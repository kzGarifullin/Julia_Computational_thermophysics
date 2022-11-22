#! 1. return -- оператор, а не функция, не нужно ставить у него скобки
#! 2. Слишком жёсткие ограничения на типы аргументов. Если матрица не из Int-ов, методы не сработают?

using LinearAlgebra


function forwardsub(L::AbstractMatrix, b::Vector{Float64})         

	x = similar(b, Float64) # память заранее
	#x = zeros(length(b))
	#x = Vector{Float64}([])
	n,m = size(L)

	#! 1. Такие if можно писать лаконичнее, через short-circuit операторы && и ||.
	#! 2. Почему return? Очевидно, что есть ошибка. Надо сообщить пользователю, чтобы он искал её тут, а не где-то дальше по коду!
	#!
	if(n != m || n != length(b)) && error("ВСЁ ПЛОХО, ПОЛЬЗОВАТЕЛЬ, ОШИБКА ИМЕННО ЗДЕСЬ!")
	end
	
	#if(n != m || n != length(b))
	#	return x
	#end

	#! Для проверки можно не выделять объект, полистайте документацию к LinearAlgebra.
	#B = LowerTriangular(L)
	#! Аналогично if-у выше.
	
	if(L != LowerTriangular(L)) && error("Ошибка! Матрица не треугольная")
        end
	

	#if(B != L)
	#	return ("Error! Matrix must be LowerTriangular")
	#end

	#! Эта операция дороже, чем всё остальное.
	if(det(L)==0) &&error("Ошибка! вырождена")
	end
	#if(det(L) == 0)
	#	return("Error! Singular matrix")
	#end
	k=0
	for i in 1:length(b)
		sum = 0 
		for j in 1:k
			sum = sum + L[i, j]*x[j]
		end
		k+=1
		#push!(x, (b[i] - sum)/L[i, i])
		x[i] = (b[i] - sum)/L[i, i] 	
	end
	return x
end



#!
function backwardsub(U::AbstractMatrix, b::Vector{Float64})
        
	x = zeros(length(b))
        n,m = size(U)
        #if(n != m || n != length(b))
        #        return("No solutions!")
        #end
        
	if(n != m || n != length(b)) && error("ВСЁ ПЛОХО, ПОЛЬЗОВАТЕЛЬ, ОШИБКА ИМЕННО ЗДЕСЬ!!!!!!!!!!!!!")        end
	
	#B = UpperTriangular(U)
        #if(B != U)
        #        return ("Error! Matrix must be UpperTriangular")
        #end
        
	if(U != UpperTriangular(U)) && error("Ошибка! Матрица не треугольная!!!!!!!!!!!!!!!!!!!!!")        end
	#if(det(U) == 0)
        #        return("Error! Singular matrix")
        #end
        if(det(U)==0) &&error("Ошибка! вырождена!!!!!!!!!!!!")        end

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
        return x
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

function solution_c()
	A = [1 8 -3 9; 0 4 10 -2; 8 2 -5 1; 3 1 6 12]
	b = [3.0, 6.0, 1.0, 4.0] 
        #! Ок, оно отработало, но в общем случае всё-такие PLU надо использовать.
    F = lu(A)
	z = forwardsub(F.L, b[F.p])
	x = backwardsub(F.U, z)
	return x, b - A * x
end
println(solution_a())
println(solution_b())
println(solution_c())
L = [4.0 0.0 0.0 0.0; 1.0 2.0 0.0 0.0; 1.0 2.0 1.0 0.0; 1.0 3.0 2.0 1.0]
L_err = [4 0 0 0; 1 1 0 0; 1 2 1 0]
L_err_nm = [4 0 0 0; 4 5 0 0] 
U = [1 0 0 1; 0 2 0 0; 0 0 1 1; 0 0 0 6]
U = [1.0 5.0 6.0 1.0; 0.0 2.0 7.0 8.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 6.0]
b = [4.0, 2.0, 3.0, 6.0]
L_err_vir = [4 0 0 0; 1 2 0 0; 1 3 2 1 ; 1 3 2 1]
L_err_nottriag = [4 0 0 0; 1 2 0 1; 1 2 1 0; 1 3 2 1]
U_err_nm = [1 0 0 1; 0 0 1 1; 0 0 0 6]
U_err = [1.0 5.0 6.0 1.0; 2.0 2.0 7.0 8.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 6.0]
println(forwardsub(L, b))
println(backwardsub(U, b))
println(backwardsub(U_err, b))
