using LinearAlgebra
function forwardsub(L::Matrix{Int64}, b::Vector{Float64})         
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

function backwardsub(U::Matrix{Int64}, b::Vector{Float64})
        
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
L = [4 0 0 0; 1 2 0 0; 1 2 1 0; 1 3 2 1]
L_err = [4 0 0 0; 1 1 0 0; 1 2 1 0]
U = [1 0 0 1; 0 2 0 0; 0 0 1 1; 0 0 0 6]
U = [1 5 6 1; 0 2 7 8; 0 0 1 1; 0 0 0 6]
b = [4.0, 2.0, 3.0, 6.0]

println(forwardsub(L, b))
println(backwardsub(U, b))
