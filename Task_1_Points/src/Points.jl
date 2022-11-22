module Points
using LinearAlgebra

export Point, neighbors, Circle, Square, center  # публичные имена

struct Point #Точка на декартовой плоскости
    x
    y
end

struct Square{T <: Point}
    o::T   
    side
end
struct Circle{T <: Point} 
    o::T
    radius
end

LinearAlgebra.dot(p1::Point, p2::Point) = p1.x * p2.x + p1.y * p2.y
LinearAlgebra.norm(p::Point) = sqrt(p.x^2 + p.y^2)
Base.:+(p1::Point, p2::Point) = Point(p1.x + p2.x, p1.y + p2.y)
Base.:-(p1::Point, p2::Point) = Point(p1.x - p2.x, p1.y - p2.y)
Base.:-(p::Point) = Point(-p.x, -p.y)
Base.:*(α::Number, p::Point) = Point(α * p.x, α * p.y)
Base.:*(p::Point, α::Number) = α * p
Base.:/(p::Point, a::Number) = Point(p.x/a, p.y/a)
Base.:in(p::Point, c::Circle) = if (((p.x-c.o.x)^2 + (p.y-c.o.y)^2) <= c.radius^2)return true else return false end 
Base.:in(p::Point, s::Square) = if( s.side/2>= (s.o.x - p.x)  &&   s.side/2 >=( p.x - s.o.x)  && s.side/2 >=( s.o.y - p.y)  &&  s.side/2 >= (p.y - s.o.y)) return true else return false end 

center(points) = sum(points)/length(points)	 
"""
    center(points) -> Point
Центр "масс" точек.
"""

function center(points, c::Circle)
	cm = Point(0, 0)
	len_new = 0
	len = length(points)
	for i in 1:len
		if points[i] in c
			cm = cm + points[i]
			len_new = len_new + 1
    			end
	end
	return cm/len_new
end
	
function center(points, s::Square) 
	cm = Point(0, 0)   
	len_new = 0 
	len = length(points)
	for i in 1:len
		if points[i] in s     
			cm = cm + points[i]    
			len_new = len_new + 1
			end    
	end	
	return cm/len_new
end

function neighbors(points, origin::Point, k::Int64)
	if k>=length(points)
		arr = Vector{Point}([])
		for i in 1:length(points)
			if points[i] != origin 
				
				push!(arr, points[i])
				
			end
		end
		return arr
	end
	if k<=0 
		return [] end
	if k>0
		len = length(points)
		dist = Vector{Float64}([])
		for i in 1:len
			if points[i] != origin
				d = (origin.x - points[i].x)^2 + (origin.y - points[i].y)^2
				if d != 0
					push!(dist, ((origin.x - points[i].x)^2 + (origin.y - points[i].y)^2))
				end
			end
		end
		dist = sort(dist)
		dist = dist[1:k]
		arr = Vector{Point}([])
		len = length(dist)
		count=0
		for i in 1:length(points)
			if points[i] != origin
				if (((origin.x - points[i].x)^2 + (origin.y - points[i].y)^2) in dist)
					push!(arr, points[i])
					count+=1
					if count==k break end
				end
			end
		end
		return arr
	end
end


"""
    center(points, area) -> Point
Центр масс точек `points`, принадлежащих области `area` (`Circle` или `Square`).
"""

end # module
