x = [1, 1, 1, 2, 2, 2]
n = length(x) 
sum(x) 
μ = sum(x) / n 
sum_2 = sum((xi - μ)^2 for xi in x)/n 
sqrt(sum_2) 


function std_dev(x)
    n = length(x) 
    μ = sum(x) / n
    sum_2 = sum((xi - μ)^2 for xi in x)/n 
    return sqrt(sum_2)
end

println(std_dev(x)) 