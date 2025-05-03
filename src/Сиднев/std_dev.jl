x = Float64[12, 23, 45, 33, 65, 54, 54]

function std_dev(x)
    n = length(x)
    mu = sum(x) / n
    sigma = (sum((x.-mu).^2) / n).^(1/2)
    return sigma
end

println(std_dev(x))


