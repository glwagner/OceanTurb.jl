using Printf, StaticArrays, BenchmarkTools

Ns = (100,)
nts = (1,)

function many_inplace_multiplies!(a, b, n)
    for j = 1:n
        for i in eachindex(a)
            @inbounds a[i] *= b[i]
            @inbounds b[i] *= a[i]
        end
    end
    return nothing
end


function many_multiplies!(a, b, c, n)
    for j = 1:n
        for i in eachindex(a)
            @inbounds a[i] = b[i] * c[i]
            @inbounds b[i] = a[i] * c[i]
        end
    end
    return nothing
end

function many_adds!(a, b, c, n)
    for j = 1:n
        for i in eachindex(a)
            @inbounds a[i] = b[i] + c[i]
            @inbounds b[i] = a[i] + c[i]
        end
    end
    return nothing
end

for N = (1, 10, 100, 1000, 10000)
    nt = 100

    a = rand(N)
    b = rand(N)
    c = rand(N)

    println("N = $N")

    @printf "\nMany inplace multiplies:"
    t = @belapsed many_inplace_multiplies!($a, $b, $nt)
    @printf "\nt/N = %.1e s" t/N

    @printf "\nMany multiplies:"
    t = @belapsed many_multiplies!($a, $b, $c, $nt)
    @printf "\nt/N = %.1e s" t/N

    @printf "\nMany adds:"
    t = @belapsed many_adds!($a, $b, $c, $nt)
    @printf "\nt/N = %.1e s" t/N

    println("")
    println("")
end
