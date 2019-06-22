using PyPlot

d = 0:0.01:1

function kp(d, A...)
    s = d * (1-d)
    r = 0.0
    for (i, a) in enumerate(A)
        r += a * (1-d)^(i-1)
    end
    return s * r
end

fig, axs = subplots()
plot(d, kp.(d, 0.0, 1.0), "-")
plot(d, kp.(d, 1.0, 0.0), "--")
plot(d, kp.(d, 0.1, 1.0), "-.")
plot(d, kp.(d, 1.0, -1.0), ":")
gcf()
