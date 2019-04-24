using Pkg; Pkg.activate("..")

using OceanTurb, Printf

for stepper in (:ForwardEuler, :BackwardEuler)
    model = Diffusion.Model(N=100, stepper=stepper)

    iterate!(model, 0.1)

    # Measure memory allocation
    @printf "\nTesting Diffusion with %s...\n" stepper
    for nt = (10, 100, 1000, 10000)
        @printf "nt: % 6d" nt
        @time iterate!(model, 0.1, nt)
    end
end

for stepper in (:ForwardEuler, :BackwardEuler)
    model = KPP.Model(N=100, stepper=stepper)

    iterate!(model, 0.1)

    # Measure memory allocation
    @printf "\nTesting KPP with %s...\n" stepper
    for nt = (10, 100, 1000, 10000)
        @printf "nt: % 6d" nt
        @time iterate!(model, 0.1, nt)
    end
end
