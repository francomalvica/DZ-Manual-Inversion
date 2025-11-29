using DelimitedFiles

function main()
    # Read RA.dat file
    data = readdlm("RAo.dat")
    ab = data[:, 1]
    RAo = data[:, 2]
    i = length(ab)

    # Initialization
    RV = copy(RAo)
    z = copy(ab)
    x = zeros(Float64, i)
    x[1] = z[1]
    for j in 2:i
        x[j] = z[j] - z[j - 1]
    end

    min_thickness = 0.1    # ← Adjust the minimum desired thickness here (e.g. 0.1 or 1.0)

    aux2 = 1000.0
    rms  = 0.0

    # First loop: adjust thicknesses
    while true
        RAc = TRS(ab, x, RV)
        rms = sqrt(sum(((RAo[j] - RAc[j]) / RAo[j])^2 for j in 1:i) / i) * 100

        if rms < aux2
            aux2 = rms
            # Reduce layers 2..i, but never below min_thickness
            x[2:end] .= clamp.(0.9 .* x[2:end], min_thickness, Inf)

            # Rebuild cumulative depth array z
            z[1] = x[1]
            for j in 2:i
                z[j] = z[j-1] + x[j]
            end
        else
            break
        end
    end

    aux2 = 2 * rms
    tol = 0.2

    # Second loop: adjust resistivities
    while true
        RAc = TRS(ab, x, RV)
        sum_sq = sum(((RAo[j] - RAc[j]) / RAo[j])^2 for j in 1:i)
        rms = sqrt(sum_sq / i) * 100

        if (rms < aux2) && (rms > tol)
            aux2 = rms
            for j in 1:i
                RV[j] *= RAo[j] / RAc[j]
            end
        else
            break
        end
    end

    RAc = TRS(ab, x, RV)

    # Write results
    open("RVz.dat", "w") do f
        write(f, "Depth  RV\n")
        write(f, "0.1 $(RV[1])\n")
        for j in 1:(i - 1)
            write(f, "$(z[j]) $(RV[j])\n")
            write(f, "$(z[j]) $(RV[j + 1])\n")
        end
        write(f, "$(z[i]) $(RV[i])\n")
        write(f, "999.0 $(RV[i])\n\n")
    end

    open("RAo_RAc.dat", "w") do f    
        write(f, "AB   RAo  RAc\n")
        for j in 1:i
            write(f, "$(ab[j]) $(RAo[j]) $(RAc[j])\n")
        end
    end

    println("\nDONE!")
    println("--------------------------------------------------")
    println(" The following fitting error was obtained between  ")
    println(" the observed curve and the calculated curve:")
    println("--------------------------------------------------\n")
    println("RMS =", round(rms, digits=4) , "%")
    println("\nGenerated file: RVz.txt")
end


# TRS function (resistivity transform)
function TRS(y, e, p)
    h = length(y)
    RAc = zeros(Float64, h)

    a = [0.003042, -0.001198, 0.001284, 0.02350, 0.08688, 0.2374, 0.6194,
         1.1817, 0.4248, -3.4507, 2.7044, -1.1324, 0.3930, -0.1436,
         0.05812, -0.02521, 0.01125, -0.004978, 0.002072, -0.000318]
    s = -0.13069
    npd = 6
    npf = 20
    dx = log(10.0) / npd
    npfm1 = npf + 1

    for k in 1:h
        x = [exp(log(y[k]) + s + (j - 15) * dx) for j in 1:npf]
        T = zeros(Float64, npf)

        # Last layer
        for j in 1:npf
            L = (p[h] - p[h - 1]) / (p[h] + p[h - 1])
            aux1 = L * exp(-2 / x[j] * e[h - 1])
            M = (1 + aux1) / (1 - aux1)

            # Inner layers
            for i in 2:(h - 1)
                aux2 = p[h - i + 1] * M
                L = (aux2 - p[h - i]) / (aux2 + p[h - i])
                aux3 = L * exp(-2 / x[j] * e[h - i])
                M = (1 + aux3) / (1 - aux3)
            end

            T[j] = p[1] * M
        end

        RAc[k] = sum(a[npfm1 - h1] * T[h1] for h1 in 1:npf)
    end

    return RAc
end

# Run the program
main()
