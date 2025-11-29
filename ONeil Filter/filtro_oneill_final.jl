using DelimitedFiles

# Display program information
println("_____________________________________________________")
println("   The program computes Apparent Resistivities       ")
println("    using O'Neill's filter and Sunde's algorithm     ")
println("           to calculate the TDR (FTR).               ")
println("_____________________________________________________\n")

# O'Neill Filter (coefficients)
a = [
    0.003042, -0.001198, 0.01284, 0.02350, 0.08688, 0.2374,
    0.6194, 1.1817, 0.4248, -3.4507, 2.7044, -1.1324,
    0.3930, -0.14362, 0.05812, -0.02521, 0.01125, -0.004978,
    0.002072, -0.000318
]

# Filter parameters
s = -0.13069
npd = 6
npf = 20
dx = log(10.0) / npd
npfm1 = npf + 1

println("-----------------------------------------------------")
println(" How would you like to input the AB/2 points for the calculation?")
println(" 1 - Read them from a file named 'Abs.dat'")
println("     (it must be in the same folder as this script)")
println(" 2 - Generate them automatically (logarithmic distribution)")
println("-----------------------------------------------------")
println("Enter 1 or 2 and press Enter:")
modo_ingreso = parse(Int, readline())

if modo_ingreso == 1
    # Read from Abs.dat
        try
        global ab = readdlm("Abs.dat")[:, 1]
        global n = length(ab)
        println("Read $n points from 'Abs.dat'")
        global ab_min = ab[1]
        global ab_max = ab[n]
    catch err
        println("\n[ERROR] Could not read the file 'Abs.dat'.")
        println("Check that the file exists and contains a column of numbers.")
        println("Exiting program...")
        exit(1)
    end

elseif modo_ingreso == 2
    # Generate automatically
    println("Enter the number of points per decade:")
    puntos_por_decada = parse(Int, readline())

    println("Enter the minimum AB/2 value (for example, 1.0):")
    ab_min = parse(Float64, readline())

    println("Enter the maximum AB/2 value (for example, 100.0):")
    ab_max = parse(Float64, readline())

                    # -------------------------------------------------------
                    # CALCULATION OF LOG-EQUISPACED POINTS
                    # -------------------------------------------------------
                    # Decade calculation
                    log10_min = log10(ab_min)
                    log10_max = log10(ab_max)
                    num_decadas = round(Int, log10_max - log10_min)

                    # Result vector
                    ab = Float64[]

                    for d in 0:num_decadas-1
                        # Log-equispaced range for one decade (includes endpoints)
                        puntos = 10 .^ LinRange(0, 1, puntos_por_decada)

                        # Scale the decade
                        puntos_d = ab_min * 10.0^d .* puntos

                        if d == 0
                            append!(ab, puntos_d)
                        else
                            # Skip first point (already included)
                            append!(ab, puntos_d[2:end])
                        end
                    end

                    # Safety filter (in case of rounding)
                    ab = filter(x -> x >= ab_min && x <= ab_max, ab)
                    global n = length(ab)

else
    error("Invalid option. You must enter 1 or 2.")
end

# Request geoelectric section parameters
println("+++++++++++++++++++++++++++++++++++++++++++++++++++++")
println(" You must now enter the parameters of the            ")
println(" geoelectric section:                                ")
println("   - Number of layers (max = 10)                     ")
println("   - Resistivities of each layer                     ")
println("   - Thicknesses of each layer except the last       ")
println("+++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

println("Enter number of layers:")
nc = parse(Int, readline())

p = Float64[]
e = Float64[]

for i in 1:nc
    println("Layer $i - Resistivity [ohm.m]:")
    push!(p, parse(Float64, readline()))
end

for i in 1:(nc - 1)
    println("Layer $i - Thickness [m]:")
    push!(e, parse(Float64, readline()))
end

# Function that computes the resistivity transform (Sunde algorithm)
function trs(x, nc, p, e)
    h = length(x)
    T = zeros(Float64, h)
    for j in 1:h
        L = (p[nc] - p[nc - 1]) / (p[nc] + p[nc - 1])
        aux1 = L * exp(-2 / x[j] * e[nc - 1])
        M = (1 + aux1) / (1 - aux1)
        if nc - 1 ≥ 2
            for i in 2:(nc - 1)
                aux2 = p[nc - i + 1] * M
                L = (aux2 - p[nc - i]) / (aux2 + p[nc - i])
                aux3 = L * exp(-2 / x[j] * e[nc - i])
                M = (1 + aux3) / (1 - aux3)
            end
        end
        T[j] = p[1] * M
    end
    return T
end

# Open output file and compute apparent resistivities
open("RA.dat", "w") do io
    for k in 1:n
        # Generate log-spaced points centered on log(ab[k])
        x = [exp(log(ab[k]) + s + (j - 15) * dx) for j in 1:npf]

        # Compute transform
        T = trs(x, nc, p, e)

        # Convolution with O’Neill filter
        suma = sum(a[npfm1 - h] * T[h] for h in 1:npf)

        # Save result
        println(io, ab[k], " ", suma)
    end
end

println("Generated RA.dat file with the apparent resistivities.")
