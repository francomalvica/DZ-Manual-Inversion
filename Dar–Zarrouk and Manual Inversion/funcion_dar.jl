using Printf
using Statistics
include("funcion_CRV.jl")

# -------------------------------------------------------------------------------------------
# FUNCTION TO MERGE LAYERS USING DAR ZARROUK ------------------------------------------------
# -------------------------------------------------------------------------------------------

function fusionar_capas_DZ!(e::Vector{Float64}, rv::Vector{Float64}, capas::Vector{Int})
    if length(capas) < 2
        println("⚠️ You must select at least two layers to merge.")
        return false
    end

    # Check that the layers are within the valid range (except the last infinite one)
    if maximum(capas) >= length(rv)
        println("⚠️ You cannot merge the last infinite layer.")
        return false
    end

    # Compute aggregated T and S
    T_total = sum(e[capas] .* rv[capas])
    S_total = sum(e[capas] ./ rv[capas])

    e_eq = sqrt(T_total * S_total)
    rho_eq = sqrt(T_total / S_total)

    # Remove the selected layers
    deleteat!(e, capas)
    deleteat!(rv, capas)

    # Insert new merged layer
    insert!(e, minimum(capas), e_eq)
    insert!(rv, minimum(capas), rho_eq)

    # Ensure rv has one more element than e
    if length(rv) == length(e)
        push!(rv, rv[end])
    elseif length(rv) > length(e) + 1
        pop!(rv)
    end

    println("✔️ Layers $(capas) merged → new layer:")
    println("    - thickness = $(round(e_eq, digits=2)) m")
    println("    - resistivity = $(round(rho_eq, digits=2)) Ω·m")

    return true
end


# -------------------------------------------------------------------------------------------
# FUNCTION THAT CALLS THE PREVIOUS ONE FOR ITERATIVE LAYER MERGING --------------------------
# -------------------------------------------------------------------------------------------
function interactivo_DZ!(e::Vector{Float64}, rv::Vector{Float64}, z_c::Vector{Float64}, CRA_original::Vector{Float64}; umbral=3.0)
    linea_anterior_CRV = nothing
    linea_anterior_CRA = nothing
    ion()
    rms_final = 0.0  # initial value

    while true
        # Ensure rv has one more layer than e
        if length(rv) == length(e)
            push!(rv, rv[end])
        elseif length(rv) > length(e) + 1
            pop!(rv)
        end

        # Show model table
        println("\n=== Current model: $(length(rv)) layers ===")
        println(rpad("Layer",6), rpad("Thickness (m)",14), rpad("Resistivity (Ω·m)",22), "Depth (m)")
        println("────────────────────────────────────────────────────────────")
        z = cumsum(e)
        for i in 1:length(e)
            @printf("%5d     %10.2f        %14.2f         %10.2f\n", i, e[i], rv[i], z[i])
        end
        @printf("%5d           ---         %14.2f                ∞\n", length(rv), rv[end])

        println("\nSelect an option:")
        println("[1] Merge layers (Dar-Zarrouk)")
        println("[2] Modify a layer manually")
        println("[ENTER] to exit")
        opcion = readline()
        isempty(opcion) && break

        if opcion == "1"
            println("Enter layers to merge (e.g., 2 3 4):")
            input = readline()
            capas = sort(parse.(Int, split(input)))
            ok = fusionar_capas_DZ!(e, rv, capas)
            if !ok continue end

        elseif opcion == "2"
            println("Enter the layer number to modify (1-$(length(rv))):")
            idx = parse(Int, readline())
            if idx < 1 || idx > length(rv)
                println("⚠️ Layer index out of range.")
                continue
            end

            println("What do you want to modify?")
            println("[1] Thickness")
            println("[2] Resistivity")
            subop = readline()

            if subop == "1" && idx <= length(e)
                println("Enter the new thickness for layer $idx:")
                nuevo_e = parse(Float64, readline())

                if iszero(nuevo_e)
                    println("⚠️ Thickness entered is 0. Removing layer $idx.")

                    deleteat!(e, idx)
                    deleteat!(rv, idx)

                    # Fix dimensions if needed
                    if length(rv) == length(e)
                        push!(rv, rv[end])
                    end

                else
                    delta = nuevo_e - e[idx]

                    if idx <= length(e) - 1
                        if e[idx+1] - delta <= 0
                            println("⚠️ The new thickness completely consumes the next layer.")
                            println("→ Removing layer $(idx+1).")
                            deleteat!(e, idx+1)
                            deleteat!(rv, idx+1)
                        else
                            e[idx+1] -= delta
                        end
                    elseif idx == length(e)
                        e[idx] = nuevo_e
                    else
                        println("⚠️ Cannot modify this layer.")
                        continue
                    end

                    e[idx] = nuevo_e
                end

            elseif subop == "2"
                println("Enter the new resistivity for layer $idx:")
                nuevo_rho = parse(Float64, readline())

                fusion_siguiente = (idx < length(rv) && round(nuevo_rho, digits=2) == round(rv[idx+1], digits=2))
                fusion_anterior  = (idx > 1 && round(nuevo_rho, digits=2) == round(rv[idx-1], digits=2))

                rv[idx] = nuevo_rho

                if fusion_siguiente
                    println("⚠️ Matches next layer. Merging layers $idx and $(idx+1).")
                    e[idx] += e[idx+1]
                    deleteat!(e, idx+1)
                    deleteat!(rv, idx+1)
                    println("    → New layer $idx: thickness = $(round(e[idx],digits=2)) m, ρ = $(round(rv[idx],digits=2)) Ω·m")

                elseif fusion_anterior
                    println("⚠️ Matches previous layer. Merging layers $(idx-1) and $idx.")
                    e[idx-1] += e[idx]
                    deleteat!(e, idx)
                    deleteat!(rv, idx)
                    println("    → New layer $(idx-1): thickness = $(round(e[idx-1],digits=2)) m, ρ = $(round(rv[idx-1],digits=2)) Ω·m")

                else
                    println("✔️ Resistivity of layer $idx updated to $(round(rv[idx],digits=2)) Ω·m.")
                end

            else
                println("⚠️ Invalid option or infinite layer without thickness.")
                continue
            end
        else
            println("⚠️ Invalid option.")
            continue
        end

        # Compute and show reduced AAR + RMS
        CRA_reducida = TRS(z_c, e, rv)
        rms = sqrt(mean(((CRA_original .- CRA_reducida) ./ CRA_original).^2)) * 100
        rms_final = rms
        println("📉 RMS relative to original curve: $(round(rms, digits=2)) %")

        # Plot CRV and new AAR
        if linea_anterior_CRV !== nothing
            linea_anterior_CRV.remove()
        end
        if linea_anterior_CRA !== nothing
            linea_anterior_CRA.remove()
        end

        linea_anterior_CRV = graficar_CRV(rv, e)
        linea_anterior_CRA = plot(z_c, CRA_reducida, "-*", label="Reduced AAR", color=:blue, lw=1)[1]
        legend()
        PyPlot.draw()
        sleep(0.1)
    end

    return e, rv, rms_final

end




# -------------------------------------------------------------------------------------------
# RESISTIVITY TRANSFORM FUNCTION TO COMPUTE APPARENT RESISTIVITY ----------------------------
# -------------------------------------------------------------------------------------------
function TRS(y, e, p)
    npts = length(y)
    ncapas = length(p)
    nesp = length(e)
    RAc = zeros(npts)

    # O’Neill filter
    a = [0.003042, -0.001198, 0.001284, 0.02350, 0.08688, 0.2374, 0.6194,
         1.1817, 0.4248, -3.4507, 2.7044, -1.1324, 0.3930, -0.1436,
         0.05812, -0.02521, 0.01125, -0.004978, 0.002072, -0.000318]
    npf, s, npd = 20, -0.13069, 6
    dx = log(10.0)/npd
    npfm1 = npf + 1

    for k in 1:npts
        x = [exp(log(y[k]) + s + (j - 15) * dx) for j in 1:npf]
        T = zeros(npf)

        for j in 1:npf
            # Last layer
            L = (p[end] - p[end - 1]) / (p[end] + p[end - 1])
            aux1 = L * exp(-2 / x[j] * e[end])
            M = (1 + aux1) / (1 - aux1)

            for i in 2:(ncapas - 1)
                idx = ncapas - i
                if idx > nesp
                    continue
                end
                aux2 = p[idx + 1] * M
                L = (aux2 - p[idx]) / (aux2 + p[idx])
                aux3 = L * exp(-2 / x[j] * e[idx])
                M = (1 + aux3) / (1 - aux3)
            end

            T[j] = p[1] * M
        end

        RAc[k] = sum(a[npfm1 - h1] * T[h1] for h1 in 1:npf)
    end

    return RAc
end

