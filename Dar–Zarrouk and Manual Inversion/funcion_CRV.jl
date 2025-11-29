#=
    graficar_CRV(rhos::Vector{Float64}, E::Vector{Float64}, zmax::Float64)

    Function to plot a true resistivity model in depth (TRM/CRV).
    The generated curve represents a layered structure, each one with a uniform resistivity.

    Parameters:
    - rhos: Vector of resistivities [Ohm·m] for each layer. Must have length n.
    - E: Vector of thicknesses [m] for the first n–1 layers. The last layer is considered to have infinite thickness.

    Restriction:
    - length(E) == length(rhos) - 1
    - The function requires REAL values, so if they are INTEGERS, write them as N.0, for example 1.0, 2.0, etc.

    Usage:
    - To use this function, we must call it in the code as: include("funcion_CRV.jl")
    - Then simply call the function followed by the parameters you want to plot:
    - "graficar_CRV(rhos, E, zmax)"
=#


function graficar_CRV(rhos::Vector{Float64}, E::Vector{Float64})
    # Validation: the number of thicknesses (E) must be one less than the number of resistivities
    if length(rhos) != length(E) + 1
        error("The number of thicknesses must be equal to the number of resistivities minus one.")
    end

    # PyPlot.cla()
    # Initialization of depth and resistivity vectors
    z = [0.0]                  # start at the surface
    rho = [rhos[1]]            # initial corresponding resistivity

    # Construction of the stepwise curve
    for i in 1:length(E)
        z_fin = z[end] + E[i]                       # final depth of the current layer
        append!(z, [z_fin, z_fin])                  # duplication to create a horizontal segment
        append!(rho, [rhos[i], rhos[i+1]])          # constant resistivity until the next layer
    end

    append!(z, zmax)                                # extend the last layer down to zmax
    append!(rho, rhos[end])                         # constant resistivity in the infinite layer

    # Final plot
    linea = plot(z, rho, "--", label="Reduced TRM", color=:blue)[1]    # plots the stepwise curve
    legend()
    return linea  # <--- return the object so it can be removed later


end
