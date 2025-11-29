using PyPlot
include("filtro_oneill_final.jl")


# PLOT -------------------------------------------------------------------------

# ARV (with the input data) ----------------------------------------------------
include("funcion_CRV.jl") # Add the function I made to plot the ARV curves
graficar_CRV(p, e, ab_min, ab_max)
# ------------------------------------------------------------------------------
# AAR (with the computed data) -------------------------------------------------
RA = readdlm("RA.dat")
AB2 = RA[:,1] 
AAR = RA[:,2]
plot(AB2, AAR, "-o", label="AAR", color=:blue)
legend()
# ------------------------------------------------------------------------------


# Save the plot
println("\nDo you want to save the plot as an image? (y/n):")
respuesta = readline()

if lowercase(respuesta) == "y"
    println("Enter the file name (without extension):")
    nombre_archivo = readline()
    savefig(nombre_archivo * ".png")
    println("Plot saved as $(nombre_archivo).png")
else
    println("The plot was not saved.")
end
# -------------------------------------------------------------------------------
