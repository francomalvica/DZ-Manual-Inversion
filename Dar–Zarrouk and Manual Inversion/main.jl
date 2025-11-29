# IMPORTANT VARIABLES TO CHECK BEFORE RUNNING THE CODE:
# zmin and zmax
# skipstart = "number of header rows"

println("╔════════════════════════════════════════════════════════════════════╗")
println("║                 Modeling of True Resistivities                     ║")
println("╚════════════════════════════════════════════════════════════════════╝")
println()
println("This program performs the following steps:")
println("  1. Reads the file 'RAo.dat' containing observed apparent resistivities.")
println("  2. Computes a layered model using the Zohdy algorithm.")
println("  3. Allows reducing the number of layers using:")
println("       a) Automatic Dar–Zarrouk method")
println("       b) Manual layer selection")
println("IMPORTANT")
println("   - The file 'RAo.dat' must be located in the program directory.")
println("   - If 'RAo.dat' contains a header, modify the `skipstart` variable")
println("     in the code to specify how many rows should be skipped.")
println("   - If you wish to modify the minimum and maximum depth values of the plot")
println("     (zmin and zmax), you must do so manually in the source code.")
println("   - It is recommended to keep the RMS error below 3% between the observed and modeled")
println("     apparent resistivity curves to ensure a good-quality fit.")
println()
print("Press ENTER to begin... ")
readline()



using PyPlot
using DelimitedFiles

zmax = 1000.0   # maximum depth
zmin = 1        # minimum depth/minimum thickness

# Open the figure to plot the curves
figure(figsize=(10, 6))                                 # figure size
xscale("log")                                           # logarithmic scale on x-axis
yscale("log")                                           # logarithmic scale on y-axis
xlabel("Depth [m] and AB/2 [m]")                        # x-axis label
ylabel("Resistivity [Ohm·m]")                           # y-axis label
title("Resistivity")                                    # plot title
grid(true, which="both", linestyle="--", linewidth=0.5) # grid
xlim(zmin, zmax)                                        # x-axis limits


# ARC (observed) -------------------------------------------------------------------
RAo = readdlm("RAo.dat", skipstart=0)
z_o = RAo[:,1]
CRAo = RAo[:,2]
plot(z_o, CRAo, "-*", label="Observed ARC", color=:red, lw=2)
legend()
# ---------------------------------------------------------------------------------

# TRC (Zohdy algorithm) – INVERSE PROBLEM (TRC_model from ARC observed) -----------
include("algoritmo_zohdy.jl")

RVz = readdlm("RVz.dat", skipstart=1)
thicknesses_zohdy = RVz[:,1]
TRC_zohdy = RVz[:,2]
plot(thicknesses_zohdy, TRC_zohdy, "--", label="TRC by Zohdy", color=:red)
legend()
# ---------------------------------------------------------------------------------

####################################################################################
####################################################################################
################################ DAR–ZARROUK #######################################
####################################################################################
####################################################################################
include("funcion_dar.jl")     # interactive function

# Preprocessing:
e_aux = diff(thicknesses_zohdy)   # compute differences because thicknesses are duplicated
ind = findall(!=(0), e_aux)       # filter duplicates
e = e_aux[ind]                    # define correct thickness sequence
rv = TRC_zohdy[ind]               # define correct resistivities
ARC_original = CRAo               # save the original ARC

# Call the interactive mode from funcion_dar.jl
e, rv, rms_final = interactivo_DZ!(e, rv, z_o, ARC_original)
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

# ---------------------------------------------------------------------------------
# Save the model
println("\nDo you want to save the final model to a file? (y/n)")
answer = lowercase(readline())

if answer == "y"
    open("final_model.txt", "w") do io
        println(io, "=== Final model: $(length(rv)) layers ===")
        println(io, rpad("Layer",6), rpad("Thickness (m)",14),
                rpad("Resistivity (Ω·m)",22), "Depth (m)")
        println(io, "────────────────────────────────────────────────────────────")

        z = cumsum(e)
        k = length(e)
        for i in 1:k
            @printf(io, "%5d     %10.2f        %14.2f         %10.2f\n",
                    i, e[i], rv[i], z[i])
        end
        @printf(io, "%5d           ---         %14.2f                ∞\n",
                length(rv), rv[end])

        println(io, "\n📉 RMS relative to the original curve: $(round(rms_final, digits=2)) %")
    end
    println("💾 Final model saved as 'final_model.txt'")
else
    println("❎ No file was saved.")
end
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# Save the plot
println("\nDo you want to save the plot as an image? (y/n):")
answer = readline()

if lowercase(answer) == "y"
    println("Enter the filename (without extension):")
    filename = readline()
    savefig(filename * ".png")
    println("💾 Plot saved as $(filename).png")
else
    println("❎ Plot was not saved.")
end
# ---------------------------------------------------------------------------------

PyPlot.close("all")  # Close all open figures

println()
println("╔════════════════════════════════════════════════════════════════════╗")
println("║                          Program finished                          ║")
println("╚════════════════════════════════════════════════════════════════════╝")
println()
println("Franco Malvica - 2025")
println("Electrical Prospecting Methods")
println("Faculty of Astronomical and Geophysical Sciences - UNLP")
println()
