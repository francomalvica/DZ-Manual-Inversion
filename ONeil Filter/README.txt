To compute apparent resistivity curves from real resistivity values
(i.e., our model), you must run the code "correr_filtro.jl" in Julia. That is:

julia> include("correr_filtro.jl")


This code uses the Oneill filter (from the file "filtro_oneill.jl") and I also added a call to
"funcion_CRV.jl", which contains the function "graficar_CRV", used to plot the curves
(or more precisely, the straight lines) corresponding to the values entered when executing the
Oneill filter.

Additionally, I included a routine to plot the apparent resistivity curves (ARCs) from the file
that is generated (RA.dat), overlaying them on the plot produced by the graficar_CRV function.
This way, when running this Julia program, I obtain the plot with both the RC (true resistivity)
and ARC (apparent resistivity) curves, as well as the file containing the apparent resistivity
calculated for each point.

⚠️ IMPORTANT!
For filtro_oneill_final.jl to work, funcion_CRV.jl must be located in the same directory.