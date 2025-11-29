#=
    graficar_CRV(rhos::Vector{Float64}, E::Vector{Float64}, zmax::Float64)

    Funcion para graficar un modelo de resistividades verdaderas en profundidad (CRV).
    La curva generada representa una estructura en capas, cada una con una resistividad uniforme.

    Parametros:
    - rhos: Vector de resistividades [Ohm·m] para cada capa. Debe tener longitud n.
    - E: Vector de espesores [m] para las primeras n-1 capas. La ultima capa se considera de espesor infinito.
    - zmax: Profundidad maxima [m] para visualizar la ultima capa extendida.
    - zmin: Profundidad minima [m] para visualizar la ultima capa extendida.

    Restriccion:
    - length(E) == length(rhos) - 1
    - La funcion pide valores REALES, entonces si son ENTEROS los escribo N.0, por ejemplo 1.0, 2.0, etc.

    Uso: 
    - Para usar esta funcion la tenemos que llamar en el codigo como "include("funcion_CRV.jl")"
    - Luego simplemente colocar el nombre de la funcion seguido de los parametros que se desean graficar
    - "graficar_CRV(rhos, E, zmax)"
=#

using PyPlot

function graficar_CRV(rhos::Vector{Float64}, E::Vector{Float64}, zmin::Float64, zmax::Float64)
    # Validacion: la cantidad de espesores (E) debe ser una menos que la de resistividades
    if length(rhos) != length(E) + 1
        error("El numero de espesores debe ser igual al numero de resistividades menos uno.")
    end

    # Configuracion del grafico
    figure(figsize=(10, 6))                          # tamano del grafico
    xscale("log")                                    # escala logaritmica en eje x
    yscale("log")                                    # escala logaritmica en eje y
    xlabel("Profundidad [m] y AB/2 [m]")             # etiqueta eje x
    ylabel("Resistividad [Ohm·m]")                   # etiqueta eje y
    title("Resistividad")                            # titulo del grafico
    grid(true, which="both", linestyle="--", linewidth=0.5)  # grilla
    xlim(zmin, zmax)                                 # limites del eje x 
    # Inicializacion de vectores de profundidad y resistividad
    z = [0.0]                  # inicio en superficie
    rho = [rhos[1]]            # resistividad inicial correspondiente

    # Construccion de la curva por capas
    for i in 1:length(E)
        z_fin = z[end] + E[i]                       # profundidad final de la capa actual
        append!(z, [z_fin, z_fin])                  # duplicacion para crear segmento horizontal
        append!(rho, [rhos[i], rhos[i+1]])          # resistividad constante hasta la siguiente capa
    end

    append!(z, zmax)                                # extension de la ultima capa hasta zmax
    append!(rho, rhos[end])                         # resistividad constante en la capa infinita

    # Grafico final
    plot(z, rho, "--", label="CRV", color=:green)   # traza la curva en forma escalonada
    legend()
end