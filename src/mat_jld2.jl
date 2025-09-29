"""
convert_mat_to_jld2 recibe dos rutas: la primera es la ruta al .mat a transformar; la segunda, una ruta que aun no existe y será la salida del archivo .jld2
"""

function convert_mat_to_jld2(ruta_mat::String, ruta_jld2::String)

    matfile = matopen(ruta_mat)
        variables = names(matfile)

        # Dict{ClaveTipo, ValorTipo}(Contenido del Dict)
        contenido = Dict{String, Any}()
        for var in variables
            # string(var): convierte en string a "var"
            contenido[string(var)] = read(matfile, var)
        end
    close(matfile)

    @save ruta_jld2 contenido #los ... hace que no se guarde una variable "contenido" y que se guarde en su lugar lo que hay dentro de "contenido".
    
    println("$(length(variables)) variable(s) guardadas en $ruta_jld2")
end