using TFG




using JLD2
@load "data/tbm/ns146raw_eo0/7x7_eo0.jld2" contenido
println(keys(contenido))
contenido["Mec_new"]
println(contenido["g"])
println(contenido["xc"])

39609+147
39714+42