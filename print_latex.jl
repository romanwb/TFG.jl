using Latexify, Format
Latexify.set_default(; starred=true)

what = latexify((ξ); fmt=x->format(round(x, sigdigits=2), autoscale=:metric))
println(what)            # lo copias del REPL
