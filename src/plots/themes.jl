using CairoMakie  # usa CairoMakie para permitir guardar en PDF
function apply_latex_style!()
    CairoMakie.activate!()

    set_theme!(Theme(
        font = "Latin Modern Roman",  # o "Latin Modern Roman" si lo prefieres
        fontsize = 18,
        linewidth = 2,
        Axis = (
            xlabelpadding = 10,
            ylabelpadding = 10,
            xticklabelsize = 18,
            yticklabelsize = 18,
            xticklabelfont = "Latin Modern Roman",
            yticklabelfont = "Latin Modern Roman",
            xticksize = 5,
            yticksize = 5,
            xtickwidth = 1,
            ytickwidth = 1,
            spinewidth = 0.8
        )
    ))
end
