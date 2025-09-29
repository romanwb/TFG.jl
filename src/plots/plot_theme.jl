using GLMakie

function use_formal_theme!(; kwargs...) # ; convierte argumentos en keyword-only. kwargs... agrupa las palabras clave pasadas a la funcion en un NamedTuple
    theme = formal_theme(; kwargs...)
    return set_theme!(; theme...) # function defined by Makie.jl
end

function formal_theme(; fullscreen = false)
    f_padding = fullscreen ? (0, 0, 0, 0) : (10, 40, 10, 20)

    Theme(
        Figure = (resolution = DEFAULT_SIZE,),

        # Font settings
        fonts = (
            regular = "CMU",
            bold = "CMU",
            italic = "CMU",
            monospace = "CMU",
            ),
        fontsize = FONTSIZE,
        textcolor = :black,

        # Padding settings
        figure_padding = f_padding,
        linewidth = LINEWIDTH,

        # Color settings
        palette = (color = COLORSCHEME,
            patchcolor = COLORSCHEME,
            marker = MARKERSLIST,
            strokecolor = COLORSCHEME),

        # Scatter settings
        markersize = MARKERSIZE,
        Scatter = (cycle = Cycle([:marker, :strokecolor], covary = true),
            strokewidth = 2, color = :white),
        ScatterLines = (cycle = Cycle([:marker, :color, :strokecolor], covary = true),
            markercolor = :white, strokewidth = 2),

        # Axis settings
        Axis = (
            titlesize = 24,

            # Label settings
            xlabel = L"x",
            xlabelpadding = 0.0,
            xlabelsize = 20,
            ylabel = L"y",
            ylabelpadding = 5.0,
            ylabelsize = 20,

            # Spine settings
            spinewidth = 1.5,
            topspinecolor = :black,
            rightspinecolor = :black,
            leftspinecolor = :black,
            bottomspinecolor = :black,

            # Margin settings
            xautolimitmargin = (0.05f0, 0.05f0), # Default
            yautolimitmargin = (0.05f0, 0.05f0), # Default

            # Grid & ticks settings
            xgridcolor = RGBf(0.85, 0.85, 0.85),
            xgridstyle = nothing, # Default
            xgridwidth = 1.0, # Default
            xgridvisible = GRID, # Default
            ygridcolor = RGBf(0.85, 0.85, 0.85),
            ygridstyle = nothing, # Default
            ygridwidth = 1.0, # Default
            ygridvisible = GRID, # Default
            
            xminorgridcolor = RGBf(0.95, 0.95, 0.95),
            xminorgridstyle = nothing,
            xminorgridvisible = MINORGRID,
            xminorgridwidth = 0.2,
            yminorgridcolor = RGBf(0.95, 0.95, 0.95),
            yminorgridstyle = nothing,
            yminorgridvisible = MINORGRID,
            yminorgridwidth = 0.2,

            xminortickalign = 0.0, # Ticks inside 
            xminortickcolor = :black,
            xminorticks = IntervalsBetween(2),
            xminorticksize = -5.0,
            xminorticksvisible = true,
            xminortickwidth = 0.8,
            yminortickalign = 0.0, # Ticks inside 
            yminortickcolor = :black,
            yminorticks = IntervalsBetween(2),
            yminorticksize = -5.0,
            yminorticksvisible = true,
            yminortickwidth = 0.8,

            xtickalign = 1,
            xtickcolor = :black,
            xticklabelpad = 10,
            xticksize = 10,
            xtickwidth = 1.0,
            xticksmirrored = true,
            ytickalign = 1,
            ytickcolor = :black,
            yticklabelpad = 10,
            yticksize = 10,
            ytickwidth = 1.0,
            yticksmirrored = true,
            ),
        Legend = (backgroundcolor = :white,
            framecolor = (:black, 1),
            framewidth = 1.5,
            margin = (25, 25, 25, 25),
            orientation = :vertical,
            padding = (10, 10, 8, 8),
            patchlabelgap = 15,
            patchsize = (70, 20),
            halign = :left,
            valign = :top,
            tellheight = false,
            tellwidth = false
            ),
        Colorbar = (bottomspinevisible = true,
            minorticksvisible = false,
            labelpadding = 25,
            labelsize = LABELSIZE,
            spinewidth = 1.75,
            tickcolor = :black,
            ticksize = 10,
            tickalign = 0,
            ticklabelpad = 7.5,
            tickwidth = 1.5,
            width = 30
            )
    )
end

use_formal_theme!()
