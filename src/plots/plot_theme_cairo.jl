using CairoMakie

function use_formal_theme!(; kwargs...)
    theme = formal_theme(; kwargs...)
    return set_theme!(; theme...) 
end

function formal_theme(; fullscreen = false)
    f_padding = fullscreen ? (0, 0, 0, 0) : (10, 40, 10, 20)

    Theme(
        Figure = (resolution = DEFAULT_SIZE,),  # tamaño fijo

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
        palette = (
            color = COLORSCHEME,
            patchcolor = COLORSCHEME,
            marker = MARKERSLIST,
            strokecolor = COLORSCHEME
        ),

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
            xautolimitmargin = (0.05f0, 0.05f0),
            yautolimitmargin = (0.05f0, 0.05f0),

            # Grid & ticks settings
            xgridcolor = RGBAf(0, 0, 0, 0.1),
            xgridstyle = nothing,
            xgridwidth = 1.0,
            xgridvisible = GRID,
            ygridcolor = RGBAf(0, 0, 0, 0.1),
            ygridstyle = nothing,
            ygridwidth = 1.0,
            ygridvisible = GRID,

            xminorgridcolor = RGBAf(0, 0, 0, 0.15),
            xminorgridstyle = nothing,
            xminorgridvisible = MINORGRID,
            xminorgridwidth = 0.2,
            yminorgridcolor = RGBAf(0, 0, 0, 0.15),
            yminorgridstyle = nothing,
            yminorgridvisible = MINORGRID,
            yminorgridwidth = 0.2,

            # Ticks mayores/minores fijados
            xticks = LinearTicks(7),   # o pon 0:5:30 si quieres fijo
            yticks = LinearTicks(5),
            xminorticks = IntervalsBetween(1),
            yminorticks = IntervalsBetween(1),

            xminortickalign = 0.0,
            xminortickcolor = :black,
            xminorticksize = -5.0,
            xminorticksvisible = true,
            xminortickwidth = 0.8,
            yminortickalign = 0.0,
            yminortickcolor = :black,
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

        Legend = (
            backgroundcolor = :white,
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

        Colorbar = (
            bottomspinevisible = true,
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
