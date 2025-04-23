using CairoMakie

function save_figure_pdf(filename::String, fig::Figure)
    save(filename, fig)
end