"Interseccion de circunferencia con circunferencias de diferente radio"

λ_vals = range(5, 20, length=10)  # Más valores de lambda para más circunferencias

circumference = Figure(resolution = (800, 800))
ax = Axis(circumference[1, 1], xlabel="x", ylabel="y", title="Intersección de Recta y Circunferencias",
        limits = (-5, 5, -5, 5))

x_range = range(-10, 10, length=100)  # Rango de x para la recta

# Inicializar listas para guardar los puntos de intersección
intersection_x = []
intersection_y = []

for λ in λ_vals
    θ = range(0, 2π, length=200)
    xc = sqrt(λ) * cos.(θ)  # Coordenadas x de la circunferencia
    yc = sqrt(λ) * sin.(θ)  # Coordenadas y de la circunferencia
    
    # Dibujar la circunferencia en gris
    lines!(ax, xc, yc, color=:gray, linewidth=1, alpha=0.5)
    
    # Resolver la intersección con la recta y = x + 5
    a, b, c = 2, 10, (25 - λ)
    Δ = b^2 - 4a*c  # Discriminante de la ecuación cuadrática
    
    if Δ ≥ 0  
        x1 = (-b + sqrt(Δ)) / (2a)
        x2 = (-b - sqrt(Δ)) / (2a)
        y1 = x1 + 5
        y2 = x2 + 5
        
        # Guardar los puntos de intersección
        push!(intersection_x, x1)
        push!(intersection_y, y1)
        push!(intersection_x, x2)
        push!(intersection_y, y2)

        # Marcar los puntos de intersección en la circunferencia
        scatter!(ax, [x1, x2], [y1, y2], color=:red, markersize=8)
    end
end

# Dibujar la recta en negro
lines!(ax, x_range, x_range .+ 5, color=:black, linewidth=1, linestyle=:dash, label="y = x + 5")

# Dibujar la curva de los puntos de intersección
lines!(ax, intersection_x, intersection_y, color=:blue, linewidth=2, label="Curva de intersección")

axislegend(ax)
circumference
