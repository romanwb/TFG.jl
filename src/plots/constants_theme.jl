"""



"""

struct CircularContainer
    v::Vector
end
vector = Vector([1,2,3])
c = CircularContainer([1,2,3])

Base.length(c::CircularContainer) = length(c.v)
Base.getindex(c::CircularContainer, i::Int) = c.v[mod1(i, length(c))]
Base.iterate(c::CircularContainer, state = 1) = iterate(c.v, state)
Base.firstindex(c::CircularContainer) = 1
Base.lastindex(c::CircularContainer) = length(c)

const BLUE = "#4263eb"
const RED = "#f03e3e"
const GREEN = "#0a9a84"
const DARKBLUE = "#191e44"
const ORANGE = "#f76707"
const PURPLE = "#7143e0"
const YELLOW = "#f2cc35"
const DARKPINK = "#791457"
const LIGHTBLUE = "#6adad3"
const PINK = "#cf6db0"


const COLORSCHEME = [
    BLUE, RED, GREEN, DARKBLUE, ORANGE, PURPLE, YELLOW, DARKPINK, LIGHTBLUE, PINK
]



const MARKERSLIST = [
    :circle,
    :utriangle,
    :rect,
    :dtriangle,
    :diamond,
    :xcross,
    :hexagon,
    :star5
]

const GREYPALETTE = [
    "#BBB",
    "#999",
    "#666",
    "#333",
    "#000"
]

const COLORS = CircularContainer(COLORSCHEME)
const GREYS = CircularContainer(GREYPALETTE)
const MARKERS = CircularContainer(MARKERSLIST)

const FONTSIZE = 20
const LABELSIZE = 20
const MARKERSIZE = 10
const LINEWIDTH = 3
const DEFAULT_SIZE = (900, 675)
const FONT = "CMU"
const MINORGRID = false
const GRID = false

