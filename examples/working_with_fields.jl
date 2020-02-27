using OceanTurb, Printf

@use_pyplot_utils

# # Operations on fields
#
# Generate a grid with 100 points and height H=10. By default the
# base of the grid is at z=-H and the top of the grid is at z=0.
grid = UniformGrid(100, 10)

# Define a function that will generate our field's data.
c0(z) = exp(1.2z)

# Create a field of zeros on the grid.
c = CellField(grid)

# Set c to the function c0.
OceanTurb.set!(c, c0)

# Calculate the derivative of c and set it to cz.
cz = ∂z(c)

# Plot the results
plot(c, label=L"c = e^{1.2z}")
plot(cz, label=L"\partial_z c = 1.2 e^{1.2z}")
legend()

title("A field and it's derivative")
xlabel("\$ c \$ and \$ \\partial_z c \$")
ylabel(L"z")


# # Interpolation
#
# Next we demonstrate how to interpolate from one grid to another, 
# respecting the integral budget of a quantity.

grid1 = UniformGrid(30, 10)
grid2 = UniformGrid(10, 10)
grid3 = UniformGrid(3, 10)

c1 = CellField(grid1)
c2 = CellField(grid2)
c3 = CellField(grid3)

OceanTurb.set!(c1, c)
OceanTurb.set!(c2, c)
OceanTurb.set!(c3, c)

@show integral(c)
@show integral(c1)
@show integral(c2)
@show integral(c3)

plot(c, ".", label="high-res \$ c \$")
plot(c1, "^", label="medium-res \$ c \$")
plot(c2, "x", label="low-res \$ c \$")
plot(c3, "+", label="crazy low-res \$ c \$")

legend()

title("\$ e^{1.2z} \$ at different resolutions")
xlabel(L"c")
ylabel(L"z")
