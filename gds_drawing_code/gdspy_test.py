


import gdspy
import numpy

# The GDSII file is called a library, which contains multiple cells.
lib = gdspy.GdsLibrary()

# Geometry must be placed in cells.
cell = lib.new_cell('FIRST')

# Create the geometry (a single rectangle) and add it to the cell.
rect = gdspy.Rectangle((0, -2), (2, -1))
# cell.add(rect)

# Manually connect the hole to the outer boundary
cutout = gdspy.Polygon(
    [(0, 0), (5, 0), (5, 5), (0, 5), (0, 0), (2, 2), (2, 3), (3, 3), (3, 2), (2, 2)]
)
# cell.add(cutout)

# Circle centered at (0, 0), with radius 2 and tolerance 0.1
circle = gdspy.Round((0, 0), 2, tolerance=1e-3)
# cell.add(circle)

# Circular arc example
arc = gdspy.Round(
    (2, 4),
    2,
    inner_radius=1,
    initial_angle=-0.2 * numpy.pi,
    final_angle=1.2 * numpy.pi,
    tolerance=1e-3,
)
# cell.add(arc)

cell.add(rect)

cell2 = lib.new_cell('SECOND')
cell2.add(cutout)

# Save the library in a file called 'first.gds'.
lib.write_gds('first.gds')

# Display all cells using the internal viewer.
gdspy.LayoutViewer()