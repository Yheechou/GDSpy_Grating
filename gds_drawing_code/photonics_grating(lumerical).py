######################################################################
#                                                                    #
#  Copyright 2009-2019 Lucas Heitzmann Gabrielli.                    #
#  This file is part of gdspy, distributed under the terms of the    #
#  Boost Software License - Version 1.0.  See the accompanying       #
#  LICENSE file or <http://www.boost.org/LICENSE_1_0.txt>            #
#                                                                    #
######################################################################

import numpy
import gdspy
from d2nn_construct import d2nn_construct


def grating_demo(
    period,
    number_of_teeth,
    fill_frac,
    width,
    position,
    direction,
    lda=1,
    sin_theta=0,
    focus_distance=-1,
    focus_width=-1,
    tolerance=0.001,
    layer=0,
    datatype=0,
):
    """
    Straight or focusing grating.

    period          : grating period
    number_of_teeth : number of teeth in the grating
    fill_frac       : filling fraction of the teeth (w.r.t. the period)
    width           : width of the grating
    position        : grating position (feed point)
    direction       : one of {'+x', '-x', '+y', '-y'}
    lda             : free-space wavelength
    sin_theta       : sine of incidence angle
    focus_distance  : focus distance (negative for straight grating)
    focus_width     : if non-negative, the focusing area is included in
                      the result (usually for negative resists) and this
                      is the width of the waveguide connecting to the
                      grating
    tolerance       : same as in `path.parametric`
    layer           : GDSII layer number
    datatype        : GDSII datatype number

    Return `PolygonSet`
    """
    if focus_distance < 0:
        p = gdspy.L1Path(
            (
                position[0] - 0.5 * width,
                position[1] + 0.5 * (number_of_teeth - 1 + fill_frac) * period,
            ),
            "+x",
            period * fill_frac,
            [width],
            [],
            number_of_teeth,
            period,
            layer=layer,
            datatype=datatype,
        )
    else:
        neff = lda / float(period) + sin_theta
        qmin = int(focus_distance / float(period) + 0.5)
        p = gdspy.Path(period * fill_frac, position)
        c3 = neff ** 2 - sin_theta ** 2
        w = 0.5 * width
        for q in range(qmin, qmin + number_of_teeth):
            c1 = q * lda * sin_theta
            c2 = (q * lda) ** 2
            p.parametric(
                lambda t: (
                    width * t - w,
                    (c1 + neff * numpy.sqrt(c2 - c3 * (width * t - w) ** 2)) / c3,
                ),
                tolerance=tolerance,
                max_points=0,
                layer=layer,
                datatype=datatype,
            )
            p.x = position[0]
            p.y = position[1]
        sz = p.polygons[0].shape[0] // 2
        if focus_width == 0:
            p.polygons[0] = numpy.vstack((p.polygons[0][:sz, :], [position]))
        elif focus_width > 0:
            p.polygons[0] = numpy.vstack(
                (
                    p.polygons[0][:sz, :],
                    [
                        (position[0] + 0.5 * focus_width, position[1]),
                        (position[0] - 0.5 * focus_width, position[1]),
                    ],
                )
            )
        p.fracture()
    if direction == "-x":
        return p.rotate(0.5 * numpy.pi, position)
    elif direction == "+x":
        return p.rotate(-0.5 * numpy.pi, position)
    elif direction == "-y":
        return p.rotate(numpy.pi, position)
    else:
        return p


def grating_lumerical(
    period,
    number_of_teeth,
    fill_frac,
    width,
    position,
    direction,
    lda=1,
    sin_theta=0,
    focus_distance=-1,
    focus_width=-1,
    tolerance=0.001,
    layer=0,
    datatype=0,
):
    """
    Straight or focusing grating.

    period          : grating period
    number_of_teeth : number of teeth in the grating
    fill_frac       : filling fraction of the teeth (w.r.t. the period)
    width           : width of the grating
    position        : grating position (feed point)
    direction       : one of {'+x', '-x', '+y', '-y'}
    lda             : free-space wavelength
    sin_theta       : sine of incidence angle
    focus_distance  : focus distance (negative for straight grating)
    focus_width     : if non-negative, the focusing area is included in
                      the result (usually for negative resists) and this
                      is the width of the waveguide connecting to the
                      grating
    tolerance       : same as in `path.parametric`
    layer           : GDSII layer number
    datatype        : GDSII datatype number

    Return `PolygonSet`
    """
    if focus_distance < 0:
        p = gdspy.L1Path(
            (
                position[0] - 0.5 * width,
                position[1] + 0.5 * (number_of_teeth - 1 + fill_frac) * period,
            ),
            "+x",
            period * fill_frac,
            [width],
            [],
            number_of_teeth,
            period,
            layer=layer,
            datatype=datatype,
        )
    else:
        p = gdspy.Path(period * fill_frac, position)
        for q in range(number_of_teeth):
            
            c1 = q*period + focus_distance
            c2 = c1/focus_distance*focus_width/2 # tan(focusing angle) = focus_width/focus_distance
            
            p.parametric(
                lambda t: (
                    c2 * 2 * t - c2,
                    numpy.sqrt(c1**2 - (c2 * 2 * t - c2)**2),
                ),
                tolerance=tolerance,
                max_points=0,
                layer=layer,
                datatype=datatype,
            )
            p.x = position[0]
            p.y = position[1]
        sz = p.polygons[0].shape[0] // 2
        # if focus_width == 0:
        #     p.polygons[0] = numpy.vstack((p.polygons[0][:sz, :], [position]))
        # elif focus_width > 0:
        #     p.polygons[0] = numpy.vstack(
        #         (
        #             p.polygons[0][:sz, :],
        #             [
        #                 (position[0] + 0.5 * focus_width, position[1]),
        #                 (position[0] - 0.5 * focus_width, position[1]),
        #             ],
        #         )
        #     )
        p.fracture()
    if direction == "-x":
        return p.rotate(0.5 * numpy.pi, position)
    elif direction == "+x":
        return p.rotate(-0.5 * numpy.pi, position)
    elif direction == "-y":
        return p.rotate(numpy.pi, position)
    else:
        return p


if __name__ == "__main__":
    # Examples
    lib = gdspy.GdsLibrary()


    # Positive resist example
    width = 0.5
    ring_radius = 20.0
    big_margin = 10.0
    small_margin = 5.0
    taper_len = 50.0
    bus_len = 2000.0#waveguide length
    
    io_gap = 500.0
    wg_gap = 20.0
    ring_gaps = [0.06 + 0.02 * i for i in range(8)]

    p_demo = gdspy.Path(
        small_margin, (0, 0), number_of_paths=2, distance=small_margin + width
    )
    p_demo.segment(21.5, "+y", final_distance=small_margin + 19)

    marker = gdspy.Path(
        small_margin, (0, 0), number_of_paths=2, distance=small_margin*3 + width
    )
    marker.segment(small_margin, "-y")


    gratsur_demo = lib.new_cell("PGRATSur_demo")
    gratsur_demo.add(p_demo)

    grat_demo = lib.new_cell("PGRAT_demo")
    #grat_demo.add(p_demo)
    #grat_demo.add(marker)#marker
    grat_demo.add(
        grating_demo(
            0.75,
            28,
            0.28,
            19,
            (0, 0),
            "+y",
            1.55,
            numpy.sin(numpy.pi * 10 / 180),
            21.5,
            tolerance=0.001,
            layer=1
        )
    )

    ##########lumerical grating
    p_lumerical = gdspy.Path(
        small_margin, (0, 0), number_of_paths=2, distance=small_margin + width
    )
    p_lumerical.segment(40, "+y", final_distance=small_margin + 40)

    gratsur_lumerical = lib.new_cell("PGratSur_lumerical")
    gratsur_lumerical.add(p_lumerical)

    marker = gdspy.Path(
        small_margin, (0, 0), number_of_paths=2, distance=small_margin*3 + width
    )
    marker.segment(small_margin, "-y")
    grat_lumerical = lib.new_cell("PGrat_lumerical")
    #grat_lumerical.add(p_lumerical)
    #grat_lumerical.add(marker)
    grat_lumerical.add(
        grating_lumerical(
            0.75,
            28,
            0.28,
            19,
            (0, 0),
            "+y",
            1.55,
            numpy.sin(numpy.pi * 10 / 180),
            21.5,
            20,
            tolerance=0.001,
            layer=1
        )
    )


    c = lib.new_cell("Positive")

    y_min = 4550 #minimum position of structure

    layer_distance = 200
    input_distance = 100
    num_layers = 5
    wg_len = 300
    x_max = 0
    y_offset = 3500
    polygon_layer=2

    y_min = 0
    filepath = '0_1\\'
    c = d2nn_construct(c, filepath, x_max, y_min, layer_distance, input_distance, num_layers, grat_lumerical, gratsur_lumerical, small_margin, wg_len, polygon_layer=polygon_layer)

    y_min = y_min + y_offset
    filepath = '1_6\\'
    c = d2nn_construct(c, filepath, x_max, y_min, layer_distance, input_distance, num_layers, grat_lumerical, gratsur_lumerical, small_margin, wg_len, polygon_layer=polygon_layer)

    y_min = y_min + y_offset
    filepath = '1_7\\'
    c = d2nn_construct(c, filepath, x_max, y_min, layer_distance, input_distance, num_layers, grat_lumerical, gratsur_lumerical, small_margin, wg_len, polygon_layer=polygon_layer)





    ########################
    ###left waveguide


    input_gap = 3000.0# distance between parallel waveguides

    x_offset = -3000-1400
    #waveguide
    for i in range(2):
        path = gdspy.FlexPath(
            [(x_offset-input_gap * i, 0)],
            width= [small_margin, small_margin],
            offset= small_margin + width,
            gdsii_path=True,
        )
        path.segment((0, bus_len), relative=True)
        c.add(path)
    #grating
    c.add(
        gdspy.CellArray(
            grat_lumerical, 1, 1, (-input_gap, 0), (x_offset, bus_len)
        )
    )
    c.add(
        gdspy.CellArray(
            grat_lumerical, 1, 1, (input_gap, 0), (x_offset, 0), 180
        )
    )
    #surreounding of grating
    c.add(
        gdspy.CellArray(
            gratsur_lumerical, 1, 1, (-input_gap, 0), (x_offset, bus_len)
        )
    )
    c.add(
        gdspy.CellArray(
            gratsur_lumerical, 1, 1, (input_gap, 0), (x_offset, 0), 180
        )
    )

    #grating
    c.add(
        gdspy.CellArray(
            grat_demo, 1, 1, (-input_gap, 0), (x_offset - input_gap, bus_len)
        )
    )
    c.add(
        gdspy.CellArray(
            grat_demo, 1, 1, (input_gap, 0), (x_offset - input_gap, 0), 180
        )
    )

    #surrounding of grating
    c.add(
        gdspy.CellArray(
            gratsur_demo, 1, 1, (-input_gap, 0), (x_offset - input_gap, bus_len)
        )
    )
    c.add(
        gdspy.CellArray(
            gratsur_demo, 1, 1, (input_gap, 0), (x_offset - input_gap, 0), 180
        )
    )
    # Save to a gds file and check out the output
    lib.write_gds("test.gds")
    gdspy.LayoutViewer(lib)
