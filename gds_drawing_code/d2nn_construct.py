import numpy
import gdspy
import scipy.io as sio
def d2nn_construct(c, filepath, x_max, y_min, layer_distance, input_distance, num_layers, grat, grat_sur, small_margin, wg_len, polygon_layer=0):
    '''
    wg_len: output vertical waveguide length
    '''
    #propagate towards left
    #add structure
    #filepath = 'H:\\tiankuang\\Projects\\chip-wavefront-shaping-D2NN\\modulator-design-200nm-lateral-range\\harmonic-testing-matlab-validation\\1209_5layer_30000pixel_sample(100uminput)\\'
    #post_widths = dict()
    x_offset = x_max-input_distance
    for i in range(num_layers):
        post = sio.loadmat(filepath+'mask_length_0_'+str(i)+'.mat')['save_mask_phase']
        x_start = -i*layer_distance + x_offset
        grid_width = 0.3
        y_offset = y_min
        for i_post in range(3000):
            post_width = post[0, i_post*10+5]
            post_width = post_width*0.05/2
            if post_width>0.01:
                c.add(gdspy.Polygon([(x_start, y_offset+grid_width*i_post+grid_width/2-post_width), (x_start, y_offset+grid_width*i_post+grid_width/2+post_width),
                    (x_start-0.4, y_offset+grid_width*i_post+grid_width/2+post_width), (x_start-0.4, y_offset+grid_width*i_post+grid_width/2-post_width),
                    (x_start, y_offset+grid_width*i_post+grid_width/2-post_width)], layer=polygon_layer))


    #input marker
    m_width = 50
    m_x = x_max + m_width 
    m_y = 0 + y_min
    c.add(gdspy.Polygon([(m_x, m_y), (m_x-m_width, m_y), (m_x-m_width, m_y-m_width), (m_x, m_y-m_width), (m_x, m_y)]))
    m_width = 50
    m_x = x_max + m_width
    m_y = 900 + m_width + y_min
    c.add(gdspy.Polygon([(m_x, m_y), (m_x-m_width, m_y), (m_x-m_width, m_y-m_width), (m_x, m_y-m_width), (m_x, m_y)]))

    m_width = 150
    mk_width = 50
    m_x = x_max + m_width
    c.add(gdspy.Polygon([(m_x, y_min-mk_width), (m_x-m_width+mk_width, y_min-mk_width), \
        (m_x-m_width+mk_width, y_min), (m_x-m_width, y_min), \
        (m_x-m_width, y_min+900), (m_x-m_width+mk_width, y_min+900), \
         (m_x-m_width+mk_width, y_min+900+mk_width), (m_x, y_min+900+mk_width), \
             (m_x, y_min-mk_width)], layer=3))


    #add waveguide
    y_offset = [600+y_min, 300+y_min]
    x_start = -num_layers*layer_distance+x_offset
    bend_radius = 150
    width = 0.5
    wg_horizon = [-300, -300]
    wg_verticl = [wg_len, -wg_len]
    for i in range(2):
        #one time
        path = gdspy.FlexPath(
            [(x_start, y_offset[i])],
            width=[small_margin, small_margin],
            offset = small_margin + width,
            corners="circular bend",
            bend_radius=bend_radius,
            gdsii_path=True,
        )
        #path.segment((0, 600 - wg_gap * i), relative=True)
        path.segment((wg_horizon[i], 0), relative=True)
        path.segment((0, wg_verticl[i]), relative=True)
        c.add(path)

        # #two time
        # path = gdspy.FlexPath(
        #     [(x_start, y_offset[i]+(small_margin + width)/2)],
        #     width=[small_margin],
        #     #offset = small_margin + width,
        #     corners="circular bend",
        #     bend_radius=bend_radius,
        #     gdsii_path=True,
        # )
        # #path.segment((0, 600 - wg_gap * i), relative=True)
        # path.segment((wg_horizon[i], 0), relative=True)
        # path.segment((0, wg_verticl[i]), relative=True)
        # c.add(path)

        # path = gdspy.FlexPath(
        #     [(x_start, y_offset[i]-(small_margin + width)/2)],
        #     width=[small_margin],
        #     #offset = small_margin + width,
        #     corners="circular bend",
        #     bend_radius=bend_radius,
        #     gdsii_path=True,
        # )
        #path.segment((0, 600 - wg_gap * i), relative=True)
        #path.segment((wg_horizon[i], 0), relative=True)
        #path.segment((0, wg_verticl[i]), relative=True)
        c.add(path)


    #grating coupler
    c.add(
        gdspy.CellArray(
            grat, 1, 1, (0, 0), (x_start+wg_horizon[0], y_offset[0]+wg_len) 
        )
    )
    #surrounding of grating
    c.add(
        gdspy.CellArray(
            grat_sur, 1, 1, (0, 0), (x_start+wg_horizon[0], y_offset[0]+wg_len) 
        )
    )
    #grating coupler
    c.add(
        gdspy.CellArray(
            grat, 1, 1, (0, 0), (x_start+wg_horizon[1], y_offset[1]-wg_len), 180
        )
    )
    #surrounding of grating
    c.add(
        gdspy.CellArray(
            grat_sur, 1, 1, (0, 0), (x_start+wg_horizon[1], y_offset[1]-wg_len), 180
        )
    )




    return c