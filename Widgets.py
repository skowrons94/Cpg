from ipywidgets import interactive_output, widgets, HBox, VBox, Layout

def createWidgets( run ):
    style = {'description_width': 'initial'}

    energy = widgets.FloatSlider(value=380.0, min=10, max=600, step=1, 
                                 continuous_update=False, 
                                 description='Beam energy (keV):', 
                                 style=style)
    
    deltaE = widgets.FloatSlider(value=150, min=10, max=5000, step=10, 
                                 continuous_update=False, 
                                 description='Target thickness (nm):',
                                 style=style)
    
    position = widgets.FloatSlider(value=1.35, min=0.1, max=20, step=0.01, 
                                   continuous_update=False, 
                                   description='Detector distance (cm):',
                                   style=style)
                  
    current = widgets.FloatSlider(value=100.0, min=0, max=500, step=1, 
                                  continuous_update=False, 
                                  description='Beam current ($\mu$A):',
                                  style=style)

    time = widgets.FloatSlider(value=120, min=2, max=4320, step=1, 
                               continuous_update=False, 
                               description='Acquisition time (min):',
                               style=style)

    sigma = widgets.FloatSlider(value=3, min=2, max=40, step=0.1, 
                                continuous_update=False, 
                                description='Detector $\sigma$ (keV):',
                                style=style)
                   
    background = widgets.Select(options=['Underground Shielded',
                                         'Underground Unshielded',
                                         'Surface Unshielded'],
                                value='Underground Shielded',
                                description='Background:',
                                disabled=False)

    scale = widgets.Select(options=['Linear',
                                    'Log'],
                           value='Linear',
                           description='Plot Scale:',
                           disabled=False)

    layout = Layout(border='2px solid grey',
                    width='1050px',
                    height='',
                    flex_flow='row',
                    display='flex',
                    justify_content = 'center')

    left_box = VBox([energy, deltaE, position, current, time, sigma])
    right_box = VBox([background, scale])

    ui = HBox([left_box, right_box], layout=layout)
    w = interactive_output(run,
                           { 
                               "energy":energy,
                               "deltaE":deltaE,
                               "position":position,
                               "current":current,
                               "time":time,
                               "background":background,
                               "scale":scale,
                               "sigma":sigma
                           })

    return ui, w
