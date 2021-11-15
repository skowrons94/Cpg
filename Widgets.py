from Plot import *
from ipywidgets import *

style = {'description_width': 'initial'}

def createWidgets( ):

    energy = widgets.FloatSlider(value=380.0, min=10, max=600, step=1, 
                                 continuous_update=False, 
                                 description='Beam energy (keV):', 
                                 style=style)

    deltaE = widgets.FloatSlider(value=15.0, min=1, max=200, step=1, 
                                 continuous_update=False, 
                                 description='Target thickness (keV):',
                                 style=style)

    position = widgets.FloatSlider(value=1.35, min=0.1, max=20, step=0.01, 
                                   continuous_update=False, 
                                   description='Detector distance (cm):',
                                   style=style)
                  
    current = widgets.FloatSlider(value=400.0, min=1, max=600, step=1, 
                                  continuous_update=False, 
                                  description='Beam current ($\mu$A):',
                                  style=style)

    time = widgets.FloatSlider(value=2, min=1, max=4320, step=1, 
                               continuous_update=False, 
                               description='Acquisition time (min):',
                               style=style)
                   
    background = widgets.Select(options=['Underground Shielded',
                                         'Underground Unshielded',
                                         'Surface Unshielded'],
                                value='Underground Shielded',
                                description='Background:',
                                disabled=False)

    layout = Layout(border='2px solid grey',
                    width='940px',
                    height='',
                    flex_flow='row',
                    display='flex',
                    justify_content = 'center')

    left_box = VBox([energy, deltaE, position, current, time])
    right_box = VBox([background])
    ui = HBox([left_box, right_box], layout=layout)

    w = interactive_output(Plot,
                           { 
                               "energy":energy,
                               "deltaE":deltaE,
                               "position":position,
                               "current":current,
                               "time":time,
                               "background":background,
                           })

    return ui, w

