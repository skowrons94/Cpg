from Simulation import *

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

sim = Simulation( )

def plot( ):

    sim.run( 380, 9, 1.35, 400, 135 )

    plt.step( sim.x, sim.spectrum, linewidth=0.5 )
    plt.savefig( "spectrum.pdf" )

plot( )
