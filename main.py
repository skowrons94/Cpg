from Simulation import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

np.seterr(all="ignore")

plt.rcParams['figure.figsize'] = [12, 12]
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)

style = {'description_width': 'initial'}

sim = Simulation( )

def findRange( array, value1, value2 ):
    idx1 = (np.abs(array[:,0] - value1)).argmin()
    idx2 = (np.abs(array[:,0] - value2)).argmin()
    return idx1, idx2
    
def findMarker( array, value ):
    idx = (np.abs(array[:,0] - value)).argmin()
    return [value, array[idx][1]]

def Plot( energy = 380, deltaE = 15, position = 1.35, current = 400, time = 2,
          background = 'Underground' ):
    
    time *= 60
    sim.setBackground( background )
    sim.run( energy, deltaE, position, current, time )
    
    gs = gridspec.GridSpec(2, 3,height_ratios=[0.7,2])
    ax2 = plt.subplot( gs[0] )
    ax3 = plt.subplot( gs[1] )
    ax4 = plt.subplot( gs[2] )
    ax1 = plt.subplot( gs[3:] )

    posMarker = findMarker( sim.effPlot, position )
    idx1, idx2 = findRange( sim.crossPlot, sim.Reaction.getCM( energy ),
                            sim.Reaction.getCM( energy - deltaE ) )

    Min = sim.roiMin
    Max = sim.roiMax
    
    ax1.step( sim.x[10:Min], sim.spectrum[10:Min], linewidth=0.9,
              color='tab:blue' )
    ax1.step( sim.x[Min-1:Max], sim.spectrum[Min-1:Max], linewidth=0.9,
              color='tab:red' )
    ax1.step( sim.x[Max-1:-10], sim.spectrum[Max-1:-10], linewidth=0.9,
              color='tab:blue')
    ax1.set_xlabel( "Energy (keV)", fontsize=20, horizontalalignment='right',
                    x=1.0 )
    ax1.set_ylabel( "Counts", fontsize=20, horizontalalignment='right', y=1.0 )
    ax1.set_ylim( bottom=0 ) 

    ax2.plot( sim.effPlot[:,0], sim.effPlot[:,1] )
    ax2.scatter( posMarker[0], posMarker[1], color='red', s=70 )
    ax2.set_xlabel( "Position (cm)", fontsize=13, horizontalalignment='right',
                    x=1.0 )
    ax2.set_ylabel( "Efficiency", fontsize=13, horizontalalignment='right',
                    y=1.0 )
    ax2.tick_params( axis = 'both', which = 'major', labelsize = 13 )
    
    ax3.plot( sim.crossPlot[:,0], sim.crossPlot[:,1] )
    ax3.plot( sim.crossPlot[idx2:idx1,0], sim.crossPlot[idx2:idx1,1],
              color='red', linewidth=4 )
    ax3.set_xlabel( "Energy (keV)", fontsize=13, horizontalalignment='right',
                    x=1.0 )
    ax3.set_ylabel( "Cross-section (b)", fontsize=13,
                    horizontalalignment='right', y=1.0 )
    ax3.tick_params( axis = 'both', which = 'major', labelsize = 13 )
    ax3.set_yscale( "log" )
                   
    ax4.plot( sim.profilePlot[:,0] - sim.Q - sim.Reaction.getCM( energy ) +
              energy, sim.profilePlot[:,1] )
    ax4.set_xlabel( "Beam energy (keV)", fontsize=13,
                    horizontalalignment='right', x=1.0 )
    ax4.set_ylabel( "Target", fontsize=13, horizontalalignment='right', y=1.0 )
    ax4.tick_params( axis = 'both', which = 'major', labelsize = 13 )
    
    plt.tight_layout( )
