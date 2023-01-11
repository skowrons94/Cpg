import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

def getProfile( sim, energy, deltaE ):
    x = [ ]
    y = [ ]
    for i in range( len( sim.profilePlot ) ):
        data = sim.profilePlot[i,0] - sim.Q - sim.Reaction.getCM( energy ) + energy
        if( ( data > energy - deltaE - 30 ) and ( data < energy + 10 ) ):
            x.append( data )
            y.append( sim.profilePlot[i,1] )

    return x, y

def findRange( array, value1, value2 ):
    idx1 = (np.abs(array[:,0] - value1)).argmin()
    idx2 = (np.abs(array[:,0] - value2)).argmin()
    return idx1, idx2
    
def findMarker( array, value ):
    idx = (np.abs(array[:,0] - value)).argmin()
    return [value, array[idx][1]]

def createPlot( sim, energy, deltaE, position, scale ):
    
    gs = gridspec.GridSpec(2, 3,height_ratios=[0.7,2])
    ax2 = plt.subplot( gs[0] )
    ax3 = plt.subplot( gs[1] )
    ax4 = plt.subplot( gs[2] )
    ax1 = plt.subplot( gs[3:] )

    posMarker = findMarker( sim.effPlot, position )
    idx1, idx2 = findRange( sim.crossPlot, sim.Reaction.getCM( energy ),
                            sim.Reaction.getCM( energy - deltaE ) )
    
    line1, = ax1.step( sim.x, sim.spectrum + sim.spectrumBack, linewidth=2,
                       color='tab:red', label="Signal" )
    line2, = ax1.step( sim.x, sim.spectrumBack, linewidth=2,
                       color='tab:blue', label="Background" )
    ax1.set_xlabel( "Energy (keV)", fontsize=20, horizontalalignment='right',
                    x=1.0 )
    ax1.set_xlim( left=0, right=6000 )
    ax1.set_ylabel( "Counts", fontsize=20, horizontalalignment='right', y=1.0 )
    if( scale == "Log" ):
        ax1.set_yscale( "log" )
        ax1.set_ylim( bottom=1e-2 )
    else:
        ax1.set_yscale( "linear" )
        ax1.set_ylim( bottom=0 )
#    ax1.legend([line1, line2], ['Signal', 'Background'], loc="upper right",
#               fontsize=20)

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

    xP, yP = getProfile( sim, energy, deltaE ) 
    
    ax4.plot( xP, yP )
    ax4.set_xlabel( "Beam Energy (keV)", fontsize=13,
                    horizontalalignment='right', x=1.0 )
    ax4.set_ylabel( "Target Profile", fontsize=13, horizontalalignment='right', y=1.0 )
    ax4.tick_params( axis = 'both', which = 'major', labelsize=13 )


    plt.tight_layout( )
    plt.show( )
