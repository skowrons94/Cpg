import numpy as np

from Reaction import *
from Efficiency import *

class Simulation( ):
    def __init__( self ):
        self.Q = 1943.5
        self.sigma = 3
        self.energy = 0
        self.deltaE = 0
        self.pos = 0
        self.eff = Efficiency( )
        self.Reaction = Reaction( )

    e = 0.3
    c = 0.5
    b = 2000
    d = 10

    back = 1
    q = 1.602176634*pow( 10, -19 )
        
    def gaussian( self, x ):
        return 1./np.sqrt( 2. * np.pi * self.sigma*self.sigma ) * np.exp( -x**2 / ( self.sigma*self.sigma*2 ) )

    def P( self, x ):
        return (1/(1 + np.exp((x - self.b)/self.e)))*(1/(1 + np.exp(( self.b - self.d - x)/self.c)))

    def createPeak( self ):        
        self.spectrum = np.zeros( len( self.x ) )
        for idx in range( len( self.x ) ):
            energyCM = self.x[idx] - self.Q
            energyLab = self.Reaction.getLab( energyCM )
            if( energyCM < 10 ):
                self.spectrum[idx] = self.back
            elif( self.x[idx] - self.b > 5 ):
                self.spectrum[idx] = self.back
            else:
                stopCM = self.Reaction.getCM( self.Reaction.stopGraph.Eval( energyLab ) )
                self.spectrum[idx] = self.P(self.x[idx])*self.Reaction.getCross( energyLab )*self.eff.effPeak( self.x[idx]/1000, self.pos )*self.Np/stopCM + self.back

        self.spectrum = np.random.poisson( self.spectrum )

        gaus = np.fromiter( ( self.gaussian( x )
                              for x in range( int( -5*self.sigma ), int( 5*self.sigma ), 1 ) ), np.float )

        self.spectrum = np.convolve( self.spectrum, gaus, mode="same" )
        for idx in range( len( self.spectrum ) ):
            self.spectrum[idx] = int( self.spectrum[idx] )

    def getCross( self ):
        self.crossPlot = np.zeros( shape=[590,2] )
        for idx in range( 590 ):
            self.crossPlot[idx][0] = idx + 10
            self.crossPlot[idx][1] = self.Reaction.getCross( idx + 10 )
            
    def getEff( self ):
        energy = self.Q + self.Reaction.getCM( self.energy ) 
        self.effPlot = np.zeros( shape=[200,2] )
        for idx in range( 200 ):
            self.effPlot[idx][0] = 0.1*idx + 0.1
            self.effPlot[idx][1] = self.eff.effPeak( energy/1000, 0.1*idx + 0.1 )
    def getProfile( self ):
        self.profilePlot = np.zeros( shape=( len( self.x ), 2 ) )
        for idx in range( len( self.x ) ):
            energyCM = self.x[idx] - self.Q
            energyLab = self.Reaction.getLab( energyCM )
            if( energyCM < 10 ):
                self.profilePlot[idx][0] = self.x[idx]
                self.profilePlot[idx][1] = 0
            elif( self.x[idx] - self.b > 5 ):
                self.profilePlot[idx][0] = self.x[idx]
                self.profilePlot[idx][1] = 0
            else:
                self.profilePlot[idx][0] = self.x[idx]
                self.profilePlot[idx][1] = self.P( self.x[idx] )
            

    def run( self, energy, deltaE, pos, current, time ):
        fRunEff = False
        fRunReaction = False
        
        if( self.pos != pos ):
            self.pos = pos
            fRunEff = True

        if( self.energy != energy ):
            self.energy = energy
            self.b = self.Reaction.getCM( energy ) + self.Q
            self.getEff( )
            self.getCross( )
            fRunReaction = True
            
        if( self.deltaE != deltaE ):
            self.deltaE = deltaE
            fRunReaction = True

        self.d = self.Reaction.getCM( deltaE )*self.Reaction.stopGraph.Eval( self.energy )/self.Reaction.stopGraph.Eval( 380 )

        if( fRunEff ):
            self.eff.run( self.pos )

        if( fRunReaction ):
            self.Reaction.run( deltaE, self.energy )

        plotMin = self.Q + self.Reaction.getCM( self.energy ) - 200
        plotMax = self.Q + self.Reaction.getCM( self.energy ) + 50
        self.x = np.linspace( plotMin, plotMax, int( plotMax - plotMin ) )
        self.Np = time*current*pow( 10, -6 )/self.q
        self.createPeak( )
        self.getProfile( )
