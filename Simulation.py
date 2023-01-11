import pickle

import numpy as np

from scipy import signal

from Edge import Edge
from Reaction import Reaction
from Efficiency import Efficiency

np.seterr( all="ignore" )

class Simulation( ):
    def __init__( self ):
        self.sigma    = 2
        self.energy   = 0
        self.deltaE   = 0
        self.pos      = 0
        self.eff      = Efficiency( )
        self.Reaction = Reaction( )
        self.getBack( )

    Q = 1943.5

    q = 1.602176634*pow( 10, -19 )
        
    e = 0.3
    c = 0.5
    b = 2000
    d = 10

    def getBack( self ):
        with open( "./data/SurfaceUnshielded.pkl", "rb" ) as f:
            self.surfaUnshield = pickle.load( f )
        with open( "./data/UndergroundUnshielded.pkl", "rb" ) as f:
            self.underUnshield = pickle.load( f )
        with open( "./data/UndergroundShielded.pkl", "rb" ) as f:
            self.underShielded = pickle.load( f )
        
    def gaussian( self, x, sigma ):
        return 1./np.sqrt( 2. * np.pi * sigma*sigma ) * np.exp( -x**2 / ( sigma*sigma*2 ) )

    def P( self, x ):
        return (1/(1 + np.exp((x - self.b)/self.e)))*(1/(1 + np.exp(( self.b - self.d - x)/self.c)))

    def setBackground( self, string ):
        if( string == "Underground Unshielded" ):
            self.envBackX = self.underUnshield["x"]
            self.envBackY = self.underUnshield["y"]
        elif( string == "Underground Shielded" ):
            self.envBackX = self.underShielded["x"]
            self.envBackY = self.underShielded["y"]
        elif( string == "Surface Unshielded" ):
            self.envBackX = self.surfaUnshield["x"]
            self.envBackY = self.surfaUnshield["y"]

    def setSigma( self, sigma ):
        self.sigma = sigma

    def setComBack( self ):
        E_g           = self.Q + self.energy
        self.compEdge = 2*pow( E_g, 2 )/( 2*E_g + 511 )
        ratio         = self.eff.effPeak( ( self.Q + self.energy )/1000, self.pos )/self.eff.effTot( ( self.Q + self.energy )/1000, self.pos )/E_g
        self.comBack  = self.Reaction.Yield*self.Np*self.eff.effPeak( ( self.Q + self.energy )/1000, self.pos )*ratio

    def getCompton( self ):
        self.spectrumComp = np.zeros( len( self.x ) )
        for idx in range( len( self.x ) ):
            if( self.x[idx] > self.Q + self.energy ):
                self.spectrumComp[idx] = 0
            elif( self.x[idx] < 50 ):
                self.spectrumComp[idx] = 0
            elif( self.x[idx] > self.compEdge - 50 ):
                self.spectrumComp[idx] = self.comBack*self.eff.effPeak( self.x[idx]/1000, self.pos )*Edge( self.x[idx], self.compEdge + 30, 30 )
            else:
                self.spectrumComp[idx] = self.comBack*self.eff.effPeak( self.x[idx]/1000, self.pos )               
        
    def createPeak( self, time ):
        self.x            = np.linspace( 0, 6000, 6000 )
        self.spectrum     = np.zeros( len( self.x ) )
        self.spectrumBack = np.zeros( len( self.x ) )

        gaus  = np.fromiter( ( self.gaussian( x, self.sigma ) for x in range( int( -5*self.sigma ), int( 5*self.sigma ), 1 ) ), float )
        gaus /= sum( gaus )

        for idx in range( len( self.x ) ):
            energyCM  = self.x[idx] - self.Q
            energyLab = self.Reaction.getLab( energyCM )
            stopCM    = self.Reaction.getCM( self.Reaction.getValue( self.Reaction.stopGraph, energyLab ) )
            
            if( energyCM < 10 ):
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                self.spectrumBack[idx] = self.envBackY[idxBack]*time

            elif( self.x[idx] - self.b > 10 ):
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                self.spectrumBack[idx] = self.envBackY[idxBack]*time

            elif( self.x[idx] > self.b - self.d/2 ):
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                self.spectrumBack[idx] = self.envBackY[idxBack]*time
                self.spectrum[idx] = self.P(self.x[idx])*self.Reaction.getCross( energyLab )*self.eff.effPeak( self.x[idx]/1000, self.pos )*self.Np/stopCM

            else:
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                self.spectrumBack[idx] = self.envBackY[idxBack]*time
                self.spectrum[idx] = self.P(self.x[idx])*self.Reaction.getCross( energyLab )*self.eff.effPeak( self.x[idx]/1000, self.pos )*self.Np/stopCM

        self.getCompton( )
                
        self.spectrum     = np.random.poisson( self.spectrum     )
        self.spectrumBack = np.random.poisson( self.spectrumBack )
        self.spectrumComp = np.random.poisson( self.spectrumComp )

        self.spectrum     = np.convolve( self.spectrum,     gaus, mode="same" )
        self.spectrumComp = np.convolve( self.spectrumComp, gaus, mode="same" )

#        if( self.sigma > 2 ):
#            gausBack = np.fromiter( ( self.gaussian( x, self.sigma - 2 ) for x in range( int( -5*self.sigma ), int( 5*self.sigma ), 1 ) ), float )
#            self.spectrumBack = np.convolve( self.spectrumBack, gausBack, mode="same" )

        for idx in range( len( self.spectrum ) ):
            self.spectrum[idx]     = int( self.spectrum[idx]                              )
            self.spectrumBack[idx] = int( self.spectrumComp[idx] + self.spectrumBack[idx] )            

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
            elif( self.x[idx] - self.b > 10 ):
                self.profilePlot[idx][0] = self.x[idx]
                self.profilePlot[idx][1] = 0
            else:
                self.profilePlot[idx][0] = self.x[idx]
                self.profilePlot[idx][1] = self.P( self.x[idx] )

    def convertDeltaE( self, deltaE ):
        conv = self.Reaction.getValue( self.Reaction.stopGraphUM, self.energy )
        return deltaE*conv
        
    def run( self, energy, deltaEUM, pos, current, time ):
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
            
        deltaE = self.convertDeltaE( deltaEUM )
        if( deltaE > energy ):
            deltaE = energy
        self.deltaE  = deltaE
        fRunReaction = True

        self.d  = self.Reaction.convertDeltaE( self.deltaE, self.energy ) 
        self.Np = time*current*pow( 10, -6 )/self.q
        
        if( fRunEff ):
            self.eff.run( self.pos )

        if( fRunReaction ):
            self.Reaction.run( self.deltaE, self.energy )

        self.setComBack( )

        self.createPeak( time )
        self.getProfile( )
