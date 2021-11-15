import pickle

import numpy as np

from Reaction import *
from Efficiency import *

np.seterr( all="ignore" )

class Simulation( ):
    def __init__( self ):
        self.Q = 1943.5
        self.sigma = 3
        self.energy = 0
        self.deltaE = 0
        self.pos = 0
        self.eff = Efficiency( )
        self.getBack( )
        self.Reaction = Reaction( )

    e = 0.3
    c = 0.5
    b = 2000
    d = 10

    comBack = 0
    envBack = pow( 10, -3 )/36
    q = 1.602176634*pow( 10, -19 )

    def getBack( self ):
        with open( "./data/SurfaceUnshielded.pkl", "rb" ) as f:
            self.surfaUnshield = pickle.load( f )
        with open( "./data/UndergroundUnshielded.pkl", "rb" ) as f:
            self.underUnshield = pickle.load( f )
        with open( "./data/UndergroundShielded.pkl", "rb" ) as f:
            self.underShielded = pickle.load( f )
        
    def gaussian( self, x ):
        return 1./np.sqrt( 2. * np.pi * self.sigma*self.sigma ) * np.exp( -x**2 / ( self.sigma*self.sigma*2 ) )

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

    def setComBack( self ):
        self.comBack = self.Reaction.Yield*self.Np*self.eff.effPeak( ( self.Q + self.energy )/1000, self.pos )*0.001

    def setRange( self, intMin = 0, intMax = 0 ):
        if( intMin == 0 and intMax == 0 ):
            plotMin = self.Q + self.Reaction.getCM( self.energy ) - 200
            plotMax = self.Q + self.Reaction.getCM( self.energy ) + 50
            self.x = np.linspace( plotMin, plotMax, int( plotMax - plotMin ) )

        else:
            self.x = np.linspace( intMin, intMax, int( intMax - intMin ) )

    def createPeak( self, time ):
        gaus = np.fromiter( ( self.gaussian( x ) for x in range( int( -5*self.sigma ), int( 5*self.sigma ), 1 ) ), np.float )
        self.spectrum = np.zeros( len( self.x ) )
        for idx in range( len( self.x ) ):
            energyCM = self.x[idx] - self.Q
            energyLab = self.Reaction.getLab( energyCM )
            stopCM = self.Reaction.getCM( self.Reaction.getValue( self.Reaction.stopGraph, energyLab ) )

            if( self.x[idx] < self.b + 10 ):
                self.roiMax = idx
                
            if( energyCM < 10 ):
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                self.spectrum[idx] = self.envBackY[idxBack]*time + self.comBack

            elif( self.x[idx] - self.b > 10 ):
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                self.spectrum[idx] = self.envBackY[idxBack]*time

            elif( self.x[idx] > self.b - self.d/2 ):
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                self.spectrum[idx] = self.P(self.x[idx])*self.Reaction.getCross( energyLab )*self.eff.effPeak( self.x[idx]/1000, self.pos )*self.Np/stopCM + self.envBackY[idxBack]*time

            else:
                idxBack = np.abs( self.envBackX - self.x[idx] ).argmin( )
                if( self.P( self.x[idx] ) < 1e-5 ):
                    self.roiMin = idx
                self.spectrum[idx] = self.P(self.x[idx])*self.Reaction.getCross( energyLab )*self.eff.effPeak( self.x[idx]/1000, self.pos )*self.Np/stopCM + self.envBackY[idxBack]*time + self.comBack
                
        self.spectrum = np.random.poisson( self.spectrum )
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
            elif( self.x[idx] - self.b > 10 ):
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

        self.d = self.Reaction.convertDeltaE( self.deltaE, self.energy ) 
        self.Np = time*current*pow( 10, -6 )/self.q
        
        if( fRunEff ):
            self.eff.run( self.pos )

        if( fRunReaction ):
            self.Reaction.run( self.deltaE, self.energy )

        self.setRange( )
        self.setComBack( )

        self.createPeak( time )
        self.getProfile( )
