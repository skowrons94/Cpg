import numpy as np
import matplotlib.pyplot as plt

class Reaction( ):
    def __init__( self ):
        self.stop = self.readStop( )
        self.sFactor = self.readSfactor( )
        self.sFactorGraph = self.createGraph( self.sFactor )
        self.stopGraph = self.createGraph( self.stop )

    M0 = 12
    M1 = 1.00727647

    def readSfactor( self ):
        dataDir = "./data/sfactor.dat"

        fIn = open( dataDir, "r" )
        Lines = fIn.readlines( )
        
        sFactor = np.zeros( shape=( len( Lines ), 3 ) )
        for idx in range( len( Lines ) ):
            l = Lines[idx].split( )
            sFactor[idx][0] = float( l[0] )
            sFactor[idx][1] = float( l[1] )

        return sFactor

    def readStop( self ):
        dataDir = "./data/c12.dat"

        fIn = open( dataDir, "r" )
        Lines = fIn.readlines( )
        
        stop = np.zeros( shape=( len( Lines ), 2 ) )
        for idx in range( len( Lines ) ):
            l = Lines[idx].split( )
            stop[idx][0] = float( l[0] )
            stop[idx][1] = float( l[1] )

        return stop

    def createGraph( self, data ):
        fig = plt.figure( )
        graph = plt.plot( data[:,0], data[:,1] )
        plt.close( fig )
        return [graph[0].get_xdata( ), graph[0].get_ydata( )]

    def getValue( self, graph, value ):
        idx = (np.abs(graph[0] - value)).argmin()
        return graph[1][idx]
        

    def getCM( self, Lab ):
        CM = Lab*( self.M0 )/( self.M1 + self.M0 )
        return CM

    def getLab( self, CM ):
        Lab = CM*( self.M1 + self.M0 )/( self.M0 )
        return Lab

    def getCross( self, energy ):
        energyCM = self.getCM( energy )
        Mr = self.M0*self.M1/( self.M0 + self.M1 )
        sFactor = self.getValue( self.sFactorGraph, energyCM )
        cross = pow(10, -6)*sFactor*np.exp( -0.989534*6*np.sqrt( Mr/( energyCM/1000 ) ) )
        cross /= energyCM
        return cross

    def convertDeltaE( self, deltaE, energy ):
        return deltaE*self.getValue( self.stopGraph, energy )/self.getValue( self.stopGraph, 380 )

    def run( self, deltaE, energy ):
        nSteps = 1000
        step = deltaE/nSteps
        EStep = energy - deltaE + step/2
        deltaE = self.convertDeltaE( deltaE, energy )
        
        integral = 0
        for idx in range( nSteps ):
            stop = 0.99*self.getValue( self.stopGraph, EStep )
            stopCM = self.getCM( stop )
            cross = self.getCross( EStep )
            integral += step*cross/stopCM
            EStep += step

        self.Yield = integral
        
