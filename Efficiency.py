import numpy as np

from numpy import exp

class Efficiency( ):
    def __init__( self ):
        self.params = self.readParams( )
        
    def A( self, pos, energy ):
        return ( 1 - exp( - ( pos + self.params['d0'])/( self.params['d1'] + self.params['d2']*np.sqrt( energy ) ) ) )/( ( pos + self.params['d0'] )**2 )

    def effPeak( self, pos, energy ):
        return self.A( pos, energy )*exp( ( self.params['a'] + self.params['b']*np.log( energy ) + self.params['c']*( np.log( energy )**2 ) ) )

    def effTot( self, pos, energy ):
        return self.effPeak( pos, energy )*exp( self.params['k1'] + self.params['k2']*np.log( energy ) + self.params['k3']*( np.log( energy ) )**2 )

    def readParams( self ):
        dataDir = "data/eff_params.dat"
        fIn = open( dataDir, "r" )
        Lines = fIn.readlines( )
        params = { }
        for line in Lines:
            l = line.split( )
            params[l[0]] = float( l[1] )
        fIn.close( )
        return params

    def run( self, pos ):
        self.data = np.zeros( shape=( 100, 2 ) )
        for idx in range( 100 ):
            self.data[idx][0] = idx*0.1 + 0.1
            self.data[idx][1] = self.effPeak( pos, idx*0.1 + 0.1 )
            
