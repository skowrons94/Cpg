import numpy as np

def Edge( x, b, a ):
    return 1/( 1 + np.exp( ( x - b )/a ) )
