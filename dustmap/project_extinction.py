#==========================================================================
# This program produces HEALPix maps as a function of distance to the Sun,
# of the Lallement+19 data cube (X, Y, Z)

# the cube's entries are dA0/ds where A0 is the extintion in MAG (N_HI=2 x 10^21 A[MAG]) and s is the unit length in parsec.
# It extends up to 600 pc from the Sun on the Galactic plane and 80 pc towards upwards/downwards (wrt the Galactic plane).
#
# z/r = co-latitudine = np.sin(b)
# y/x = np.tan2(l)
#==========================================================================

import os
import h5py as h5
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

L19_h5 = h5.File( 'map3D_GAIAdr2_feb2019.h5' )
print( L19_h5.keys() )
L19_group = L19_h5.get( 'stilism' )
print( L19_group.keys() )
L19_cube = L19_group.get( 'cube_datas' )

# convertion factor from extinction to HI column density
# From the literature: MNRAS 471, 3494–3528 (2017)
# The gas-to-extinction ratio and the gas distribution in the Galaxy; It is an assumption!
MAGtoCM2 = 2.47e21 #[mag/pc * cm-2 mag-1] = [cm-2 pc-1]

X = np.arange( np.shape( L19_cube )[0] )
Y = np.arange( np.shape( L19_cube )[1] )
Z = np.arange( np.shape( L19_cube )[2] )

# cube center coordinates
SUN = [ int( len( X ) / 2 ), int( len( Y ) / 2 ), int( len( Z ) / 2 ) ]

# normalize mesh to (0, 0, 0) at the center
mesh = np.meshgrid( X - SUN[0], Y - SUN[1], Z - SUN[2] )

# create cube containing the radial distances to the Sun
#unit_distance = 6000 / len(X) # not used
distance = np.zeros([ len(X), len(Y), len(Z) ])
distance = np.sqrt( mesh[0] ** 2 + mesh[1] **2 + mesh[2] **2 )# * unit_distance

# compute GLON, GLAT from voxel position wrt the Sun
lonlat = np.zeros([ 2, len(X), len(Y), len(Z) ])
lonlat[0] = np.arctan2( mesh[0], mesh[1] ) - np.pi / 2
lonlat[1] = np.pi / 2. - np.arcsin( mesh[2] / distance )
lonlat = np.nan_to_num( lonlat )

# step of distance ladder (physical scale: 800 pc / r_steps)
r_steps = 5
#r_max = np.nanmax( distance ) Not to change!!! This the number of pixels
r_max = 600
radius_pc = np.linspace( 0, r_max, r_steps )

for r in range( 1, len( radius_pc ) ) :
    print( 'ring ', r, ' / ', len(radius_pc) )

    # select optimal NSIDE based on the resolution of the r-th ring in the cube
    resolution_threshold = np.arcsin( 1. / radius_pc[r] )
    # initialise to coarse resolution
    resolution_hmap = resolution_threshold + 1
    e = 1
    while resolution_hmap > resolution_threshold :
        resolution_previous = resolution_hmap
        # increase resolution by a factor 2
        resolution_hmap = 4 * np.pi / hp.nside2npix( int( 2 ** e ) )
        e += 1
    # check if second last resolution is closer to the r-th ring one wrt the last selected
    if np.abs(resolution_threshold - resolution_previous) < np.abs(resolution_threshold - resolution_hmap) :
        e -= 1

    NSIDE =  int( 2 ** e ) ## use this with sufficient resolution for longer
                # distance
    print( 'using optimal NSIDE = ', NSIDE )
    NPIX = hp.nside2npix( NSIDE )
    pixRange = np.arange(NPIX)
    hmap = np.zeros([ NPIX ])

    # store the pixel number information for every voxel in the cube
    ipix = hp.ang2pix( NSIDE, theta = lonlat[1], phi = lonlat[0] )
    # flag the voxels in the right distance range
    mask = ( distance > radius_pc[r-1] ) * ( distance < radius_pc[r] )
    hmap[ ipix[ mask ] ] = L19_cube[ mask ]

    mask_0s = ( hmap == 0 )
    hmap[ mask_0s ] = hmap[ mask_0s ] * float('Nan')
    mask_nan = ( hmap != hmap )
    if 0 : #not ( mask_nan == 0 ).all() :
        vec = np.array( hp.pix2vec( NSIDE, pixRange[mask_nan] ) )
        # fill invalid pixels with average of the surroundings
        hmap[ mask_nan ] = np.nanmean( hmap[ hp.query_disc( NSIDE, *vec, np.deg2rad( 0.5 ) ) ], axis=0 )
        #hmap[ mask_nan ] = 0

    if True :
        out_folder = 'LB_extinction_HEALPix_up%.0f'%( r_max )
        #out_folder = 'LB_extinction_HEALPix_up%.0f_n%.0f'%( r_max, NSIDE )
        out_file = './' + out_folder + '/LB_extinction_HEALPix_{}_{}pc.dat'.format( '%.0f'%NSIDE, '%.0f'%( radius_pc[r] * 5 ) )
        if not os.path.exists( out_folder ) :
            os.system( 'mkdir ' + out_folder )
        np.savetxt( out_file, hmap )
        print( out_file )

    hp.mollview( np.log10( hmap ), cmap = 'viridis' )
    plt.savefig( out_file.replace( '.dat', '.png' ) )

#hp.mollview( hmap )
#plt.show()
L19_h5.close()
