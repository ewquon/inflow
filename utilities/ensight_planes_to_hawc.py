#!/usr/bin/env python
#
# Script to read "array" sampling planes from SOWFA, output in Ensight format
# and convert to a binary HAWC-style full-field file
#
# Written by Eliot Quon (eliot.quon@nrel.gov)
#
# USAGE: Process ./*/inflowPlane_03km_U.000.U for a reference velocity of 8 m/s
#   ensight_planes_to_hawc.py 'inflowPlane_03km_U' 8.0
#
from __future__ import print_function
import numpy as np

from NWTC.datatools.dataloaders import foam_ensight_array
from FAST.InflowWind import input_template
from datatools.binario import binaryfile


def generate_inflow(prefix,uref,zref=90.0,
        ufile='u.bin',vfile='v.bin',wfile='w.bin',
        inflowfile='InflowWind_from_SOWFA.dat'):
    """Writes out one binary file for each wind component in the HAWC
    format as described in the InflowWind manual, in addition to an
    InflowWind input file"""

    inflow = foam_ensight_array('.', prefix=prefix,
                                npzdata=prefix+'.npz') # auto-detect NX,NY,NZ

    t = np.array(inflow.ts.outputTimes) # detected from time directory names
    X,Y,Z,U = inflow.sliceI() # return arrays with shape (NY,NZ) or in the case of U: (Ntimes,NY,NZ,3)
    assert(np.min(X) == np.max(X)) # assume flow is in x

    nx = inflow.ts.Ntimes
    ny = inflow.NY
    nz = inflow.NZ
    y = Y[:,0]
    z = Z[0,:]
    print('x :',nx,(t-t[0])*uref)
    print('y :',ny,y)
    print('z :',nz,z)
    dx = uref*(t[1]-t[0])
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    U[:,:,:,0] -= uref # InflowWind will add this back to the x-component
    with binaryfile(ufile,'w') as f:
        #for i in range(nx)[::-1]: # indexing goes backwards (frozen turbulence)
        for i in range(nx): # last time plane is first time snapshot
            for j in range(ny)[::-1]: # backwards
                f.write_float(U[i,j,:,0]) # forward
    print('Wrote binary',ufile)

    with binaryfile(vfile,'w') as f:
        #for i in range(nx)[::-1]: # indexing goes backwards (frozen turbulence)
        for i in range(nx): # last time plane is first time snapshot
            for j in range(ny)[::-1]: # backwards
                f.write_float(U[i,j,:,1]) # forward
    print('Wrote binary',vfile)

    with binaryfile(wfile,'w') as f:
        #for i in range(nx)[::-1]: # indexing goes backwards (frozen turbulence)
        for i in range(nx): # last time plane is first time snapshot
            for j in range(ny)[::-1]: # backwards
                f.write_float(U[i,j,:,2]) # forward
    print('Wrote binary',wfile)

    with open(inflowfile,'w') as f:
        f.write(
            input_template.format(
                WindType=5,
                RefHt=zref,
                URef=uref,
                hawc_ufile=ufile,
                hawc_vfile=vfile,
                hawc_wfile=wfile,
                nx=nx, ny=ny, nz=nz,
                dx=dx, dy=dy, dz=dz
        ))
    print('Wrote',inflowfile)


#===============================================================================
if __name__ == '__main__':
    import sys
    prefix = sys.argv[1]
    Uref = float(sys.argv[2])
    generate_inflow(prefix,Uref)

