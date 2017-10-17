#!/usr/bin/env python
#
# Module for processing output from the WAsP IEC Turbulence Generator
#
import sys,os
import time
import numpy as np

from inflow.translator import InflowPlane
from binario import binaryfile

class windsimu_binary(InflowPlane):

    extension = '.bin'

    def __init__(self, fnames=None,
            inputfile=None,
            Nx=1, Ny=1, Nz=1, Ncomp=1,
            Lx=None, Ly=None, Lz=None,
            dt=None, Umean=None,
            verbose=False):
        """Processes binary output from the WaSP IEC Turbulence
        Simulator.
        """
        super(self.__class__,self).__init__(verbose)
        
        self.fnames = fnames
        self.Umean = Umean
        self.dt = dt
        if inputfile is not None:
            self.readInput(inputfile)
        else:
            self.NX = Nx
            self.NY = Ny
            self.NZ = Nz
            self.Ncomp = Ncomp
            if Lx is None: Lx = float(Nx)
            if Ly is None: Ly = float(Ny)
            if Lz is None: Lz = float(Nz)
            self.Lx = Lx
            self.Ly = Ly
            self.Lz = Lz
        self.N = self.NX # time steps equal to x planes
        self.dx = self.Lx/self.NX
        self.dy = self.Ly/self.NY
        self.dz = self.Lz/self.NZ

        iterablefnames = hasattr(self.fnames,'__iter__')
        if self.Ncomp > 1:
            assert( iterablefnames and len(self.fnames)==self.Ncomp )
        elif self.fnames is not None and not iterablefnames:
            self.fnames = [self.fnames]

        if self.verbose:
            print 'Domain extents: ',[self.Lx,self.Ly,self.Lz]
            print 'Cell spacings: ', [self.dx,self.dy,self.dz]

        if self.dt is None and self.Umean is None:
            self.dt = 1.0
            self.Umean = self.dx
        elif self.Umean is None:
            self.Umean = self.dx / self.dt
            print 'Specified dt =',self.dt
            print 'Calculated Umean =',self.Umean
        elif self.dt is None:
            self.dt = self.dx / self.Umean
            print 'Specified Umean =',self.Umean
            print 'Calculated dt =',self.dt

        self.t = np.arange(self.NX)*self.dt
        self.y = np.arange(self.NY)*self.dy
        self.z = np.arange(self.NZ)*self.dz
        if self.verbose:
            print 't range:',[np.min(self.t),np.max(self.t)]
            print 'y range:',[np.min(self.y),np.max(self.y)]
            print 'z range:',[np.min(self.z),np.max(self.z)]

        if self.fnames is not None:
            self.readField(self.fnames)


    def readInput(self,fname):
        with open(fname,'r') as f:
            def readInt(): return int(f.readline())
            def readFloat(): return float(f.readline())
            fieldDim = readInt()
            if not fieldDim==3:
                print 'fieldDim != 3 not handled'
            Ncomp = readInt()
            components = [1,2,3]
            if Ncomp == 1:
                components = [readInt()]
            elif Ncomp == 2:
                components = components.pop(readInt())
            elif Ncomp > 3:
                print 'Ncomp ==',Ncomp,'???'
            dimensions = []
            for idim in range(fieldDim):
                dimensions.append(readInt())
            extents = []
            for idim in range(fieldDim):
                extents.append(readFloat())
            turbtype = f.readline().strip()
            if turbtype == 'land':
                self.Umean = readFloat()
                self.zref = readFloat()
                self.z0 = readFloat()
            # TODO: process other turb types
            # can skip rest of file unless we need filenames
            if self.fnames is None:
                remainingLines = f.readlines()
                self.fnames = [ line.strip() for line in remainingLines[-Ncomp:] ] 
        self.Ncomp = Ncomp
        self.NX = dimensions[0]
        self.NY = dimensions[1]
        self.NZ = dimensions[2]
        self.Lx = extents[0]
        self.Ly = extents[1]
        self.Lz = extents[2]

        print 'Read input file',fname
        if self.verbose:
            print '  {:d}D field'.format(fieldDim)
            print '  velocity components:',components
            print '  domain dimensions:',dimensions
            print '  domain extents:',extents,'m'
            print '  turbulence type:',turbtype


    def readField(self,fnames):
        self.U = np.zeros((self.Ncomp,self.NX,self.NY,self.NZ))
        self.scaling = np.ones((3,self.NZ))

        for icomp,fname in enumerate(fnames):
            if not fname.endswith(self.extension):
                fname = fname + self.extension
            self._readBinary(fname, icomp)

        self.haveField = True


    def _readBinary(self,fname,icomp):
        if self.verbose:
            print 'Reading',fname
        N = self.NX * self.NY * self.NZ
        with binaryfile(fname) as f:
            data = f.read_real4(N)
            self.U[icomp,:,:,:] =  np.array(data).reshape((self.NX,self.NY,self.NZ),order='C')
        if self.verbose:
            print 'Velocity component ranges:'
            for i in range(self.Ncomp):
                print '  U'+str(i)+': ',[np.min(self.U[:,:,:,i]),np.max(self.U[:,:,:,i])]

