#!/usr/bin/env python
import sys,os
import time
import numpy as np

from inflow import base
from binario import binaryfile

class binary(base.specified_profile):

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
        if inputfile is not None:
            self.readInput(inputfile)
        else:
            self.Nx = Nx
            self.Ny = Ny
            self.Nz = Nz
            self.Ncomp = Ncomp
            if Lx is None: Lx = float(Nx)
            if Ly is None: Ly = float(Ny)
            if Lz is None: Lz = float(Nz)
            self.Lx = Lx
            self.Ly = Ly
            self.Lz = Lz
        self.dx = self.Lx/self.Nx
        self.dy = self.Ly/self.Ny
        self.dz = self.Lz/self.Nz

        iterablefnames = hasattr(self.fnames,'__iter__')
        if self.Ncomp > 1:
            assert( iterablefnames and len(self.fnames)==self.Ncomp )
        elif self.fnames is not None and not iterablefnames:
            self.fnames = [self.fnames]

        if self.verbose:
            print 'Domain extents: ',[self.Lx,self.Ly,self.Lz]
            print 'Cell spacings: ', [self.dx,self.dy,self.dz]

        if dt is None and Umean is None:
            dt = 1.0
            Umean = self.dx
        elif Umean is None:
            self.dt = dt
            self.Umean = self.dx / self.dt
            print 'Specified dt =',self.dt
            print 'Calculated Umean =',self.Umean
        elif dt is None:
            self.Umean = Umean
            self.dt = self.dx / Umean
            print 'Specified Umean =',self.Umean
            print 'Calculated dt =',self.dt

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
            turbtype = f.readline()
            # can skip rest of file unless we need filenames
            if self.fnames is None:
                remainingLines = f.readlines()
                self.fnames = [ line.strip() for line in remainingLines[-Ncomp:] ] 
        self.Ncomp = Ncomp
        self.Nx = dimensions[0]
        self.Ny = dimensions[1]
        self.Nz = dimensions[2]
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
        self.V = np.zeros((self.Ncomp,self.Nx,self.Ny,self.Nz))
        for icomp,fname in enumerate(fnames):
            if not fname.endswith(self.extension):
                fname = fname + self.extension
            self._readBinary(fname, icomp)
        self.haveField = True


    def _readBinary(self,fname,icomp):
        if self.verbose:
            print 'Reading',fname
        N = self.Nx * self.Ny * self.Nz
        with binaryfile(fname) as f:
            data = f.read_real4(N)
            self.V[icomp,:,:,:] =  np.array(data).reshape((self.Nx,self.Ny,self.Nz),order='C')
        if self.verbose:
            print 'Velocity component ranges:'
            for i in range(self.Ncomp):
                print '  V'+str(i)+': ',[np.min(self.V[:,:,:,i]),np.max(self.V[:,:,:,i])]

