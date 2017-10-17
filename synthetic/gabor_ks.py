#!/usr/bin/env python
#
# Module for processing output from GaborKS
# (Kinematic Simulation using Fourier-Gabor modes developed by Aditya Ghate)
#
import sys,os
import time
import numpy as np

from inflow.translator import InflowPlane

class GaborKS(InflowPlane):

    def __init__(self, prefix=None,
            tidx=0,
            dt=None, Umean=None,
            potentialTemperature=None,
            verbose=True):
        """Processes binary output from Gabor KS.
        """
        super(self.__class__,self).__init__(verbose)
        
        fieldnames = ['uVel','vVel','wVel']
        self.Ncomp = 3
        if potentialTemperature is not None:
            self.Ncomp += 1
            print 'Note: Potential temperature is not currently handled!'
            fieldnames.append('potT')

        self.fnames = [ '{}_{}_t{:06d}.out'.format(prefix,fieldvar,tidx) for fieldvar in fieldnames ]
        self.infofile = '{}_info_t{:06d}.out'.format(prefix,tidx)
        self.Umean = Umean
        self.dt = dt

        self.readInfo(self.infofile)

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
        else:
            if self.verbose:
                print 'Specified Umean, dt =',self.Umean,self.dt

        self.t = np.arange(self.NX)*self.dt
        self.y = np.arange(self.NY)*self.dy
        self.z = np.arange(self.NZ)*self.dz
        if self.verbose:
            print 't range:',[np.min(self.t),np.max(self.t)]
            print 'y range:',[np.min(self.y),np.max(self.y)]
            print 'z range:',[np.min(self.z),np.max(self.z)]

        if self.fnames is not None:
            self.readField(self.fnames)


    def readInfo(self,fname):
        info = np.genfromtxt(fname, dtype=None)
        self.t0 = info[0]
        self.NX = int(info[1])
        self.NY = int(info[2])
        self.NZ = int(info[3])
        self.Lx = info[4]
        self.Ly = info[5]
        self.Lz = info[6]
        self.N = self.NX # time steps equal to x planes
        self.dx = self.Lx/self.NX
        self.dy = self.Ly/self.NY
        self.dz = self.Lz/self.NZ

        self.xG,self.yG,self.zG = np.meshgrid(
                np.linspace(0,self.Lx-self.dx,self.NX),
                np.linspace(0,self.Ly-self.dy,self.NY),
                np.linspace(self.dz/2,self.Lz-(self.dz/2),self.NZ),
                indexing='ij')

        print 'Read info file',fname
        if self.verbose:
            print '  domain dimensions:',[self.NX,self.NY,self.NZ]
            print '  domain extents:',[self.Lx,self.Ly,self.Lz],'m'


    def readField(self,fnames):
        self.U = np.zeros((self.Ncomp,self.NX,self.NY,self.NZ))
        self.scaling = np.ones((3,self.NZ))

        for icomp,fname in enumerate(self.fnames):
            tmpdata = np.fromfile(fname,dtype=np.dtype(np.float64),count=-1)
            self.U[icomp,:,:,:] = tmpdata.reshape((self.NX,self.NY,self.NZ),order='F')

        self.haveField = True

