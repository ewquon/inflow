#!/usr/bin/env python
#
# Module for loading turbsim data
#
# written by Eliot Quon (eliot.quon@nrel.gov)
#
import sys,os
import time
import numpy as np

from inflow.general import specified_profile
from binario import binaryfile

#from memory_profiler import profile #-- THIS IS SLOW
# faster to uncomment the @profile lines and then run from the command line:
#   mprof run turbsim_bts.py
#   mprof plot

class bts(specified_profile):

    extension = '.bts'

    def __init__(self, fname=None, Umean=None, verbose=False):
        """Processes binary full-field time series output from TurbSim.

        Tested with TurbSim v2.00.05c-bjj, 25-Feb-2016
        Tested with pyTurbsim, 10-07-2017
        """
        super(self.__class__,self).__init__(verbose)
        self.Umean = Umean

        if fname is not None:
            self.readField(fname)


    def readField(self,fname):
        if not fname.endswith(self.extension):
            fname = fname + self.extension
        self._readBTS(fname)
        self.haveField = True


    #@profile
    def _readBTS(self,fname):
        """ Process AeroDyn full-field files. Fluctuating velocities and
        coordinates (y & z) are calculated.

        V.shape = (3,NY,NZ,N)  # N: number of time steps
        """
        with binaryfile(fname) as f:
            #
            # read header info
            #
            if self.verbose: print 'Reading header information from',fname

            ID = f.read_int2()
            assert( ID==7 or ID==8 )
            if ID==7: filetype = 'non-periodic'
            elif ID==8: filetype = 'periodic'
            else: filetype = 'UNKNOWN'
            if self.verbose: print '  id= {:d} ({:s})'.format(ID,filetype)

            # - read resolution settings
            self.NZ = f.read_int4()
            self.NY = f.read_int4()
            self.Ntower = f.read_int4()
            if self.verbose:
                print '  NumGrid_Z,_Y=',self.NZ,self.NY
                print '  ntower=',self.Ntower
            self.N = f.read_int4()
            self.dz = f.read_float(dtype=self.realtype)
            self.dy = f.read_float(dtype=self.realtype)
            self.dt = f.read_float(dtype=self.realtype)
            self.period  = self.realtype(self.N * self.dt)
            self.Nsize = 3*self.NY*self.NZ*self.N
            if self.verbose:
                print '  nt=',self.N
                print '  (problem size: {:d} points)'.format(self.Nsize)
                print '  dz,dy=',self.dz,self.dy
                print '  TimeStep=',self.dt
                print '  Period=',self.period

            # - read reference values
            self.uhub = f.read_float(dtype=self.realtype)
            self.zhub = f.read_float(dtype=self.realtype) # NOT USED
            self.zbot = f.read_float(dtype=self.realtype)
            if self.Umean is None:
                self.Umean = self.uhub
                if self.verbose:
                    print '  Umean = uhub =',self.Umean,'(for calculating fluctuations)'
            else: # user-specified Umean
                if self.verbose:
                    print '  Umean =',self.Umean,'(for calculating fluctuations)'
                    print '  uhub=',self.uhub,' (NOT USED)'
            if self.verbose:
                print '  HubHt=',self.zhub,' (NOT USED)'
                print '  Zbottom=',self.zbot

            # - read scaling factors
            self.Vslope = np.zeros(3,dtype=self.realtype)
            self.Vintercept = np.zeros(3,dtype=self.realtype)
            for i in range(3):
                self.Vslope[i] = f.read_float(dtype=self.realtype)
                self.Vintercept[i] = f.read_float(dtype=self.realtype)
            if self.verbose:
                # output is float64 precision by default...
                #print '  Vslope=',self.Vslope
                #print '  Vintercept=',self.Vintercept
                print '  Vslope=',[self.Vslope[i] for i in range(3)]
                print '  Vintercept=',[self.Vintercept[i] for i in range(3)]

            # - read turbsim info string
            nchar = f.read_int4()
            version = f.read(N=nchar)
            if self.verbose: print version

            #
            # read normalized data
            #
            # note: need to specify Fortran-order to properly read data using np.nditer
            t0 = time.clock()
            if self.verbose: print 'Reading normalized grid data'

            self.U = np.zeros((3,self.NY,self.NZ,self.N),order='F',dtype=self.realtype)
            if self.verbose:
                print '  U size :',self.U.nbytes/1024.**2,'MB'
            for val in np.nditer(self.U, op_flags=['writeonly']):
                val[...] = f.read_int2()
            self.U = self.U.swapaxes(3,2).swapaxes(2,1) # new shape: (3,self.N,self.NY,self.NZ)

            if self.Ntower > 0:
                if self.verbose: print 'Reading normalized tower data'
                self.Utow = np.zeros((3,self.Ntower,self.N),order='F',dtype=self.realtype)
                if self.verbose: print '  Utow size :',self.Utow.nbytes/1024.**2,'MB'
                for val in np.nditer(self.Utow, op_flags=['writeonly']):
                    val[...] = f.read_int2()

            if self.verbose: print '  Read velocitiy fields in',time.clock()-t0,'s'
                            
            #
            # calculate dimensional velocity
            #
            if self.verbose: print 'Calculating velocities from normalized data'
            for i in range(3):
                self.U[i,:,:,:] -= self.Vintercept[i]
                self.U[i,:,:,:] /= self.Vslope[i]
                if self.Ntower > 0:
                    self.Utow[i,:,:] -= self.Vintercept[i]
                    self.Utow[i,:,:] /= self.Vslope[i]
            self.U[0,:,:,:] -= self.Umean # uniform inflow w/ no shear assumed

            print '  u min/max [',np.min(self.U[0,:,:,:]),np.max(self.U[0,:,:,:]),']'
            print '  v min/max [',np.min(self.U[1,:,:,:]),np.max(self.U[1,:,:,:]),']'
            print '  w min/max [',np.min(self.U[2,:,:,:]),np.max(self.U[2,:,:,:]),']'

            self.scaling = np.ones((3,self.NZ))

            #
            # calculate coordinates
            #
            if self.verbose: print 'Calculating coordinates'
            #self.y = -0.5*(self.NY-1)*self.dy + np.arange(self.NY,dtype=self.realtype)*self.dy
            self.y =             np.arange(self.NY,dtype=self.realtype)*self.dy
            self.z = self.zbot + np.arange(self.NZ,dtype=self.realtype)*self.dz
            #self.ztow = self.zbot - np.arange(self.NZ,dtype=self.realtype)*self.dz #--NOT USED

            self.t = np.arange(self.N,dtype=self.realtype)*self.dt
            if self.verbose:
                #print 'Read times',self.t
                print 'Read times [',self.t[0],self.t[1],'...',self.t[-1],']'


