#!/usr/bin/env python
#
# Module for processing the mean flow to be used by the synthetic turbulence
# modules. This is the mean velocity on the inflow plane onto which velocity 
# fluctuations will be superimposed.
#
# Written by Eliot Quon (eliot.quon.nrel.gov) - 2017-10-17
#
import os

import numpy as np

from datatools.timeseries import TimeSeries
import datatools.SOWFA.timeVaryingMappedBC as bc

class InletPlane(object):
    """A general 2-D representation fo the mean flow, i.e., U(y,z) for flow in x
    """

    def __init__(self,y=[0],z=[0],sourceDir=None,tstart=None):
        """Initialize a mean inflow plane with dimensions of len(y) and len(z)
        """
        self.y = y # these are the dimensions from the synthetic simulation
        self.z = z # these are the dimensions from the synthetic simulation
        self.NY = len(y)
        self.NZ = len(z)
        self.tstart = tstart

        # mean profile
        self.z_profile = None
        self.U_profile = None
        self.V_profile = None
        self.W_profile = None
        self.T_profile = None

        self.uu_profile = None
        self.vv_profile = None
        self.ww_profile = None

        # mean inflow plane
        self.y_planar = None # 1D array
        self.z_planar = None # 1D array
        self.U_planar = None
        self.V_planar = None
        self.W_planar = None
        self.T_planar = None

        self.inflowIs2D = False
        self.inflowSourceDir = None
        if sourceDir is not None:
            # Inflow format is assumed to be OpenFOAM's timeVaryingMapped BC;
            # a 'points' file will be read
            self.inflowSourceDir = sourceDir
            y,z = bc.readBoundaryPoints(os.path.join(sourceDir,'points'))
            self.y_planar = y
            self.z_planar = z
            if not np.all(self.y == self.y_planar):
                print 'Note: Mean inflow plane initialized with mismatched boundary points in y'
                print '  fluctuations y-grid:',self.y
                print '     mean flow y-grid:',self.y_planar
            if not np.all(self.z == self.z_planar):
                print 'Note: Mean inflow plane initialized with mismatched boundary points in z'
                print '  fluctuations z-grid:',self.z
                print '     mean flow z-grid:',self.z_planar
            self.inflowIs2D = True

            self.Useries = TimeSeries(sourceDir,filename='U')
            self.Tseries = TimeSeries(sourceDir,filename='T')
            assert(len(self.Useries) == len(self.Tseries))
            print 'U inflow series :',self.Useries
            print 'T inflow series :',self.Tseries

            self.timeseries = np.array(self.Useries.outputTimes)
            assert(np.all(self.timeseries == self.Tseries.outputTimes))
            if tstart is None:
                self.tstart = self.timeseries[0]
            self.timeseries -= self.tstart
            print 'inflow times :',self.timeseries,'( started from',self.tstart,')'

        # flow input flags
        self.haveMeanFlow = False # True after setMeanProfiles(), readMeanProfile(), readAllProfiles(), readMeanPlane(), or readAllMeanPlanes()
        self.haveVariances = False # True after readVarianceProfile() is called

        self.interpTime = False # True if readAllMeanPlanes() is called

        # inlet setup flags
        self.meanFlowSet = False # True after setup() is called
        self.tkeProfileSet = False # True after setTkeProfile() is called (not used)


    def setMeanProfiles(self, z,
            U=None, V=None, W=None, T=None,
            uu=None,vv=None,ww=None):
        """Sets the mean velocity and temperature profiles (and,
        optionally, the variance profiles as well) from user-specified
        np.ndarrays.

        Calls setupInterpolated() to set up interpolation
        functions.
        """
        self.z_profile = np.array(z)
        Nz = len(z)
        if U is None:
            self.U_profile = np.zeros(Nz)
        else:
            assert(len(U) == Nz)
            self.U_profile = np.array(U)
        if V is None:
            self.V_profile = np.zeros(Nz)
        else:
            assert(len(V) == Nz)
            self.V_profile = np.array(V)
        if W is None:
            self.W_profile = np.zeros(Nz)
        else:
            assert(len(W) == Nz)
            self.W_profile = np.array(W)
        if T is None:
            self.T_profile = np.zeros(Nz)
        else:
            assert(len(T) == Nz)
            self.T_profile = np.array(T)
        if uu is None:
            self.uu_profile = np.zeros(Nz)
        else:
            assert(len(uu) == Nz)
            self.uu_profile = np.array(uu)
        if vv is None:
            self.vv_profile = np.zeros(Nz)
        else:
            assert(len(vv) == Nz)
            self.vv_profile = np.array(vv)
        if ww is None:
            self.ww_profile = np.zeros(Nz)
        else:
            assert(len(ww) == Nz)
            self.ww_profile = np.array(ww)

        self.meanFlowRead = True
        self.variancesRead = True

        self.setupInterpolated()


    def setTkeProfile(self,k_profile=lambda z:0.0):
        """Sets the mean TKE profiles, used for output from writeMappedBC.
        Note: This is required by the timeVaryingMappedFixedValue BC, but
        isn't actually used by windPlantSolver.*

        Can also be directly called with a user-specified analytical profile.
        """
        self.kinlet = np.zeros(self.NZ)
        for iz,z in enumerate(self.z):
            self.kinlet[iz] = k_profile(z)

        #print 'Set TKE profile:  z  k'
        #for iz,k in enumerate(self.kinlet):
        #    print self.z[iz],k

        self.tkeProfileSet = True


    def readAllProfiles(self,fname='averagingProfiles.csv',delim=','):
        """Read all mean profiles (calculated separately) from a file.
        Expected columns are:
            0  1  2  3  4   5  6  7  8  9 10   11  12  13  14  15  16
            z, U, V, W, T, uu,vv,ww,uv,uw,vw, R11,R22,R33,R12,R13,R23

        Calls setupInterpolated() to set up interpolation
        functions.
        """
        data = np.loadtxt(fname,delimiter=delim)
        self.z_profile = np.array(data[:,0])
        self.U_profile = np.array(data[:,1])
        self.V_profile = np.array(data[:,2])
        self.W_profile = np.array(data[:,3])
        self.T_profile = np.array(data[:,4])
        self.uu_profile = np.array(data[:,5])
        self.vv_profile = np.array(data[:,6])
        self.ww_profile = np.array(data[:,7])

        self.meanFlowRead = True
        self.variancesRead = True

        self.setupInterpolated()


    def readMeanProfile(self,
            Ufile='U.dat',
            Vfile=None,
            Wfile=None,
            Tfile=None,
            delim=None):
        """Read planar averages (postprocessed separately) from
        individual files.  These are saved into arrays for interpolation
        assuming that the heights in all files are the same.

        Calls setupInterpolated() to set up interpolation
        functions.
        """
        Udata = np.loadtxt(Ufile,delimiter=delim)
        hmean = Udata[:,0]
        Umean = Udata[:,1]
        if Vfile is not None:
            Vmean = np.loadtxt(Vfile,delimiter=delim)[:,1]
        else:
            Vmean = np.zeros(len(hmean))
        if Wfile is not None:
            Wmean = np.loadtxt(Wfile,delimiter=delim)[:,1]
        else:
            Wmean = np.zeros(len(hmean))
        if Tfile is not None:
            Tmean = np.loadtxt(Tfile,delimiter=delim)[:,1]
        else:
            Tmean = np.zeros(len(hmean))

        assert( len(hmean)==len(Umean)
            and len(Umean)==len(Vmean)
            and len(Vmean)==len(Wmean)
            and len(Wmean)==len(Tmean) )

        self.z_profile = np.array(hmean)
        self.U_profile = np.array(Umean)
        self.V_profile = np.array(Vmean)
        self.W_profile = np.array(Wmean)
        self.T_profile = np.array(Tmean)

        self.meanFlowRead = True

        self.setupInterpolated()


    def readVarianceProfile(self,
            uufile='uu.dat',
            vvfile='vv.dat',
            wwfile='ww.dat',
            delim=None):
        """Read planar averages (postprocessed separately) from
        individual files.  These are saved into arrays for interpolation
        assuming that the heights in all files are the same.
        """
        uudata = np.loadtxt(uufile,delimiter=delim)
        hmean = uudata[:,0]
        uumean = uudata[:,1]
        vvmean = np.loadtxt(vvfile,delimiter=delim)[:,1]
        wwmean = np.loadtxt(wwfile,delimiter=delim)[:,1]

        assert( len(hmean)==len(uumean)
            and len(uumean)==len(vvmean)
            and len(vvmean)==len(wwmean) )
        if self.meanFlowRead:
            assert( np.all(np.array(hmean) == self.z_profile) )

        self.uu_profile = np.array( uumean )
        self.vv_profile = np.array( vvmean )
        self.ww_profile = np.array( wwmean )

        #from scipy import interpolate
        #self.uu_fn = interpolate.interp1d(hmean, uumean,
        #    kind='linear',fill_value='extrapolate')
        #self.vv_fn = interpolate.interp1d(hmean, vvmean,
        #    kind='linear',fill_value='extrapolate')
        #self.ww_fn = interpolate.interp1d(hmean, wwmean,
        #    kind='linear',fill_value='extrapolate')

        self.variancesRead = True

    def readMeanPlane(self,itime,readU=True,readT=True):
        """Reads the mean inflow on a plane from data in OpenFOAM's
        timeVaryingMapped boundary condition format. 
        """
        assert(self.inflowSourceDir is not None)

        NY = len(self.y_planar)
        NZ = len(self.z_planar)
        self.U_planar = np.zeros((NY,NZ))
        self.V_planar = np.zeros((NY,NZ))
        self.W_planar = np.zeros((NY,NZ))
        self.T_planar = np.zeros((NY,NZ))

        if readU:
            vdata = bc.readVectorData(self.Useries[itime],NY=NY,NZ=NZ)
            self.U_planar[:,:] = vdata[0,:,:]
            self.V_planar[:,:] = vdata[1,:,:]
            self.W_planar[:,:] = vdata[2,:,:]
        if readT:
            self.T_planar[:,:] = bc.readScalarData(self.Tseries[itime],NY=NY,NZ=NZ)

        self.setupInterpolated2D()

        self.meanFlowRead = True

    def readAllMeanPlanes(self):
        """Reads all mean inflow planes from data in OpenFOAM's
        timeVaryingMapped boundary condition format. 
        """
        assert(self.inflowSourceDir is not None)

        Ntimes = len(self.timeseries)
        NY = len(self.y_planar)
        NZ = len(self.z_planar)
        self.U_planar = np.zeros((Ntimes,NY,NZ))
        self.V_planar = np.zeros((Ntimes,NY,NZ))
        self.W_planar = np.zeros((Ntimes,NY,NZ))
        self.T_planar = np.zeros((Ntimes,NY,NZ))

        for itime,ti in enumerate(self.timeseries):
            print 'Processing mean inflow at t=',ti,':'
            vdata = bc.readVectorData(self.Useries[itime],NY=NY,NZ=NZ)
            self.U_planar[itime,:,:] = vdata[0,:,:]
            self.V_planar[itime,:,:] = vdata[1,:,:]
            self.W_planar[itime,:,:] = vdata[2,:,:]
            self.T_planar[itime,:,:] = bc.readScalarData(self.Tseries[itime],NY=NY,NZ=NZ)

        self.interpTime = True
        self.setupInterpolated3D()

        self.meanFlowRead = True


    #def applyMeanProfiles(self,
    def setup(self,
            Uinput=None, #lambda z: [0.0,0.0,0.0],
            Tinput=None  #lambda z: 0.0
        ):
        """Sets the mean velocity and temperature profiles (which
        affects output from writeMappedBC and writeVTK).  Called by
        readMeanProfile after reading in post-processed planar averages.

        Can also be directly called with a user-specified analytical
        profile.
        """
        if self.meanFlowSet:
            print 'Note: Mean profiles have already been set and will be overwritten'

        self.Uinlet = np.zeros((3,self.NY,self.NZ))
        self.Tinlet = np.zeros((self.NY,self.NZ))

        if Uinput is not None:
            if self.inflowIs2D:
                for iz,z in enumerate(self.z):
                    for iy,y in enumerate(self.y):
                        self.Uinlet[:,iy,iz] = Uinput(y,z)
            else:
                for iz,z in enumerate(self.z):
                    Uz = Uinput(z)
                    for iy,y in enumerate(self.y):
                        self.Uinlet[:,iy,iz] = Uz

        if Tinput is not None:
            if self.inflowIs2D:
                for iz,z in enumerate(self.z):
                    for iy,y in enumerate(self.y):
                        self.Tinlet[iy,iz]   = Tinput(y,z)
            else:
                for iz,z in enumerate(self.z):
                    Tz = Tinput(z)
                    for iy,y in enumerate(self.y):
                        self.Tinlet[iy,iz]   = Tz

       #if self.verbose:
       #    print 'Specified mean profile:  z  U  T'
       #    for iz,z in enumerate(self.z):
       #        print self.z[iz], self.Uinlet[:,0,iz], self.Tinlet[0,iz]

        self.meanFlowSet = True


    #def applyInterpolatedMeanProfile(self):
    def setupInterpolated(self):
        """Helper routine to calculate interpolation functions after
        mean profiles have been input.

        Sets Ufn, Vfn, Wfn, and Tfn that can be called at an arbitrary
        height.
        """
        from scipy.interpolate import interp1d

        self.Ufn = interp1d(self.z_profile, self.U_profile,
                kind='linear',fill_value='extrapolate') 

        if self.V_profile is not None:
            self.Vfn = interp1d(self.z_profile, self.V_profile,
                    kind='linear',fill_value='extrapolate') 
        else:
            self.Vfn = lambda z: 0.0

        if self.W_profile is not None:
            self.Wfn = interp1d(self.z_profile, self.W_profile,
                    kind='linear',fill_value='extrapolate') 
        else:
            self.Wfn = lambda z: 0.0

        if self.T_profile is not None:
            self.Tfn = interp1d(self.z_profile, self.T_profile,
                    kind='linear',fill_value='extrapolate')
        else:
            self.Tfn = lambda z: 0.0

        self.setup(
            Uinput=lambda z: [self.Ufn(z),self.Vfn(z),self.Wfn(z)],
            Tinput=lambda z: self.Tfn(z)
        )


    def setupInterpolated2D(self):
        """Helper routine to calculate interpolation functions after
        mean inflow planes have been read.

        Sets Ufn, Vfn, and Tfn that can be called at an arbitrary
        height.
        """
        from scipy.interpolate import interp2d

        self.Ufn = interp2d(self.y_planar, self.z_planar, self.U_planar,
                kind='linear',fill_value='extrapolate') 

        if self.V_planar is not None:
            self.Vfn = interp2d(self.y_planar, self.z_planar, self.V_planar,
                    kind='linear',fill_value='extrapolate') 
        else:
            self.Vfn = lambda y,z: 0.0

        if self.W_planar is not None:
            self.Wfn = interp2d(self.y_planar, self.z_planar, self.W_planar,
                    kind='linear',fill_value='extrapolate') 
        else:
            self.Wfn = lambda y,z: 0.0

        if self.T_planar is not None:
            self.Tfn = interp2d(self.y_planar, self.z_planar, self.T_planar,
                    kind='linear',fill_value='extrapolate')
        else:
            self.Tfn = lambda y,z: 0.0

        self.setup(
            Uinput=lambda y,z: [self.Ufn(y,z),self.Vfn(y,z),self.Wfn(y,z)],
            Tinput=lambda y,z: self.Tfn(y,z)
        )


    def setupInterpolated3D(self):
        """Helper routine to calculate interpolation functions after
        mean inflow planes have been read (2D in space + time).

        Sets Ufn, Vfn, and Tfn that can be called at an arbitrary
        height.
        """
        from scipy.interpolate import RegularGridInterpolator #trilinear w/o triangulation

        self.Ufn = RegularGridInterpolator(
                (self.timeseries, self.y_planar, self.z_planar), self.U_planar,
                method='linear', fill_value=None) 

        if self.V_planar is not None:
            self.Vfn = RegularGridInterpolator(
                    (self.timeseries, self.y_planar, self.z_planar), self.V_planar,
                    method='linear', fill_value=None) 
        else:
            self.Vfn = lambda t,y,z: 0.0

        if self.W_planar is not None:
            self.Wfn = RegularGridInterpolator(
                    (self.timeseries, self.y_planar, self.z_planar), self.W_planar,
                    method='linear', fill_value=None) 
        else:
            self.Wfn = lambda t,y,z: 0.0

        if self.T_planar is not None:
            self.Tfn = RegularGridInterpolator(
                    (self.timeseries, self.y_planar, self.z_planar), self.T_planar,
                    method='linear', fill_value=None)
        else:
            self.Tfn = lambda t,y,z: 0.0

        self.setupForTime(0)

    def setupForTime(self,t):
        # Note: RegularGridInterpolator expects input as f([t,y,z]) as opposed
        #       to f(t,y,z) and returns a list of values as opposed to a single
        #       scalar.
        assert(self.interpTime == True)
        self.setup(
            Uinput=lambda y,z: [self.Ufn([t,y,z])[0],
                                self.Vfn([t,y,z])[0],
                                self.Wfn([t,y,z])[0]],
            Tinput=lambda y,z: self.Tfn([t,y,z])[0]
        )


    def addUmean(self,u):
        #for k in range(self.NZ):
        #    for j in range(self.NY):
        #        u[:,j,k] += self.Uinlet[:,j,k]
        return self.Uinlet + u


    def addTmean(self,the):
        #for k in range(self.NZ):
        #    for j in range(self.NY):
        #        the[j,k] += self.Tinlet[j,k]
        return self.Tinlet + the


#class profile(plane):
#    """A 1-D representation of the mean flow, i.e., U(z) for flow in x
#    This is based on the general 2-D inflow plane.
#    """
#    def __init__(self):
#        super(self.__class__,self).__init__(*args,**kwargs)

