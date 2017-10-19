#!/usr/bin/env python
#
# Module for processing the mean flow to be used by the synthetic turbulence
# modules. This is the mean velocity on the inflow plane onto which velocity 
# fluctuations will be superimposed.
#
# Written by Eliot Quon (eliot.quon.nrel.gov) - 2017-10-17
#
import numpy as np

class InletPlane(object):
    """A general 2-D representation fo the mean flow, i.e., U(y,z) for flow in x
    """

    def __init__(self,y,z):
        """Initialize a mean inflow plane with dimensions of len(y) and len(z)
        """
        self.y = y
        self.z = z
        self.NY = len(y)
        self.NZ = len(z)

        self.z_profile = None
        self.U_profile = None
        self.V_profile = None
        self.W_profile = None
        self.T_profile = None

        self.uu_profile = None
        self.vv_profile = None
        self.ww_profile = None

        # flow input flags
        self.haveMeanFlow = False # True after setMeanProfiles(), readMeanProfile(), or readAllProfiles() is called
        self.haveVariances = False # True after readVarianceProfile() is called

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


    #def applyMeanProfiles(self,
    def setup(self,
            Uprofile=lambda z: [0.0,0.0,0.0],
            Tprofile=lambda z: 0.0
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
        for iz,z in enumerate(self.z):
            for iy,y in enumerate(self.y):
                self.Uinlet[:,iy,iz] = Uprofile(z)
                self.Tinlet[iy,iz]   = Tprofile(z)

       #if self.verbose:
       #    print 'Specified mean profile:  z  U  T'
       #    for iz,z in enumerate(self.z):
       #        print self.z[iz], self.Uinlet[:,0,iz], self.Tinlet[0,iz]

        self.meanFlowSet = True


    #def applyInterpolatedMeanProfile(self):
    def setupInterpolated(self):
        """Helper routine to calculate interpolation functions after
        mean profiles have been input.

        Sets Ufn, Vfn, and Tfn that can be called at an arbitrary
        height.
        """
        from scipy import interpolate

        self.Ufn = interpolate.interp1d(self.z_profile, self.U_profile,
                kind='linear',fill_value='extrapolate') 

        if self.V_profile is not None:
            self.Vfn = interpolate.interp1d(self.z_profile, self.V_profile,
                    kind='linear',fill_value='extrapolate') 
        else:
            self.Vfn = lambda z: 0.0

        if self.W_profile is not None:
            self.Wfn = interpolate.interp1d(self.z_profile, self.W_profile,
                    kind='linear',fill_value='extrapolate') 
        else:
            self.Wfn = lambda z: 0.0

        if self.T_profile is not None:
            self.Tfn = interpolate.interp1d(self.z_profile, self.T_profile,
                    kind='linear',fill_value='extrapolate')
        else:
            self.Tfn = lambda z: 0.0

        self.setup(
            Uprofile=lambda z: [self.Ufn(z),self.Vfn(z),self.Wfn(z)],
            Tprofile=lambda z: self.Tfn(z)
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

