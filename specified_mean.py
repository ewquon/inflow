#!/usr/bin/env python
#
# Module for processing the mean flow to be used by the synthetic turbulence
# modules. This is the mean velocity on the inflow plane onto which velocity 
# fluctuations will be superimposed.
#
# Written by Eliot Quon (eliot.quon.nrel.gov) - 2017-10-17
#

class plane(object):
    """A general 2-D representation fo the mean flow, i.e., U(y,z) for flow in x
    """

    def __init__(self):
        self.z_profile = None
        self.U_profile = None
        self.V_profile = None
        self.W_profile = None
        self.T_profile = None
        self.haveMeanFlow = False # True after setMeanProfiles(), readMeanProfile(), or readAllProfiles() is called
        self.haveVariances = False # True  after readVarianceProfile is called

    def setMeanProfiles(self, z,
            U, V, T,
            uu=None,vv=None,ww=None):
        """Sets the mean velocity and temperature profiles (and,
        optionally, the variance profiles as well) from user-specified
        np.ndarrays.

        Calls applyInterpolatedMeanProfile() to set up interpolation
        functions.
        """
        self.z_profile = np.array(z)
        self.U_profile = np.array(U)
        self.V_profile = np.array(V)
        self.T_profile = np.array(T)
        if uu is None:
            self.uu_profile = np.zeros(len(z))
        else:
            self.uu_profile = np.array(uu)
        if vv is None:
            self.vv_profile = np.zeros(len(z))
        else:
            self.vv_profile = np.array(vv)
        if ww is None:
            self.ww_profile = np.zeros(len(z))
        else:
            self.ww_profile = np.array(ww)
        meanNZ = len(z)
        assert(len(U_profile) == meanNZ)
        assert(len(V_profile) == meanNZ)
        assert(len(T_profile) == meanNZ)
        assert(len(uu_profile) == meanNZ)
        assert(len(vv_profile) == meanNZ)
        assert(len(ww_profile) == meanNZ)

        self.meanFlowRead = True

        if self.haveField:
            self.applyInterpolatedMeanProfile()
        else:
            print 'Note: Interpolated mean profile has not been set up since'
            print '      inflow data have not been read.'



    def readAllProfiles(self,fname='averagingProfiles.csv',delim=','):
        """Read all mean profiles (calculated separately) from a file.
        Expected columns are:
            0  1  2  3  4   5  6  7  8  9 10   11  12  13  14  15  16
            z, U, V, W, T, uu,vv,ww,uv,uw,vw, R11,R22,R33,R12,R13,R23

        Calls applyInterpolatedMeanProfile() to set up interpolation
        functions.
        """
        data = np.loadtxt(fname,delimiter=delim)
        self.z_profile = np.array(data[:,0])
        self.U_profile = np.array(data[:,1])
        self.V_profile = np.array(data[:,2])
        self.T_profile = np.array(data[:,4])
        self.uu_profile = np.array(data[:,5])
        self.vv_profile = np.array(data[:,6])
        self.ww_profile = np.array(data[:,7])

        self.meanFlowRead = True
        self.variancesRead = True

        if self.haveField:
            self.applyInterpolatedMeanProfile()
        else:
            print 'Note: Interpolated mean profile has not been set up since'
            print '      inflow data have not been read.'


    def readMeanProfile(self,
            Ufile='U.dat',
            Vfile='V.dat',
            Tfile='T.dat',
            delim=None):
        """Read planar averages (postprocessed separately) from
        individual files.  These are saved into arrays for interpolation
        assuming that the heights in all files are the same.

        Calls applyInterpolatedMeanProfile() to set up interpolation
        functions.
        """
        Udata = np.loadtxt(Ufile,delimiter=delim)
        hmean = Udata[:,0]
        Umean = Udata[:,1]
        Vmean = np.loadtxt(Vfile,delimiter=delim)[:,1]
        Tmean = np.loadtxt(Tfile,delimiter=delim)[:,1]
        assert( len(hmean)==len(Umean)
            and len(Umean)==len(Vmean)
            and len(Vmean)==len(Tmean) )
        self.z_profile = np.array(hmean)
        self.U_profile = np.array(Umean)
        self.V_profile = np.array(Vmean)
        self.T_profile = np.array(Tmean)

        self.meanFlowRead = True

        if self.haveField:
            self.applyInterpolatedMeanProfile()
        else:
            print 'Note: Interpolated mean profile has not been set up since'
            print '      inflow data have not been read.'


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


    def applyMeanProfiles(self,
            Uprofile=lambda z:[0.0,0.0,0.0],
            Tprofile=lambda z:0.0
        ):
        """Sets the mean velocity and temperature profiles (which
        affects output from writeMappedBC and writeVTK).  Called by
        readMeanProfile after reading in post-processed planar averages.

        Can also be directly called with a user-specified analytical
        profile.
        """
        if self.meanFlowSet:
            print 'Note: Mean profiles have already been set and will be overwritten'

        self.Uinlet = np.zeros((self.NZ,3))
        self.Tinlet = np.zeros(self.NZ)
        for iz,z in enumerate(self.z):
            self.Uinlet[iz,:] = Uprofile(z)
            self.Tinlet[iz]   = Tprofile(z)

        if self.verbose:
            print 'Specified mean profile:  z  U  T'
            for iz,U in enumerate(self.Uinlet):
                print self.z[iz],U,self.Tinlet[iz]

        self.meanFlowSet = True


class profile(plane):
    """A 1-D representation of the mean flow, i.e., U(z) for flow in x
    This is based on the general 2-D inflow plane.
    """
    def __init__(self):
        super(self.__class__,self).__init__(*args,**kwargs)

