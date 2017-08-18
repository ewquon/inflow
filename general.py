#!/usr/bin/env python
#
# Generalized inflow preprocessing module for numerical ABL inflow experiments
# written by Eliot Quon (eliot.quon.nrel.gov) - 2017-07-10
#
import sys,os
import numpy as np
import VTKwriter

import time

class specified_profile(object):

    realtype = np.float32


    def __init__(self, verbose=False):
        """Defaults are set here

        After initialization, the following variables should be set:
        * Dimensions: NY, NZ (horizontal, vertical)
        * Number of time snapshots: N
        * Spacings/step size: dt, dy, dz
        * Rectilinear grid: y, z
        * Sampling times: t
        * Velocity field: U (with shape==(3,Ntimes,NY,NZ))
        * Scaling function: scaling (shape==(3,NZ))

        Optionally, the following parameters may be set:
        * Reference velocity: Umean
        """
        self.verbose = verbose
        self.Umean = None

        self.haveField = False
        self.meanProfilesSet = False
        self.meanProfilesRead = False
        self.variancesRead = False
        self.tkeProfileSet = False


    def createEmptyField(self, Ly, Lz, Ny, Nz):
        """Create field with no fluctuations, for development and
        testing (e.g., to create a smooth inflow)
        """
        self.N = 3
        self.NY = Ny
        self.NZ = Nz

        self.t = [0, 1000.0, 2000.0]
        self.y = np.linspace(0, Ly, Ny)
        self.z = np.linspace(0, Lz, Nz)

        self.dt = self.t[1] - self.t[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        self.U = np.zeros((3,self.N,self.NY,self.NZ))
        self.scaling = np.ones((3,self.NZ))

        self.haveField = True
    
    def readField(self):
        """Stub to read inflow field"""
        print 'This function should be overridden for each inflow class...'
        print 'No inflow data were read.'


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

        self.meanProfilesRead = True

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

        self.meanProfilesRead = True

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

        self.meanProfilesRead = True

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
        if self.meanProfilesRead:
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


    def calculateRMS(self,output=None):
        """Calculate root-mean square or standard deviation of the
        fluctuating velocities.  Output is the square root of the
        average of the fluctuations, i.e. the root-mean-square or
        standard deviation, which should match the output in PREFIX.sum.
        """
        self.uu = self.U[0,:,:,:]**2
        self.vv = self.U[1,:,:,:]**2
        self.ww = self.U[2,:,:,:]**2
        self.uu_tavg = np.mean(self.uu,2) # time averages
        self.vv_tavg = np.mean(self.vv,2)
        self.ww_tavg = np.mean(self.ww,2)
        self.uu_mean = np.mean( self.uu_tavg ) # space/time (ensemble) average
        self.vv_mean = np.mean( self.vv_tavg )
        self.ww_mean = np.mean( self.ww_tavg )

        print 'Spatial average of <u\'u\'>, <v\'v\'>, <w\'w\'> :',self.uu_mean,self.vv_mean,self.ww_mean

        if output is not None:
            with open(output,'w') as f:
                f.write('   Height   Standard deviation at grid points for the u component:\n')
                for i,zi in enumerate(self.z):
                        f.write('z= {:.1f} : {}\n'.format(zi,np.sqrt(self.uu_tavg[:,i])))
                f.write('   Height   Standard deviation at grid points for the v component:\n')
                for i,zi in enumerate(self.z):
                        f.write('z= {:.1f} : {}\n'.format(zi,np.sqrt(self.vv_tavg[:,i])))
                f.write('   Height   Standard deviation at grid points for the w component:\n')
                for i,zi in enumerate(self.z):
                        f.write('z= {:.1f} : {}\n'.format(zi,np.sqrt(self.ww_tavg[:,i])))
            print 'Wrote out',output


    def tileY(self,ntiles,mirror=False):
        """Duplicate field in lateral direction
        'ntiles' is the final number of panels including the original

        Set 'mirror' to True to flip every other tile
        """
        ntiles = int(ntiles)
        print 'Creating',ntiles,'horizontal tiles'
        print '  before:',self.U.shape
        if mirror:
            # [0 1 2] --> [0 1 2 1 0 1 2 .. ]
            NYnew = (self.NY-1)*ntiles + 1
            Unew = np.zeros((3,self.N,NYnew,self.NZ))
            Unew[:,:,:self.NY,:] = self.U[:,:,:self.NY,:]
            delta = self.NY - 1
            flipped = True
            for i in range(1,ntiles):
                if flipped:
                    Unew[:,:,i*delta+1:(i+1)*delta+1,:] = self.U[:,:,delta-1::-1,:]
                else:
                    Unew[:,:,i*delta+1:(i+1)*delta+1,:] = self.U[:,:,1:,:]
                flipped = not flipped
            self.U = Unew
        else:
            # [0 1 2] --> [0 1 0 1 .. 0 1 2]
            self.U = np.tile(self.U[:,:,:-1,:],(1,1,ntiles,1))
            plane0 = np.zeros((3,self.N,1,self.NZ))
            plane0[:,:,0,:] = self.U[:,:,-1,:]
            self.U = np.concatenate((self.U,plane0),axis=1)
        print '  after :',self.U.shape

        self.NY = NYnew
        assert( self.U.shape == (3,self.N,self.NY,self.NZ) )
        self.y = np.arange(self.NY,dtype=self.realtype)*self.dy


    def resizeY(self,yMin=None,yMax=None,dryrun=False):
        """Resize inflow domain to fit LES boundary and update NY.
        Min(y) will be shifted to coincide with yMin.
        """
        if yMin is None:
            yMin = self.y[0]
        if yMax is None:
            yMax = self.y[-1]
        Ly_specified = yMax - yMin
        Ly = self.y[-1] - self.y[0]
        if Ly_specified > Ly:
            print 'Specified y range', (yMin,yMax), \
                    'greater than', (self.y[0],self.y[-1])
            return

        if dryrun: sys.stdout.write('(DRY RUN) ')
        print 'Resizing fluctuations field in y-dir from [', self.y[0],self.y[-1],'] to [',yMin,yMax,']'
        print '  before:',self.U.shape
        
        newNY = int(np.ceil(Ly_specified/Ly * self.NY))
        Unew = self.U[:,:,:newNY,:]
        print '  after:',Unew.shape
        if not dryrun:
            self.U = Unew
            self.NY = newNY

        ynew = yMin + np.arange(newNY,dtype=self.realtype)*self.dy
        if not dryrun:
            print 'Updating y coordinates'
            self.y = ynew
        else:
            print '(DRY RUN) y coordinates:',ynew

        # flag update for mean profile
        self.meanProfilesSet = False


    def resizeZ(self,zMin=None,zMax=None,shrink=False,dryrun=False):
        """Set/extend inflow domain to fit LES boundary and update NZ.
        Values between zMin and min(z) will be duplicated from
        V[:3,y,z=min(z),t], whereas values between max(z) and zMax will
        be set to zero.

        By default, this function will not resize inflow plane to a
        smaller domain; to override this, set shrink to True. (NOT YET
        TESTED)
        """
        if zMin is None:
            zMin = self.z[0]
        if zMax is None:
            zMax = self.z[-1]
        if not shrink:
            if zMin > self.z[0]:
                print 'zMin not changed from',self.z[0],'to',zMin
                return
            if zMax < self.z[-1]:
                print 'zMax not changed from',self.z[-1],'to',zMax
                return

        self.zbot = zMin

        imin = int(zMin/self.dz)
        imax = int(np.ceil(zMax/self.dz))
        zMin = imin*self.dz
        zMax = imax*self.dz
        ioff = int((self.z[0]-zMin)/self.dz)
        if dryrun: sys.stdout.write('(DRY RUN) ')
        print 'Resizing fluctuations field in z-dir from [', self.z[0],self.z[-1],'] to [',zMin,zMax,']'
        print '  before:',self.U.shape
        
        newNZ = imax-imin+1
        Unew = np.zeros( (3,self.N,self.NY,newNZ) )
        for iz in range(ioff):
            Unew[:,:,:,iz] = self.U[:,:,:,0]
        if not shrink:
            Unew[:,:,:,ioff:ioff+self.NZ] = self.U
        else:
            iupper = np.min(ioff+self.NZ, newNZ)
            Unew[:,:,:,ioff:iupper] = self.U[:,:,:,:iupper-ioff]
        print '  after:',Unew.shape
        if not dryrun:
            self.U = Unew
            self.NZ = newNZ

        znew = self.zbot + np.arange(newNZ,dtype=self.realtype)*self.dz
        if not dryrun:
            print 'Updating z coordinates'
            self.z = znew
        else:
            print '(DRY RUN) z coordinates:',znew

        if not dryrun:
            print 'Resetting scaling function'
            self.scaling = np.ones((3,newNZ))

        # flag update for mean profile
        self.meanProfilesSet = False


    def applyInterpolatedMeanProfile(self):
        """Helper routine to calculate interpolation functions after
        mean profiles have been input.

        Sets Ufn, Vfn, and Tfn that can be called at an arbitrary
        height.
        """
        from scipy import interpolate
        self.Ufn = interpolate.interp1d(self.z_profile, self.U_profile,
                kind='linear',fill_value='extrapolate') 
        self.Vfn = interpolate.interp1d(self.z_profile, self.V_profile,
                kind='linear',fill_value='extrapolate') 
        self.Tfn = interpolate.interp1d(self.z_profile, self.T_profile,
                kind='linear',fill_value='extrapolate')

        self.applyMeanProfiles(
            Uprofile=lambda z: [self.Ufn(z),self.Vfn(z),0.0],
            Tprofile=lambda z: self.Tfn(z)
        )


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
        if self.meanProfilesSet:
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

        self.meanProfilesSet = True


    def setTkeProfile(self,k_profile=lambda z:0.0):
        """Sets the mean TKE profiles (affects output from writeMappedBC)

        Can also be directly called with a user-specified analytical profile.
        """
        self.kinlet = np.zeros(self.NZ)
        for iz,z in enumerate(self.z):
            self.kinlet[iz] = k_profile(z)

        #print 'Set TKE profile:  z  k'
        #for iz,k in enumerate(self.kinlet):
        #    print self.z[iz],k

        self.tkeProfileSet = True


    def setScaling(self,
            tanh_z90=0.0,
            tanh_z50=0.0,
            max_scaling=1.0,
            output=''):
        """Set scaling of fluctuations with height.  The scaling
        function ranges from 0 to max_scaling.  The heights at which the
        fluctuation magnitudes are decreased by 90% and 50% (tanh_z90
        and tanh_z50, respectively) are specified to scale the
        hyperbolic tangent function; tanh_z90 should be set to
        approximately the inversion height:
            f = max_scaling * 0.5( tanh( k(z-z_50) ) + 1 )
        where
            k = arctanh(0.8) / (z_90-z_50)
        Note: If extendZ is used, that should be called to update the z
        coordinates prior to using this routine.

        max_scaling may be:
        1) a constant, equal for the x, y, and z directions; 
        2) a list or nd.array of scalars; or
        3) a list of lambda functions for non-tanh scaling.

        Note: The scaled perturbations is no longer conservative, i.e., the
        field is not divergence free. The cells adjacent to the inflow boundary
        will make the inflow field solenoidal after the first pressure solve.
        """
        evalfn = False
        if isinstance(max_scaling,(list,tuple,np.ndarray)):
            assert( len(max_scaling) == 3 )
            if any( [ hasattr(f, '__call__') for f in max_scaling ] ): evalfn = True
        else:
            if hasattr(max_scaling,'__call__'): evalfn = True
            max_scaling = [max_scaling,max_scaling,max_scaling]

        if evalfn: print 'Using custom scaling function instead of tanh'
        else:
            assert( tanh_z90 > 0 and tanh_z50 > 0 )
            k = np.arctanh(0.8) / (tanh_z90-tanh_z50)

        for i in range(3):
            if evalfn:
                self.scaling[i,:] = max_scaling[i](self.z)
            else:
                self.scaling[i,:] = max_scaling[i] * 0.5*(np.tanh(-k*(self.z-tanh_z50)) + 1.0)
            fmin = np.min(self.scaling[i,:])
            fmax = np.max(self.scaling[i,:])
            #assert( fmin >= 0. and fmax <= max_scaling[i] )
            assert(fmax <= max_scaling[i])
            if fmin < 0:
                print 'Attempting to correct scaling function with fmin =',fmin
                self.scaling = np.maximum(self.scaling,0)
                fmin = 0
            print 'Updated scaling range (dir={}) : {} {}'.format(i,fmin,fmax)
        
        if output:
            with open(output,'w') as f:
                if evalfn:
                    f.write('# custom scaling function\n')
                else:
                    f.write('# tanh scaling parameters: z_90={:f}, z_50={:f}, max_scaling={}\n'.format(
                        tanh_z90,tanh_z50,max_scaling))
                f.write('# z  f_u(z)  f_v(z)  f_w(z)\n')
                for iz,z in enumerate(self.z):
                    f.write(' {:f} {f[0]:g} {f[1]:g} {f[2]:g}\n'.format(z,f=self.scaling[:,iz]))
            print 'Wrote scaling function to',output


    def writeMappedBC(self,
            outputdir='boundaryData',
            interval=1,
            Tstart=0., Tend=None,
            xinlet=0.0,
            bcname='inlet',
            LESyfac=None, LESzfac=None,
            writeU=True, writeT=True, writek=True,
            stdout='verbose'):
        """For use with OpenFOAM's timeVaryingMappedFixedValue boundary
        condition.  This will create a points file and time directories
        in 'outputdir', which should be placed in
            constant/boundaryData/<patchname>.

        The output interval is in multiples of the inflow time step,
        with output up to time Tend.  Inlet location and bcname should
        correspond to the LES inflow plane.
        
        LESyfac and LESzfac specify refinement in the microscale domain.
        """
        if not os.path.isdir(outputdir):
            print 'Creating output dir :',outputdir
            os.makedirs(outputdir)

        if not self.meanProfilesSet:
            print 'Note: Mean profiles have not been set or read from files'
            self.applyMeanProfiles() # set up inlet profile functions
        if writek and not self.meanProfilesSet:
            print 'Note: Mean TKE profile has not been set'
            self.setTkeProfile()

        if writeU:
            # for scaling fluctuations
            up = np.zeros((3,1,self.NY,self.NZ))

        if LESyfac >= 1 and LESzfac >= 1:
            NY = int( LESyfac*(self.NY-1) ) + 1
            NZ = int( LESzfac*(self.NZ-1) ) + 1
            print 'LES y resolution increased from',self.NY,'to',NY
            print 'LES z resolution increased from',self.NZ,'to',NZ
            jidx = np.zeros(NY,dtype=np.int)
            kidx = np.zeros(NZ,dtype=np.int)
            for j in range(NY): jidx[j] = int(j/LESyfac)
            for k in range(NZ): kidx[k] = int(k/LESzfac)
            y =             np.arange(NY,dtype=self.realtype)*self.dy/LESyfac
            z = self.zbot + np.arange(NZ,dtype=self.realtype)*self.dz/LESzfac
            print 'refined y range :',np.min(y),np.max(y)
            print 'refined z range :',np.min(z),np.max(z)
        else:
            NY = self.NY
            NZ = self.NZ
            jidx = np.arange(NY,dtype=np.int)
            kidx = np.arange(NZ,dtype=np.int)
            y = self.y
            z = self.z

        #
        # write points file
        #
        pointshdr = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant/boundaryData/{patchName}";
    object      points;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"""
        fname = outputdir + os.sep + 'points'
        print 'Writing',fname
        with open(fname,'w') as f:
            f.write(pointshdr.format(patchName=bcname))
            f.write('{:d}\n(\n'.format(NY*NZ))
            for k in range(NZ):
                for j in range(NY):
                    f.write('({:f} {:f} {:f})\n'.format(xinlet,y[j],z[k]))
            f.write(')\n')

        #
        # write time dirs
        #
        datahdr = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.x                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {patchType}AverageField;
    location    "constant/boundaryData/{patchName}/{timeName}";
    object      values;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Average
{avgValue}\n\n"""
        if Tend is None: 
            Tend = self.t[-1]
        istart = int(self.realtype(Tstart) / self.dt)
        iend = int(self.realtype(Tend) / self.dt)
        print 'Outputting time length',(iend-istart)*self.dt
        #
        # time-step loop
        #
        for i in range(istart,iend,interval):
            itime = np.mod(i-istart,self.N)
            tname = '{:f}'.format(self.realtype(i*self.dt)).rstrip('0').rstrip('.')
            try: os.mkdir(outputdir+os.sep+tname)
            except: pass

            prefix = outputdir + os.sep + tname + os.sep
            if stdout=='overwrite':
                sys.stdout.write('\rWriting {}* (itime={})'.format(prefix,itime))

            #
            # write out U
            #
            if writeU:
                fname = prefix + 'U'
                if not stdout=='overwrite':
                    sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                # scale fluctuations
                up[0,0,:,:] = self.U[0,itime,:,:]
                up[1,0,:,:] = self.U[1,itime,:,:]
                up[2,0,:,:] = self.U[2,itime,:,:]
                for iz in range(self.NZ): # note: up is the original size
                    for i in range(3):
                        up[i,0,:,iz] *= self.scaling[i,iz]
                with open(fname,'w') as f:
                    f.write(datahdr.format(patchType='vector',patchName=bcname,timeName=tname,avgValue='(0 0 0)'))
                    f.write('{:d}\n(\n'.format(NY*NZ))
                    for k in range(NZ):
                        for j in range(NY):
                            f.write('({v[0]:f} {v[1]:f} {v[2]:f})\n'.format(
                                v=self.Uinlet[kidx[k],:] + up[:,0,jidx[j],kidx[k]]
                            ))
                    f.write(')\n')

            #
            # write out T
            #
            if writeT:
                fname = prefix + 'T'
                if not stdout=='overwrite':
                    sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                with open(fname,'w') as f:
                    f.write(datahdr.format(patchType='scalar',patchName=bcname,timeName=tname,avgValue='0'))
                    f.write('{:d}\n(\n'.format(NY*NZ))
                    for k in range(NZ):
                        for j in range(NY):
                            f.write('{s:f}\n'.format(s=self.Tinlet[kidx[k]]))
                    f.write(')\n')

            #
            # write out k
            #
            if writek:
                fname = prefix + 'k'
                if not stdout=='overwrite':
                    sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                with open(fname,'w') as f:
                    f.write(datahdr.format(patchType='scalar',patchName=bcname,timeName=tname,avgValue='0'))
                    f.write('{:d}\n(\n'.format(NY*NZ))
                    for k in range(NZ):
                        for j in range(NY):
                            f.write('{s:f}\n'.format(s=self.kinlet[kidx[k]]))
                    f.write(')\n')

        # end of time-step loop
        if stdout=='overwrite':
            sys.stdout.write('\n')


    def writeVTK(self, fname,
            itime=None,
            output_time=None,
            scaled=True,
            stdout='overwrite'):
        """Write out binary VTK file with a single vector field for a
        specified time index or output time.
        """
        if not self.meanProfilesSet: self.applyMeanProfiles()

        if output_time:
            itime = int(output_time / self.dt)
        if itime is None:
            print 'Need to specify itime or output_time'
            return
	if stdout=='overwrite':
            sys.stdout.write('\rWriting time step {:d} :  t= {:f}'.format(
                itime,self.t[itime]))
	else: #if stdout=='verbose':
            print 'Writing out VTK for time step',itime,': t=',self.t[itime]

        # scale fluctuations
        up = np.zeros((1,self.NY,self.NZ))
        wp = np.zeros((1,self.NY,self.NZ))
        vp = np.zeros((1,self.NY,self.NZ))
        up[0,:,:] = self.U[0,itime,:,:]
        vp[0,:,:] = self.U[1,itime,:,:]
        wp[0,:,:] = self.U[2,itime,:,:]
        if scaled:
            for iz in range(self.NZ):
                up[0,:,iz] *= self.scaling[0,iz]
                vp[0,:,iz] *= self.scaling[1,iz]
                wp[0,:,iz] *= self.scaling[2,iz]

        # calculate instantaneous velocity
        U = up.copy()
        V = vp.copy()
        W = wp.copy()
        for iz in range(self.NZ):
            U[0,:,iz] += self.Uinlet[iz,0]
            V[0,:,iz] += self.Uinlet[iz,1]
            W[0,:,iz] += self.Uinlet[iz,2]

        # write out VTK
        VTKwriter.vtk_write_structured_points( open(fname,'wb'), #binary mode
            1, self.NY, self.NZ,
            [ U,V,W, up,vp,wp ],
            datatype=['vector','vector'],
            dx=1.0, dy=self.dy, dz=self.dz,
            dataname=['U','u\''], #['fluctuations'], #dataname=['inflow_velocity'],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')


    def writeVTKSeries(self,
            outputdir='.',
            prefix='inflow',
            step=1,
            scaled=True,
            stdout='overwrite'):
        """Driver for writeVTK to output a range of times"""
        if not os.path.isdir(outputdir):
            print 'Creating output dir :',outputdir
            os.makedirs(outputdir)

        for i in range(0,self.N,step):
            fname = outputdir + os.sep + prefix + '_' + str(i) + '.vtk'
            self.writeVTK(fname,itime=i,scaled=scaled,stdout=stdout)
	if stdout=='overwrite': sys.stdout.write('\n')


    def writeVTKSeriesAsBlock(self,
            fname='frozen_block.vtk',
            outputdir=None,
            step=1,
            scaled=True):
        """Write out a 3D block wherein the x planes are comprised of
        temporal snapshots spaced (Umean * step * dt) apart.

        This invokes Taylor's frozen turbulence assumption.
        """
        if outputdir is None:
            outputdir = '.'
        elif not os.path.isdir(outputdir):
            print 'Creating output dir :',outputdir
            os.makedirs(outputdir)

        fname = os.path.join(outputdir,fname)
        print 'Writing VTK block',fname

        if self.Umean is not None:
            Umean = self.Umean
        else:
            Umean = 1.0

        # scale fluctuations
        Nt = self.N / step
        up = np.zeros((Nt,self.NY,self.NZ))
        vp = np.zeros((Nt,self.NY,self.NZ))
        wp = np.zeros((Nt,self.NY,self.NZ))
        up[:,:,:] = self.U[0,:Nt*step:step,:,:]
        vp[:,:,:] = self.U[1,:Nt*step:step,:,:]
        wp[:,:,:] = self.U[2,:Nt*step:step,:,:]
        if scaled:
            for iz in range(self.NZ):
                up[:,:,iz] *= self.scaling[0,iz]
                vp[:,:,iz] *= self.scaling[1,iz]
                wp[:,:,iz] *= self.scaling[2,iz]

        # write out VTK
        VTKwriter.vtk_write_structured_points( open(fname,'wb'), #binary mode
            Nt, self.NY, self.NZ,
            [ up,vp,wp ],
            datatype=['vector'],
            dx=step*Umean*self.dt, dy=self.dy, dz=self.dz,
            dataname=['u\''],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')


    #****************#
    # define aliases #
    #****************#
    writeVTKBlock = writeVTKSeriesAsBlock

