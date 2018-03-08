#!/usr/bin/env python
#
# Generalized inflow preprocessing module for numerical ABL inflow experiments
# written by Eliot Quon (eliot.quon.nrel.gov) - 2017-07-10
#
from __future__ import print_function
import sys,os
import time

import numpy as np

import inflow.time_varying_mapped
from datatools.vtkTools import vtk_write_structured_points
import datatools.SOWFA.timeVaryingMappedBC as bc

class InflowPlane(object):

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
        * Potential temperature field: T (with shape==(Ntimes,NY,NZ))
        * Scaling function: scaling (shape==(3,NZ))

        Optionally, the following parameters may be set:
        * Reference velocity: Umean
        """
        self.verbose = verbose
        self.Umean = None # reference velocity
        self.inletMean = None # a time_varying_mapped object, which may be a
                              # profile or (time-varying) inflow plane
        self.haveField = False # True after the velocity field has been read

        self.inflowSourceDir = None
        self.pointsFile = None
        self.needUpdateMean = False # set to true update the mean inflow at every time step. 
        self.timeseries = None # used if needUpdateMean is True
        self.Useries = None # used if needUpdateMean is True
        self.Tseries = None # used if needUpdateMean is True

        # set by calculateRMS
        self.uu_mean = None
        self.vv_mean = None
        self.ww_mean = None

        # for backwards compatibility (set by readAllProfiles or readVarianceProfile)
        self.z_profile = None
        self.uu_profile = None
        self.vv_profile = None
        self.ww_profile = None


    def readField(self):
        """Stub to read inflow field"""
        print('This function should be overridden for each inflow class...')
        print('No inflow data were read.')


    def createEmptyField(self, Ly, Lz, Ny, Nz, times=[0,1000.0,2000.0]):
        """Create field with no fluctuations, for development and
        testing (e.g., to create a smooth inflow)
        """
        self.N = 3
        self.NY = Ny
        self.NZ = Nz

        self.t = times
        self.y = np.linspace(0, Ly, Ny)
        self.z = np.linspace(0, Lz, Nz)

        self.dt = self.t[1] - self.t[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        self.U = np.zeros((3,self.N,self.NY,self.NZ))
        self.T = np.zeros((self.N,self.NY,self.NZ))
        self.scaling = np.ones((3,self.NZ))

        self.haveField = True
    

    def calculateRMS(self,output=None):
        """Calculate root-mean square or standard deviation of the
        fluctuating velocities.  Output is the square root of the
        average of the fluctuations, i.e. the root-mean-square or
        standard deviation, which should match the output in PREFIX.sum.
        """
        self.uu = self.U[0,:,:,:]**2
        self.vv = self.U[1,:,:,:]**2
        self.ww = self.U[2,:,:,:]**2
        self.uu_tavg = np.mean(self.uu,0) # time averages
        self.vv_tavg = np.mean(self.vv,0)
        self.ww_tavg = np.mean(self.ww,0)
        self.uu_mean = np.mean( self.uu_tavg ) # space/time average
        self.vv_mean = np.mean( self.vv_tavg )
        self.ww_mean = np.mean( self.ww_tavg )

        print('Spatial average of <u\'u\'>, <v\'v\'>, <w\'w\'> :',self.uu_mean,self.vv_mean,self.ww_mean)

        if output is not None:
            with open(output,'w') as f:
                f.write('Spatial average of <u\'u\'>, <v\'v\'>, <w\'w\'> : {} {} {}\n'.format(self.uu_mean,self.vv_mean,self.ww_mean))
                f.write('\n   Height   Standard deviation at grid points for the u component:\n')
                for i,zi in enumerate(self.z):
                        f.write('z= {:.1f} : {}\n'.format(zi,np.sqrt(self.uu_tavg[:,i])))
                f.write('\n   Height   Standard deviation at grid points for the v component:\n')
                for i,zi in enumerate(self.z):
                        f.write('z= {:.1f} : {}\n'.format(zi,np.sqrt(self.vv_tavg[:,i])))
                f.write('\n   Height   Standard deviation at grid points for the w component:\n')
                for i,zi in enumerate(self.z):
                        f.write('z= {:.1f} : {}\n'.format(zi,np.sqrt(self.ww_tavg[:,i])))
            print('Wrote out',output)


    #===========================================================================
    #
    # Domain manipulation
    #
    #===========================================================================

    def tileY(self,ntiles,mirror=False):
        """Duplicate field in lateral direction
        'ntiles' is the final number of panels including the original

        Set 'mirror' to True to flip every other tile
        """
        ntiles = int(ntiles)
        print('Creating',ntiles,'horizontal tiles')
        print('  before:',self.U.shape)
        if mirror:
            # [0 1 2] --> [0 1 2 1 0 1 2 .. ]
            NYnew = (self.NY-1)*ntiles + 1
            Unew = np.zeros((3,self.N,NYnew,self.NZ))
            Tnew = np.zeros((  self.N,NYnew,self.NZ))
            Unew[:,:,:self.NY,:] = self.U[:,:,:self.NY,:]
            Tnew[  :,:self.NY,:] = self.T[  :,:self.NY,:]
            delta = self.NY - 1
            flipped = True
            for i in range(1,ntiles):
                if flipped:
                    Unew[:,:,i*delta+1:(i+1)*delta+1,:] = self.U[:,:,delta-1::-1,:]
                    Tnew[  :,i*delta+1:(i+1)*delta+1,:] = self.T[  :,delta-1::-1,:]
                else:
                    Unew[:,:,i*delta+1:(i+1)*delta+1,:] = self.U[:,:,1:,:]
                    Tnew[  :,i*delta+1:(i+1)*delta+1,:] = self.T[  :,1:,:]
                flipped = not flipped
            self.U = Unew
            self.T = Tnew
        else:
            # [0 1 2] --> [0 1 0 1 .. 0 1 2]
            self.U = np.tile(self.U[:,:,:-1,:],(1,1,ntiles,1))
            self.T = np.tile(self.T[  :,:-1,:],(  1,ntiles,1))
            Uplane0 = np.zeros((3,self.N,1,self.NZ))
            Tplane0 = np.zeros((  self.N,1,self.NZ))
            Uplane0[:,:,0,:] = self.U[:,:,-1,:]
            Tplane0[  :,0,:] = self.T[  :,-1,:]
            self.U = np.concatenate((self.U,Uplane0),axis=1)
            self.T = np.concatenate((self.T,Tplane0),axis=1)
        print('  after :',self.U.shape)

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
            print('Specified y range', (yMin,yMax),
                    'greater than', (self.y[0],self.y[-1]))
            return

        if dryrun: sys.stdout.write('(DRY RUN) ')
        print('Resizing fluctuations field in y-dir from [',
                self.y[0],self.y[-1],'] to [',yMin,yMax,']')
        print('  before:',self.U.shape)
        
        newNY = int(np.ceil(Ly_specified/Ly * self.NY))
        Unew = self.U[:,:,:newNY,:]
        Tnew = self.T[  :,:newNY,:]
        print('  after:',Unew.shape)
        if not dryrun:
            self.U = Unew
            self.T = Tnew
            self.NY = newNY

        ynew = yMin + np.arange(newNY,dtype=self.realtype)*self.dy
        if not dryrun:
            print('Updating y coordinates')
            self.y = ynew
        else:
            print('(DRY RUN) y coordinates:',ynew)

        # flag update for mean profile
       #self.inletMean.meanFlowSet = False


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
                print('zMin not changed from',self.z[0],'to',zMin)
                return
            if zMax < self.z[-1]:
                print('zMax not changed from',self.z[-1],'to',zMax)
                return

        self.zbot = zMin

        imin = int((zMin-self.z[0])/self.dz)
        imax = int(np.ceil((zMax-self.z[0])/self.dz))
        zMin = imin*self.dz + self.z[0]
        zMax = imax*self.dz + self.z[0]
        ioff = int((self.z[0]-zMin)/self.dz)
        if dryrun: sys.stdout.write('(DRY RUN) ')
        print('Resizing fluctuations field in z-dir from [',
                self.z[0],self.z[-1],'] to [',zMin,zMax,']')
        print('  before:',self.U.shape)
        
        newNZ = imax-imin+1
        Unew = np.zeros((3,self.N,self.NY,newNZ))
        Tnew = np.zeros((  self.N,self.NY,newNZ))
        for iz in range(ioff):
            Unew[:,:,:,iz] = self.U[:,:,:,0]
            Tnew[  :,:,iz] = self.T[  :,:,0]
        if not shrink:
            Unew[:,:,:,ioff:ioff+self.NZ] = self.U
            Tnew[  :,:,ioff:ioff+self.NZ] = self.T
        else:
            iupper = np.min((ioff+self.NZ, newNZ))
            Unew[:,:,:,ioff:iupper] = self.U[:,:,:,:iupper-ioff]
            Tnew[  :,:,ioff:iupper] = self.T[  :,:,:iupper-ioff]
        print('  after:',Unew.shape)
        if not dryrun:
            self.U = Unew
            self.T = Tnew
            self.NZ = newNZ

        znew = self.zbot + np.arange(newNZ,dtype=self.realtype)*self.dz
        if not dryrun:
            print('Updating z coordinates')
            self.z = znew
        else:
            print('(DRY RUN) z coordinates:',znew)

        if not dryrun:
            print('Resetting scaling function')
            self.scaling = np.ones((3,newNZ))

        # flag update for mean profile
       #self.inletMean.meanFlowSet = False


    #===========================================================================
    #
    # 1D mean flow set up
    #
    #===========================================================================

    def readAllProfiles(self,*args,**kwargs):
        """Automatically create a new time_varying_mapped instance and call
        the appropriate input function"""
        if self.inletMean is None:
            self.inletMean = inflow.time_varying_mapped.Inlet(self.y,self.z)
        self.inletMean.readAllProfiles(*args,**kwargs)
        # for backwards compatibility:
        self.z_profile  = self.inletMean.z_profile
        self.uu_profile = self.inletMean.uu_profile
        self.vv_profile = self.inletMean.vv_profile
        self.ww_profile = self.inletMean.ww_profile
        self.uv_profile = self.inletMean.uv_profile
        self.uw_profile = self.inletMean.uw_profile
        self.vw_profile = self.inletMean.vw_profile

    def readMeanProfile(self,*args,**kwargs):
        """Automatically create a new time_varying_mapped instance and call
        the appropriate input function"""
        if self.inletMean is None:
            self.inletMean = inflow.time_varying_mapped.Inlet(self.y,self.z)
        self.inletMean.readMeanProfile(*args,**kwargs)
        # for backwards compatibility:
        self.z_profile  = self.inletMean.z_profile

    def readVarianceProfile(self,*args,**kwargs):
        """Automatically create a new time_varying_mapped instance and call
        the appropriate input function"""
        if self.inletMean is None:
            self.inletMean = inflow.time_varying_mapped.Inlet(self.y,self.z)
        self.inletMean.readVarianceProfile(*args,**kwargs)
        # for backwards compatibility:
        self.z_profile  = self.inletMean.z_profile
        self.uu_profile = self.inletMean.uu_profile
        self.vv_profile = self.inletMean.vv_profile
        self.ww_profile = self.inletMean.ww_profile

    def setMeanProfile(self,*args,**kwargs):
        """Automatically create a new time_varying_mapped instance and call
        the appropriate input function"""
        if self.inletMean is None:
            self.inletMean = inflow.time_varying_mapped.Inlet(self.y,self.z)
        self.inletMean.setMeanProfiles(*args,**kwargs)
        # for backwards compatibility:
        self.z_profile  = self.inletMean.z_profile

    def setTkeProfile(self,*args,**kwargs):
        """Automatically create a new time_varying_mapped instance and call
        the appropriate input function"""
        if self.inletMean is None:
            self.inletMean = inflow.time_varying_mapped.Inlet(self.y,self.z)
        self.inletMean.setTkeProfile(*args,**kwargs)
        # for backwards compatibility:
        self.z_profile  = self.inletMean.z_profile


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

        if evalfn: print('Using custom scaling function instead of tanh')
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
                print('Attempting to correct scaling function with fmin =',fmin)
                self.scaling = np.maximum(self.scaling,0)
                fmin = 0
            print('Updated scaling range (dir={}) : {} {}'.format(i,fmin,fmax))
        
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
            print('Wrote scaling function to',output)


    #===========================================================================
    #
    # 2D mean flow set up
    #
    #===========================================================================

    def setInflowSourceDirectory(self,dpath,tstart=None):
        """This sets inflowSourceDir (after checking if a 'points' file exists)
        which will be passed to time_varying_mapped. Specification of a source dir
        will set the 2D inflow flag.

        If tstart is not specified, then the starting time corresponding to t=0
        in the time series of the fluctuations, will be assumed to be the first
        time directory in the source directory. 
        """
        pointsFile = os.path.join(dpath,'points')
        if os.path.isfile(pointsFile):
            self.inletMean = inflow.time_varying_mapped.Inlet(self.y,self.z,dpath,tstart)
        else:
            print('Error:',pointsFile,'does not exist!')
        self.needUpdateMean = True

    def readInflowFromBC(self,itime=None,*args,**kwargs):
        """Reads inflow at time index itime from Useries and Tseries (setup after
        a time_varying_mapped object is initialized from setupInflowSourceDirectory).

        U_planar, V_planar, W_planar, and T_planar are updated using the
        datatools/SOWFA/timeVaryingMappedBC module. These data are interpolated to
        Uinlet and Tinlet at y and z.

        If a time index, itime, is not specified then all inflow planes are read
        simultaneously and interpolation is performed to get Uinlet and Tinlet at
        a particular time.
        """
        if itime is not None:
            self.inletMean.readMeanPlane(itime,*args,**kwargs)
        else:
            self.inletMean.readAllMeanPlanes(*args,**kwargs)


    #===========================================================================
    #
    # Output routines
    #
    #===========================================================================

    def write_mapped_BC(self,
                        outputdir='boundaryData',
                        time_varying_input=None,
                        bcname='west',
                        xinlet=0.0):
        """For use with OpenFOAM's timeVaryingMappedFixedValue boundary
        condition.  This will create a points file and time directories
        in 'outputdir', which should be placed in
            constant/boundaryData/<patchname>.

        time_varying_input should be a (NT, NY, NZ, 3) array which
        shoud be aligned with the loaded data in terms of (dy, dz, dt,
        and NT)
        """
        dpath = os.path.join(outputdir, bcname)
        if not os.path.isdir(dpath):
            print('Creating output dir :',dpath)
            os.makedirs(dpath)

#        if (self.inletMean is not None) and not self.inletMean.meanFlowSet:
#            print('Note: Mean profiles have not been set or read from files')
#            self.inletMean.setup() # set up inlet profile functions
#        if not self.inletMean.tkeProfileSet:
#            print('Note: Mean TKE profile has not been set')
#            self.setTkeProfile()

        # TODO: check time-varying input
        assert(time_varying_input is not None) # TODO: GENERALIZE THIS LATER
        #NT, NY, NZ, _ = time_varying_input['U'].shape
        Uinput = time_varying_input['U']
        Tinput = time_varying_input['T']
        kinput = time_varying_input['k']
        NT, NY, NZ, _ = Uinput.shape
        u = np.zeros((NY,NZ)) # working array
        v = np.zeros((NY,NZ)) # working array
        w = np.zeros((NY,NZ)) # working array
        T = np.zeros((NY,NZ)) # working array

        # write points
        fname = os.path.join(dpath,'points')
        print('Writing',fname)
        with open(fname,'w') as f:
            f.write(bc.pointsheader.format(patchName=bcname,N=NY*NZ))
            for k in range(NZ):
                for j in range(NY):
                    f.write('({:f} {:f} {:f})\n'.format(xinlet,
                                                        self.y[j],
                                                        self.z[k]))
            f.write(')\n')

        # begin time-step loop
        for itime in range(NT):
            tname = '{:f}'.format(self.realtype(itime*self.dt)).rstrip('0').rstrip('.')

            prefix = os.path.join(dpath,tname)
            if not os.path.isdir(prefix):
                os.makedirs(prefix)

            # scale fluctuations
            u[:,:] = self.U[0,itime,:NY,:NZ] # self.U.shape == (3, self.NT, self.NY, self.NZ)
            v[:,:] = self.U[1,itime,:NY,:NZ] # self.U.shape == (3, self.NT, self.NY, self.NZ)
            w[:,:] = self.U[2,itime,:NY,:NZ] # self.U.shape == (3, self.NT, self.NY, self.NZ)
            T[:,:] = self.T[itime,:NY,:NZ] # self.T.shape == (self.NT, self.NY, self.NZ)
            for iz in range(NZ): # note: u is the original size
                u[:,iz] *= self.scaling[0,iz]
                v[:,iz] *= self.scaling[1,iz]
                w[:,iz] *= self.scaling[2,iz]

            # superimpose inlet snapshot
            u[:,:] += Uinput[itime,:,:,0]
            v[:,:] += Uinput[itime,:,:,1]
            w[:,:] += Uinput[itime,:,:,2]
            T[:,:] += Tinput[itime,:,:]

            # write out U
            fname = os.path.join(prefix,'U')
            print('Writing out',fname)
            bc.write_data(fname,
                          np.stack((u.ravel(order='F'),
                                    v.ravel(order='F'),
                                    w.ravel(order='F'))),
                          patchName=bcname,
                          timeName=tname,
                          avgValue='(0 0 0)')

            # write out T
            fname = os.path.join(prefix,'T')
            print('Writing out',fname)
            bc.write_data(fname,
                          T.ravel(order='F'),
                          patchName=bcname,
                          timeName=tname,
                          avgValue='0')

            # write out k
            fname = os.path.join(prefix,'k')
            print('Writing out',fname)
            bc.write_data(fname,
                          kinput[itime,:,:].ravel(order='F'),
                          patchName=bcname,
                          timeName=tname,
                          avgValue='0')


    def writeMappedBC(self,
            outputdir='boundaryData',
            interval=1,
            Tstart=0., Tend=None,
            xinlet=0.0,
            bcname='inlet',
            LESyfac=None, LESzfac=None,
            writePoints=True, writeU=True, writeT=True, writek=True,
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
            print('Creating output dir :',outputdir)
            os.makedirs(outputdir)

        if (self.inletMean is not None) and not self.inletMean.meanFlowSet:
            print('Note: Mean profiles have not been set or read from files')
            self.inletMean.setup() # set up inlet profile functions
        if writek and not self.inletMean.tkeProfileSet:
            print('Note: Mean TKE profile has not been set')
            self.setTkeProfile()

        # figure out indexing in the case that the LES mesh has some (integer)
        #   refinement factor
        if LESyfac >= 1 and LESzfac >= 1:
            NY = int( LESyfac*(self.NY-1) ) + 1
            NZ = int( LESzfac*(self.NZ-1) ) + 1
            print('LES y resolution increased from',self.NY,'to',NY)
            print('LES z resolution increased from',self.NZ,'to',NZ)
            jidx = np.zeros(NY,dtype=np.int)
            kidx = np.zeros(NZ,dtype=np.int)
            for j in range(NY): jidx[j] = int(j/LESyfac)
            for k in range(NZ): kidx[k] = int(k/LESzfac)
            y =             np.arange(NY,dtype=self.realtype)*self.dy/LESyfac
            z = self.zbot + np.arange(NZ,dtype=self.realtype)*self.dz/LESzfac
            print('refined y range :',np.min(y),np.max(y))
            print('refined z range :',np.min(z),np.max(z))
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
        if writePoints:
            fname = outputdir + os.sep + 'points'
            print('Writing',fname)
            with open(fname,'w') as f:
                f.write(bc.pointsheader.format(patchName=bcname,N=NY*NZ))
                for k in range(NZ):
                    for j in range(NY):
                        f.write('({:f} {:f} {:f})\n'.format(xinlet,y[j],z[k]))
                f.write(')\n')

        #
        # write time dirs
        #
        if Tend is None: 
            Tend = self.t[-1]
        istart = int(self.realtype(Tstart) / self.dt)
        iend = int(self.realtype(Tend) / self.dt)
        print('Outputting time length',(iend-istart)*self.dt)

        if writeU:
            u = np.zeros((3,self.NY,self.NZ))
        if writeT:
            theta = np.zeros((self.NY,self.NZ))

        # begin time-step loop
        for i in range(istart,iend,interval):
            itime = np.mod(i-istart,self.N) # enable looping
            tname = '{:f}'.format(self.realtype(i*self.dt)).rstrip('0').rstrip('.')
            try: os.mkdir(outputdir+os.sep+tname)
            except: pass

            if self.inletMean.interpTime:
                self.inletMean.setupForTime(i*self.dt) # update Uinlet, Tinlet

            prefix = outputdir + os.sep + tname + os.sep
            if stdout=='overwrite':
                sys.stdout.write('\rWriting {}* (itime={})'.format(prefix,itime))

            # - write out U
            if writeU:
                fname = prefix + 'U'
                if not stdout=='overwrite':
                    sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                # scale fluctuations
                u[:,:,:] = self.U[:,itime,:,:]
                for iz in range(self.NZ): # note: u is the original size
                    for i in range(3):
                        u[i,:,iz] *= self.scaling[i,iz]
                u = self.inletMean.addUmean(u)

                with open(fname,'w') as f:
                    f.write(bc.dataheader.format(patchType='vector',
                                                 patchName=bcname,
                                                 timeName=tname,
                                                 avgValue='(0 0 0)',
                                                 N=NY*NZ))
                    for k in range(NZ):
                        for j in range(NY):
                            f.write('({v[0]:f} {v[1]:f} {v[2]:f})\n'.format(v=u[:,jidx[j],kidx[k]]))
                    f.write(')\n')

            # - write out T
            if writeT:
                fname = prefix + 'T'
                if not stdout=='overwrite':
                    sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))

                theta[:,:] = self.T[itime,:,:]
                theta = self.inletMean.addTmean(theta)

                with open(fname,'w') as f:
                    f.write(bc.dataheader.format(patchType='scalar',
                                                 patchName=bcname,
                                                 timeName=tname,
                                                 avgValue='0',
                                                 N=NY*NZ))
                    for k in range(NZ):
                        for j in range(NY):
                            f.write('{s:f}\n'.format(s=theta[jidx[j],kidx[k]]))
                    f.write(')\n')

            # - write out k
            if writek:
                fname = prefix + 'k'
                if not stdout=='overwrite':
                    sys.stdout.write('Writing {} (itime={})\n'.format(fname,itime))
                with open(fname,'w') as f:
                    f.write(bc.dataheader.format(patchType='scalar',
                                                 patchName=bcname,
                                                 timeName=tname,
                                                 avgValue='0',
                                                 N=NY*NZ))
                    for k in range(NZ):
                        for j in range(NY):
                            f.write('{s:f}\n'.format(s=self.inletMean.kinlet[kidx[k]]))
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
        if (self.inletMean is not None) and not self.inletMean.meanFlowSet:
            self.inletMean.setup()

        if output_time:
            itime = int(output_time / self.dt)
        if itime is None:
            print('Need to specify itime or output_time')
            return
        if stdout=='overwrite':
            sys.stdout.write('\rWriting time step {:d} :  t= {:f}'.format(
                itime,self.t[itime]))
        else: #if stdout=='verbose':
            print('Writing out VTK for time step',itime,': t=',self.t[itime])

        if self.inletMean.interpTime:
            self.inletMean.setupForTime(self.t[itime]) # update Uinlet, Tinlet

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
        if self.inletMean is not None:
            for iz in range(self.NZ):
                U[0,:,iz] += self.inletMean.Uinlet[0,:,iz]
                V[0,:,iz] += self.inletMean.Uinlet[1,:,iz]
                W[0,:,iz] += self.inletMean.Uinlet[2,:,iz]

        # write out VTK
        vtk_write_structured_points( open(fname,'wb'), #binary mode
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
            print('Creating output dir :',outputdir)
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
            print('Creating output dir :',outputdir)
            os.makedirs(outputdir)

        fname = os.path.join(outputdir,fname)
        print('Writing VTK block',fname)

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
        vtk_write_structured_points( open(fname,'wb'), #binary mode
            Nt, self.NY, self.NZ,
            [ up,vp,wp ],
            datatype=['vector'],
            dx=step*Umean*self.dt, dy=self.dy, dz=self.dz,
            dataname=['u\''],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')


    def writeVTK_yslice(self, fname, idx=0, scaled=True):
        """Write out binary VTK file with a single vector field at a
        specified vertical index.
        """
        if (self.inletMean is not None) and not self.inletMean.meanFlowSet:
            self.inletMean.setup()

        print('Writing out VTK slice',idx,'at y=',self.y[idx],'to',fname)

        # scale fluctuations
        up = np.zeros((self.N,1,self.NZ))
        wp = np.zeros((self.N,1,self.NZ))
        vp = np.zeros((self.N,1,self.NZ))
        up[:,0,:] = self.U[0,:,idx,:]
        vp[:,0,:] = self.U[1,:,idx,:]
        wp[:,0,:] = self.U[2,:,idx,:]
        if scaled:
            for iz in range(self.NZ):
                up[:,0,iz] *= self.scaling[0,iz]
                vp[:,0,iz] *= self.scaling[1,iz]
                wp[:,0,iz] *= self.scaling[2,iz]

        # calculate instantaneous velocity
        U = up.copy()
        V = vp.copy()
        W = wp.copy()
        if self.inletMean is not None:
            for iz in range(self.NZ):
                U[:,0,iz] += self.inletMean.Uinlet[0,idx,iz]
                V[:,0,iz] += self.inletMean.Uinlet[1,idx,iz]
                W[:,0,iz] += self.inletMean.Uinlet[2,idx,iz]

        # write out VTK
        vtk_write_structured_points( open(fname,'wb'), #binary mode
            self.N, 1, self.NZ,
            [ U,V,W, up,vp,wp ],
            datatype=['vector','vector'],
            dx=self.dx, dy=1, dz=self.dz,
            dataname=['U','u\''],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')


    def writeVTK_zslice(self, fname, idx=0, scaled=True):
        """Write out binary VTK file with a single vector field at a
        specified vertical index.
        """
        if (self.inletMean is not None) and not self.inletMean.meanFlowSet:
            self.inletMean.setup()

        print('Writing out VTK slice',idx,'at z=',self.z[idx],'to',fname)

        # scale fluctuations
        up = np.zeros((self.N,self.NY,1))
        wp = np.zeros((self.N,self.NY,1))
        vp = np.zeros((self.N,self.NY,1))
        up[:,:,0] = self.U[0,:,:,idx]
        vp[:,:,0] = self.U[1,:,:,idx]
        wp[:,:,0] = self.U[2,:,:,idx]
        if scaled:
            up[:,:,0] *= self.scaling[0,idx]
            vp[:,:,0] *= self.scaling[1,idx]
            wp[:,:,0] *= self.scaling[2,idx]

        # calculate instantaneous velocity
        U = up.copy()
        V = vp.copy()
        W = wp.copy()
        if self.inletMean is not None:
            for iy in range(self.NY):
                U[:,iy,0] += self.inletMean.Uinlet[0,iy,idx]
                V[:,iy,0] += self.inletMean.Uinlet[1,iy,idx]
                W[:,iy,0] += self.inletMean.Uinlet[2,iy,idx]

        # write out VTK
        vtk_write_structured_points( open(fname,'wb'), #binary mode
            self.N, self.NY, 1,
            [ U,V,W, up,vp,wp ],
            datatype=['vector','vector'],
            dx=self.dx, dy=self.dy, dz=1,
            dataname=['U','u\''],
            origin=[0.,self.y[0],self.z[0]],
            indexorder='ijk')



    #===========================================================================
    #
    # Define aliases here
    #
    writeVTKBlock = writeVTKSeriesAsBlock

    def writeVTK_xslice(self, fname, idx=0, scaled=True):
        writeVTK(self, fname, itime=idx, scaled=scaled, stdout='verbose')


