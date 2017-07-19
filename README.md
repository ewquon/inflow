Inflow Processing Module
========================

This Python module should be used to prepare a time-varying Dirichlet boundary
for a high-fidelity CFD analysis (e.g., a microscale large-eddy simulation using
SOWFA). 


Usage
-----

For high-fidelity wind plant flow simulations--typically a large eddy simulation
(LES)--a time-varying velocity field is needed to define the inflow boundary.
This can come from:

* A precursor LES,
* A recycling region within the same computational domain,
* Mesoscale model output, or
* Field measurements.

In the former two cases, the turbulence is fully resolved by the LES. However,
when limited space-time data area available as is the case for thee latter two
cases, the turbulence within the LES resolved scales is not known a priori. In
this situation, it is possible to superimpose velocity fluctuations from a
synthetic turbulence model--two popular choices are Turbsim and the Mann model. 

This sample code will take a smooth (mean) inflow profile and a desired 
variance profile, and add synthetic velocity fluctuations coming from a Turbsim
simulation of a 1-km wide by 0.8-km high inflow plane and generate inflow for a
3-km by 1-km SOWFA inflow boundary; the time-varying inflow field is provided
in the Turbsim binary `bts` format.

```python
import numpy as np
from inflow.synthetic.turbsim import bts

precDir = '/Users/equon/MMC/inflowStudy/precursors/neutral/8mps_shear0.2_TI10.0_rerun2/postProcessing'

# Read plane of fluctuating velocities to be superimposed
inflowPlane = bts('turbsimBinaryData',verbose=True)

outputInterval = int(round(1.0/inflowPlane.dt))

# Calculate mean variance (and optionally output standard deviations)
inflowPlane.calculateRMS(output='rms.txt')

meanVariances = np.array([ inflowPlane.uu_mean,
                           inflowPlane.vv_mean,
                           inflowPlane.ww_mean ]) # output from inflow model

# Write out VTK set for visualizing the Turbsim solution (OPTIONAL)
#inflowPlane.writeVTKSeries(outputdir='vtkSeries',prefix='turbsim_output',step=outputInterval)

# Resize the domain (OPTIONAL)
# Note: This should be done prior to calling readMeanProfile so
#       that the mean profile is set on the new resized domain.
#inflowPlane.tileY(Npanels=3, mirror=True)
inflowPlane.resizeY(yMin=0.0, yMax=3000.0)
inflowPlane.resizeZ(zMin=0.0, zMax=1000.0)

# Set input profiles (e.g., from mesoscale model or precursor planar average)
# Note 1: This sets the *_profile properties of inflowPlane 
# Note 2: Can also call readProfile('averagingProfiles.csv',delim=',')
#         or setMeanProfile(U,V,T) and setVarianceProfile(uu,vv,ww)
#inflowPlane.readMeanProfile('U.dat','V.dat','T.dat')
#inflowPlane.readVarianceProfile('uu.dat','vv.dat','ww.dat') # OPTIONAL
#inflowPlane.setMeanProfile(U,V,T)
#inflowPlane.setMeanProfile(U,V,T,uu,vv,ww) # OPTIONAL
inflowPlane.readAllProfiles(precDir+'/averagingProfiles.csv')

# Set reference value for scaling (OPTIONAL--this data may not be available)
maxVarianceProfile = np.array([ np.max(inflowPlane.uu_profile),
                                np.max(inflowPlane.vv_profile),
                                np.max(inflowPlane.ww_profile) ])
print 'Reference TIx, TIy, TIz =',np.sqrt(maxVarianceProfile),'m/s'

# Set up the TKE profile
# Note: This is written out but not directly used by ABLSolver, windPlantSolver, etc
inflowPlane.setTkeProfile(lambda z: 0.1) # constant field

# Since the synthetic velocity fluctuations have no height dependence,
# scale the velocities so that the fluctuations decay above some height, zi
# (where zi == inversion height >= the ABL height)
# Note 1: Since the scaling function is applied to the velocities, proper
#          scaling of variances requires the square root of the scaling factor
# Note 2: 'max_scaling' can be a constant, a list of scalars, or a list of
#         lambda functions for non-tanh scaling

# Approach 1: Use fluctuations history from TurbSim but scale by the LES profile
#from scipy import interpolate
#uu_LES_norm = inflowPlane.uu_profile / inflowPlane.uu_mean
#vv_LES_norm = inflowPlane.vv_profile / inflowPlane.vv_mean
#ww_LES_norm = inflowPlane.ww_profile / inflowPlane.ww_mean
#uu_fn = interpolate.interp1d(inflowPlane.z_profile, (1.0*uu_LES_norm)**0.5, kind='linear', fill_value='extrapolate')
#vv_fn = interpolate.interp1d(inflowPlane.z_profile, (1.0*vv_LES_norm)**0.5, kind='linear', fill_value='extrapolate')
#ww_fn = interpolate.interp1d(inflowPlane.z_profile, (1.0*ww_LES_norm)**0.5, kind='linear', fill_value='extrapolate')
#inflowPlane.setScaling(max_scaling=(uu_fn, vv_fn, ww_fn), output='variance_profile_scaled.dat')

# Approach 2: Assume ~constant up to specified height; blending with tanh function
inflowPlane.setScaling(
    tanh_z90=750.0,
    tanh_z50=600.0,
    max_scaling=(0.75*maxVarianceProfile/meanVariances)**0.5
)

# Write out boundary data--fluctuations plus the mean profile--for the OpenFOAM
# 'timeVaryingMappedFixedValue' boundary typically used in SOWFA wind plant
# simulations
# Note: Output is every (interval * inflowPlane.dt) seconds, starting at Tstart
inflowPlane.writeMappedBC(
    'boundaryData/west',
    interval=outputInterval,
    Tstart=0.0, Tend=1000.0,
    bcname='west',
    LESyfac=3, LESzfac=3 # if inflow plane spacing doesn't match the actual boundary
)

# Write out VTK set for visualizing the generated inflow (OPTIONAL)
#inflowPlane.writeVTKSeries(outputdir='vtkSeries',prefix='inflow',step=outputInterval)
inflowPlane.writeVTKSeries(outputdir='vtkSeries',prefix='scaled',step=outputInterval)

# Write out VTK 3D domain (dx=U*dt) for visualizing the generated inflow (OPTIONAL)
step = int(round(inflowPlane.dz/(inflowPlane.Umean*inflowPlane.dt))) # try to match dz spacing
inflowPlane.writeVTKSeriesAsBlock('refScaled.vtk',step=step)
```

