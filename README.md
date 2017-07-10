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
from inflow.synthetic.turbsim import bts

# Read plane of fluctuating velocities to be superimposed
inflowPlane = bts('class3C.bts')

# Calculate mean variance (and optionally output standard deviations)
inflowPlane.calculateRMS(output='rms.txt')

# Resize the domain (OPTIONAL)
# Note: This should be done prior to calling readMeanProfile so
#       that the mean profile is set on the new resized domain.
inflowPlane.extendZ(zMin=0.0, zMax=1000.0)
inflowPlane.tileY(Npanels=3, mirror=True)

# Set up mean profiles
# Note: This sets the *_profile properties of inflowPlane 
# Note: Can also call readProfile('averagingProfiles.csv',delim=',')
#       or setMeanProfile(U,V,T) and setVarianceProfile(uu,vv,ww)
inflowPlane.readMeanProfile('U.dat','V.dat','T.dat')
inflowPlane.readVarianceProfile('uu.dat','vv.dat','ww.dat')

# Set up the TKE profile
inflowPlane.setTkeProfile(kprofile=lambda z: 0.1) # constant field

# Since the modeled velocity fluctuations have no height dependence,
# scale the velocities so that fluctuations decay above some height, zi
# (where zi == inversion height >= the ABL height)
# Note: 'max_scaling' can be a constant, a list of scalars, or a list of
        lambda functions for non-tanh scaling
maxVarianceProfile = np.array([ np.max(inflowPlane.uu_profile),
                                np.max(inflowPlane.vv_profile),
                                np.max(inflowPlane.ww_profile) ]) # input
meanVariances = np.array([ inflowPlane.uu_mean,
                           inflowPlane.vv_mean,
                           inflowPlane.ww_mean ]) # output from turbsim
inflowPlane.setScaling(
    tanh_z90=750.0,
    tanh_z50=600.0,
    max_scaling=(0.75*maxVarianceProfile/meanVariances)**0.5
)

# Write out boundary data for the OpenFOAM 'timeVaryingMappedFixedValue'
# boundary typically used in SOWFA wind plant simulations
# Note: Output is every (interval * dt) seconds 
inflowPlane.writeMappedBC(
    'boundaryData/west',
    interval=20,
    Tstart=0.0,
    bcname='west'
)

# Write out VTK set for visualizing the Turbsim solution (OPTIONAL)
inflowPlane.writeVTKSeries(Umean=8.0, step=1)
```

