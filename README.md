detectorSimulations_v10
===================

The detectorSimulations_v10 package contains the Geant4 simulations for GRIFFIN, TIGRESS, and all of their auxiliary detectors.  Please note that in order to respect our non-disclosure agreements, all source containing third party IP have been omitted from this repo, and can be obtained from your colleagues directly at the lab.

# Release

version | DOI
--------|------
1.0.0   | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1123461.svg)](https://doi.org/10.5281/zenodo.1123461)

# Setup

### Requirements
detectorSimulations is confirmed to run on geant4.10.02 and ROOT 6.04.00. Do not use geant4.10.01 as it has a bug in the gamma de-excitation (see issue #35).

There is a new branch geant4.10.04, which was used to compile the simulation sucessfully against geant4.10.01.p01 (using ROOT 6.10/08). The main (untested) change to the branch is the removal of the explicit use of the obsolete G4FermiBreakup model in PhysListHadon. It has only been tested (so far) with a simulation of one million gamma rays of 1000 keV using only GRIFFIN (with 20 mm delrin). See also issue #40.

### Getting the code

To setup the simulation package on a computer with GEANT4 already present, just copy the code to your machine:

    git clone https://github.com/GRIFFINCollaboration/detectorSimulations_v10.git
    
Then you'll need to get the files containing our NDA-protected parameters. To do this register on gitlab.com, have your account added to the GRIFFINCollaboration, and register your ssh-keys with gitlab. Then you can run the script SetupSuppressed.sh (in the detectorSimulation_v10 folder). This script can either be run as is, which will install the suppressed files in a sub-folder "suppressed" and create symbolic links in the src directory, or you can give it the path of the directory where you want the suppressed files installed (this directory has to be empty!).

### Building

The build process is pretty standard for a geant simulation; in a build directory (ie any clean new directory that isn't the source directory), do 

```
cmake path/to/detectorSimulations
make clean
make
```

Keep in mind that cmake does not regenerate all the files it uses every time it runs!  So if something changes and this build process suddenly fails, try deleting the build directory and starting over.

### Setup FAQ

- Yes, you need both the secret suppressed files AND their unsuppressed equivalents, not just either / or.

# TI-STAR Branch

The following commands/details are specific to this branch and are related to the TI-STAR simulation.

## Usage 

### TI-STAR Geometry

These commands cover the physical TI-STAR and gas target geometry

TI-STAR detector

| Command | Brief Description |
| :------ | :---------------- |
| ``` /DetSys/det/addTISTARLayer ``` | Add a single TI-STAR layer | 
| ``` /DetSys/det/setTISTARSiDimensions x y z unit ``` | Set the silicon dimensions of the TI-STAR layer | 
| ``` /DetSys/det/setTISTARPCBDimensions x y z unit ``` | Set the PCB dimensions of the TI-STAR layer | 
| ``` /DetSys/det/setTISTARSiOffsetInPCB x y z unit ``` | Set the offset of the silicon layer within the PCB of the TI-STAR layer |
| ``` /DetSys/det/setTISTARPosition x y z unit ``` | Set the position of a single TI-STAR layer (built w/ ```/DetSys/det/addTISTARLayer```) | 
| ``` /DetSys/det/setTISTARPositionOffset x y z unit ``` | Set the position offset of the 4-strip TI-STAR layer, used to line up the inner 4-strip layer with the outer 2-strip layers  | 
| ``` /DetSys/det/setTISTARRotation x y z ``` | Set the rotation of a single TI-STAR layer - rotate by the x, y, then z axes - units in deg (built w/ ```/DetSys/det/addTISTARLayer```) | 
| ``` /DetSys/det/addTISTAR2StripLayer ``` | Add a 2-strip TI-STAR layer | 
| ``` /DetSys/det/addTISTAR4StripLayer ``` | Add a 4-strip TI-STAR layer |  
| ``` /DetSys/det/setTISTARDistFromBeam distance unit ``` | Set the distance from the beam axis - used by ```/DetSys/det/addTISTAR2StripLayer``` and ```/DetSys/det/addTISTAR4StripLayer``` | 
| ``` /DetSys/det/setTISTARGapZ distance unit``` | Set the distance between the two strips on either side of the beam axis (z) - used by ```/DetSys/det/addTISTAR4StripLayer``` | 
| ``` /DetSys/det/setTISTARSiCentered bool``` | Set boolean to build the TI-STAR layer silicon-centered (true) or entire chip [silicon + PCB] (false). Mainly used to center the outer 2-strip layers of the current TI-STAR design. |

Gas target

| Command | Brief Description |
| :------ | :---------------- |
| ``` /DetSys/app/addTISTARGasTarget ``` | build the TI-STAR gas target - parameters currently set via commands in ```/DetSys/det/miniball/``` directory | 

Vacuum chamber

| Command | Brief Description |
| :------ | :---------------- |
| ``` /DetSys/app/addTISTARVacuumChamber ``` | Add the TI-STAR vacuum chamber |
| ``` /DetSys/app/setTISTARVacuumChamberShape ``` | Set the TI-STAR vacuum chamber shape (box or cylinder) |
| ``` /DetSys/app/setTISTARVacuumChamberMaterial ``` | Set the TI-STAR vacuum chamber material |
| ``` /DetSys/app/setTISTARVacuumChamberBoxDimensions ``` | Set the TI-STAR vacuum chamber dimensions (if shape = "box") |
| ``` /DetSys/app/setTISTARVacuumChamberCylinderRadius ``` | Set the TI-STAR vacuum chamber radius (if shape = "cylinder") |
| ``` /DetSys/app/setTISTARVacuumChamberCylinderZ ``` | Set the TI-STAR vacuum chamber z-dimensions (if shape = "cylinder") |
| ``` /DetSys/app/setTISTARVacuumChamberExteriorMaterial ``` | Set the TI-STAR vacuum chamber exterior material |
| ``` /DetSys/app/setTISTARVacuumChamberExteriorThickness ``` | Set the TI-STAR vacuum chamber exterior layer thickness |


### Miniball/TRex Commands

These commands are all related to the code copied over from [TI-STAR/TRex](https://github.com/VinzenzBildstein/TI-STAR) and [Miniball](https://github.com/VinzenzBildstein/Miniball).

| Command | Brief Description |
| :------ | :---------------- |
| ``` /DetSys/miniball/SetPrimaryGenerator generator_name``` | Set which TRexBaseGenerator-derived generator to use. The choices are: ```TestSource```, ```Rutherford```, ```AngularDistribution```, ```AlphaSource```, or ```BeamIn```. |
| ``` /DetSys/miniball/SimulateEjectiles bool``` | Set bool for simulating ejectiles |
| ``` /DetSys/miniball/SimulateGammas bool``` | Set bool for simulating gamma rays |
| ``` /DetSys/miniball/IncludeEnergyResolution bool``` | Set bool for including an energy resolution - NOT CURRENTLY USED |
| ``` /DetSys/miniball/IncludeVacuumChamber bool``` | Set bool for including the vacuum chamber - NOT CURRENTLY USED |
| ``` /DetSys/miniball/SetVacuumChamberType name``` | Set the the vacuum chamber type - NOT CURRENTLY USED |
| ``` /DetSys/miniball/SetVacuumChamberGas name``` | Set name of the the vacuum chamber gas - NOT CURRENTLY USED |
| ``` /DetSys/miniball/SetTestSourceEnergy energy unit``` | Set the energy for the test source - used by  ```TRexTestSource``` |
| ``` /DetSys/miniball/SetBeamEnergy energy unit``` | Set beam energy - used by ```TRexBeam``` derived sources |
| ``` /DetSys/miniball/SetBeamWidth length unit``` | Set beam width (unless reactionZ vs beam radius has been set via /DetSys/miniball/SetReactionZvsRadiusFile) - used by ```TRexBeam``` derived sources |
| ``` /DetSys/miniball/SetThetaCmMin theta``` | Set the minimum theta value, value in radians - used by ```TRexRutherford``` and ```TRexAngularDistribution``` |
| ``` /DetSys/miniball/SetProjectileName name```| Set the name of the Projectile (e.g. 132Sn) |
| ``` /DetSys/miniball/SetProjectileZ Z``` | Set Z for the Projectile |
| ``` /DetSys/miniball/SetProjectileA A``` | Set A for the Projectile |
| ``` /DetSys/miniball/SetTargetName name```| Set the name of the Target (e.g. 2H) |
| ``` /DetSys/miniball/SetTargetZ Z``` | Set Z for the Target |
| ``` /DetSys/miniball/SetTargetA A``` | Set A for the Target |
| ``` /DetSys/miniball/SetEjectileName name```| Set the name of the Ejectile (e.g. 133Sn) |
| ``` /DetSys/miniball/SetEjectileZ Z``` | Set Z for the Ejectile |
| ``` /DetSys/miniball/SetEjectileA A``` | Set A for the Ejectile |
| ``` /DetSys/miniball/SetRecoilName name```| Set the name of the Recoil (e.g. 1H) |
| ``` /DetSys/miniball/SetRecoilZ Z``` | Set Z for the Recoil |
| ``` /DetSys/miniball/SetRecoilA A``` | Set A for the Recoil |
| ``` /DetSys/miniball/SetTargetMaterialName name``` | Set the target material name (e.g. 2H) |
| ``` /DetSys/miniball/SetTargetAtomicRatio ratio``` | Set the target atomic ratio - NOT CURRENTLY USED |
| ``` /DetSys/miniball/SetTransferOrCoulexProbability prob``` | Set the transfer/coulex probability |
| ``` /DetSys/miniball/SetLevelFile filename``` | Set the name/location of the level file |
| ``` /DetSys/miniball/SetAngularDistributionFile filename``` | Set the name/location of the angular distribution file |
| ``` /DetSys/miniball/SetMassFile filename``` | Set the name/location of the mass file |
| ``` /DetSys/miniball/SetCrossSectionFile filename``` | Set the name/location of the cross-section file |
| ``` /DetSys/miniball/SetReactionZvsRadiusFile filename``` | Set the name/location of the reactionZ vs. beam radius file |
| ``` /DetSys/miniball/SetReactionZDistributionFile filename``` | Set the name/location of the reactionZ probability distribution file |
| ``` /DetSys/miniball/SetAlphaSourceDiameter diameter unit``` | Set the diameter of the alpha source - used by ```TRexAlphaSource``` |
| ``` /DetSys/miniball/SetAlphaSourceThickness thickness unit``` | Set the thickness of the alpha source - used by ```TRexAlphaSource``` |
| ``` /DetSys/miniball/SetTargetDiameter diameter unit``` | Set the diameter of the target |
| ``` /DetSys/miniball/SetTargetThickness thickness ``` | Set the thickness of the target (in mg/cm2) |
| ``` /DetSys/miniball/SetGasTargetLength length unit``` | Set the length of the gas target |
| ``` /DetSys/miniball/SetTargetPressure pressure unit``` | Set the gas target pressure |
| ``` /DetSys/miniball/SetTargetMaterialDensity density``` | Set the gas target material density (in mg/cm3) |
| ``` /DetSys/miniball/SetTargetMylarThickness thickness unit``` | Set the thickness of the target mylar foil |
| ``` /DetSys/miniball/SetTargetBeWindowThickness thickness unit``` | Set the thickness of the target Be window |
| ``` /DetSys/miniball/Print``` | Print the values stored in ```TRexSettings``` |

### Physics Lists

For TI-STAR, we typically use the QGSP_BIC physics list, but we have also added the physics list PhysListEmStandardNR from the example TestEm7 which includes G4ScreenedNuclearRecoil.

| Command | Brief Description |
| :------ | :---------------- |
| ``` /DetSys/phys/SelectPhysics physics_list ``` | Select a physics list - for TI-STAR we recommend ```QGSP_BIC``` by default, or ```standardNR``` for the G4ScreenedNuclearRecoil process  | 

# Usage

The following might be out of date with the change to geant4.10, we haven't had time yet to verify the information below.

### Particle Emission

#### General
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/gun/energy double unit ``` | Set energy of particle | 1000 keV |
| ``` /DetSys/gun/particle string ``` | Set particle type (e-, e+, gamma, proton, alpha) | gamma |
| ``` /DetSys/gun/ion Z A E ``` | Set ion type (excitation energy E in keV) | |
| ``` /DetSys/gun/direction x y z ``` | Set momentum direction | |
| ``` /DetSys/gun/position x y z  unit``` | Set particle position | 0.0 0.0 0.0 mm |
| ``` /DetSys/gun/radius r unit``` | Set source radius | 0.0 mm |
| ``` /DetSys/gun/energyrange min max step ``` | Set energy (keV) of particle, loops from min to max | |

#### Decay Schemes
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/gun/betaPlusEmission filename ``` | Simulate beta plus decay with energy distribution input file | |
| ``` /DetSys/gun/betaMinusEmission filename ``` | Simulate beta negative decay with energy distribution input file | |
| ``` /DetSys/gun/polarization double``` | Set Polarization of Nuclei (before radioactiveBetaDecay) | |
| ``` /DetSys/gun/radioactiveBetaDecay directory ``` | Simulate complete beta negative decay with simulation directory | |
| ``` /DetSys/gun/emitBetaParticle 0/1 ``` | Emit Beta Particle? True/False | |
| ``` /DetSys/gun/includeXRayInputFileKShell 0/1 ``` | Emit X-rays from K-shell vacancies using input file? True/False | |
| ``` /DetSys/gun/includeXRayInputFileLShell 0/1 ``` | Emit X-rays from L-shell vacancies using input file? True/False | |
| ``` /DetSys/gun/includeXRayInputFileMShell 0/1 ``` | Emit X-rays from M-shell vacancies using input file? True/False | |
| ``` /DetSys/gun/radioactiveDecayHalflife double ``` | Half-life of radioactive isotope simulation (seconds) | |
| ``` /DetSys/gun/numberOfRadioactiveNuclei int ``` | Set the number of radioactive nuclei | |
| ``` /DetSys/gun/radioactiveSourceDecay filename ``` | Simulate source decay with a file containing the decay data | |

#### Kinematics
Whenever an electron or positron is emitted, you can choose to simulate kinematic effects (i.e. shifting of energy due to relativistic speeds of ions).  You must specify the ion with ``` /DetSys/gun/ion Z A E ``` first.

| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/gun/simulateKinematics 0/1 ``` | Choose to simulate kinematic effects | False |
| ``` /DetSys/gun/kinematicsBetaValue double ```* | Specify beta value of incident ion | 0 |
| ``` /DetSys/gun/kinematicsIonEnergy double unit ```* | Specify energy of heavy ion | 0 MeV |

*Only one of these should be specified

### Detector Specific

#### GRIFFIN
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/addGriffinForward int ``` | Add Detection System GriffinForward |  |
| ``` /DetSys/det/addGriffinForwardDetector int ``` | Add GriffinForward Detector |  |
| ``` /DetSys/det/addGriffinBack int ``` | Add Detection System GriffinBack |  |
| ``` /DetSys/det/addGriffinBackDetector int ``` | Add GriffinBack Detector |  |
| ``` /DetSys/det/UseTIGRESSPositions ``` | Use TIGRESS detector positions rather than GRIFFIN | False |
| ------ | ---------------- | ------ |
| ``` /DetSys/det/addGriffinCustomDetector 0 ``` | Adds a detector using the paramaters specified |  |
| ``` /DetSys/det/SetCustomShieldsPresent 0/1 ``` | Selects whether or not the detector suppressors are included | True |
| ``` /DetSys/det/SetCustomRadialDistance double unit ``` | Selects the radial distance for the detector from the origin |  |
| ``` /DetSys/det/SetCustomExtensionSuppressorLocation 0/1 ``` | Selects a position for the extension suppressors. Either forward (0) or back (1) | Forward (0) |
| ``` /DetSys/det/SetCustomPosition det_num pos_num null ``` | Sets the position for the detector placed in the next call to addGriffinCustom |  |
| ``` /DetSys/det/SetCustomDeadLayer det_num cry_num dead_layer ``` | Sets the dead layer for the crystal specified (used in any following calls to addGriffinCustom) |  |
| ``` /DetSys/det/includeGriffinHevimet 0/1 ``` | Includes the Hevimet for a Griffin detector | False |
| ``` /DetSys/det/addGriffinCustom int ``` | Adds a detection system using the paramaters specified |  |

#### SPICE
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/app/addSpiceTargetChamber string ``` | Add SPICE target chamber* |  |
| ``` /DetSys/Spice/setResolution double double ``` |Set resolution of SPICE Si(Li)  |  |
| ``` /DetSys/det/addSpice int``` | Add Si(Li) detector |  |
| ``` /DetSys/det/addS3 ``` | Add SPICE S3 detector |  |

#### 8PI
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/add8pi int ``` | Add Detection System 8pi |  |
| ``` /DetSys/det/add8piDetector int ``` | Add 8pi Detector |  |
| ``` /DetSys/app/add8piVacuumChamber ``` | Add 8pi vacuum chamber |  |
| ``` /DetSys/app/add8piVacuumChamberAuxMatShell double unit ``` | Add AuxMat shell around 8pi vacuum chamber with specified thickness |  |

#### Other
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/addGammaTracking int ``` | Add Detection System GammaTracking |  |
| ``` /DetSys/det/addSodiumIodide int ``` | Add Detection System SodiumIodide |  |
| ``` /DetSys/det/addLanthanumBromide int ``` | Add Detection System LanthanumBromide |  |
| ``` /DetSys/det/addSceptar int ``` | Add Detection System Sceptar |  |
| ``` /DetSys/det/addPaces int ``` | Add Detection System Paces |  |

### Detector General

``` /DetSys/det/update ``` must be called before using ``` beamOn ``` if the geometrical values
have been altered.

#### World Volume (i.e. experimental hall)
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/world/material string ``` | Select material for the world | G4_AIR |
| ``` /DetSys/world/dimensions x y z unit ``` | Set world dimensions | 10m x 10m x 10m |
| ``` /DetSys/world/vis 0/1``` | Set world visibility (depreciated) | False |

#### Generic Target
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/app/genericTarget string ``` | Create a target with specified material |  |
| ``` /DetSys/app/genericTargetDimensions x y z unit ``` | Set target dimensions |  |
| ``` /DetSys/app/genericTargetPosition x y z unit ``` | Set target position |  |
Requires all three commands to build

#### Field Box
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/app/fieldBoxMaterial string ``` | Create a field box with specified material |  |
| ``` /DetSys/app/fieldBoxDimensions x y z unit ``` | Set field box dimensions |  |
| ``` /DetSys/app/fieldBoxPosition x y z unit ``` | Set field box position |  |
| ``` /DetSys/app/fieldBoxMagneticField x y z unit ``` | Set field box magnetic field |  |
Requires all four commands to build

#### Box
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/boxMat string ``` | Set box material | G4_WATER |
| ``` /DetSys/det/boxThickness double unit ``` | Set box thickness | 0.0 mm |
| ``` /DetSys/det/boxInnerDimensions x y z unit ``` | Set box inner dimensions | 0.0 0.0 0.0 mm |
| ``` /DetSys/det/boxColour r g b``` | Set box colour | 0.0 0.0 1.0 |
| ``` /DetSys/det/addBox ``` | Build/add box (if thickness is not 0) |  |

#### Grid
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/det/gridMat string ``` | Set grid material | G4_WATER |
| ``` /DetSys/det/gridSize double unit ``` | Set grid size | 0.0 mm |
| ``` /DetSys/det/gridDimensions x y z unit ``` | Set grid dimensions | 0.0 0.0 0.0 mm |
| ``` /DetSys/det/gridColour r g b ``` | Set grid colour | 1.0 0.0 0.0 |
| ``` /DetSys/det/addGrid ``` | Build/add grid (if grid size is not 0) |  |

#### Magnetic Fields
| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ``` /DetSys/world/magneticField x y z unit``` | Set world magnetic field (depreciated) | 0, 0, 0 |
| ``` /DetSys/world/tabMagneticField filename ``` | Set tabulated magnetic field* | Disabled |
| ``` /DetSys/gun/coneRadius radius unit ``` | Set cone radius |  |


*Used for SPICE: these files are 50Mb each and are different for each lens so until a better solution comes along they will be kept on the network at TRIUMF (email Mohamad, Lee, or Michael for precise location)


| Command | Brief Description | Default |
| :------ | :---------------- | :------ |
| ```  ``` |  |  |
| ```  ``` |  |  |
| ```  ``` |  |  |
| ```  ``` |  |  |

