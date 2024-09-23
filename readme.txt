Opticalizer code can be ran on its own or used in a job array or my perosnal favorite a meta farm.
This code takes non optical CAl files and makes a ntp 

Opticalizer.py is a python code used to transform non optical CaL file into NTP with many details behind what has caused each event. 
It is currently used mainly for high energy events E>1MeV. First the code takes each event and simulates a 13.6us DAQ trigger window.
If a particle falls within this window its origin is anlyzed. it can belong to one of 4 groups, an Inelastic event, Neutron capture, 
Source gamma, or other. (you may add your own interaction of interst seemlessly) from there it records all peritnent information about the event.
if there is more energy deposits past the 13.6 us widnow it again repeats the process until ever enegry deposit is acounted for.
the code then repeats the process for the next Cal event.

**note for nuclear recoils a quenching factor of 0.3 is used

Nomenclature
nc=neutron capture
ni=neutron inelastic
cg=capture gamma

Neutron captures:
capture element (Z,A and Q) (nc_Z,nc_A,nc_Q) , the neutron energy before capture (nc_ne),the total energy released in gammas (nc_ge), 
the total energy released (nc_e). the top 3 highest energy capture gammas energy and momentum (solid angle theta & phi) (nc_cg1,nc_cg1_theta,nc_cg1_phi). 
boolean if a Ar capture released a 4.7 MeV gamma (nc_ar47==1) else (nc_ar47==0) and similarly for Fe 7.6 MeV (nc_fe76== 0 or 1)

Neutron inelastics:
Scatter element (ni_Z,ni_A), neutorn energy before scatter (ni_ne), energy released (ni_e), energy released in gammas (ni_ge),

Important variable:
Energy :standard energy deposited in the trigger window.
EnergyRS takes a random number from a gaussian distribution defined by the response function in Carl Rethemeier's thesis. (RS stands for resolution)
qPE takes carls energy response function and transforms energy to qPE.
Energy_ni, Energy_nc, Energy_sg : this is the energy that was deposited by a neutron inelastic, neutron capture or source gamma event

fprompt: (energy by nuclear recoils)/(total energy )

DTw: delta Time of energy deposit between a neutron capture and a neutron inelastic or source gamma. time of one type is averaged by energy weight(see "weighted_average" function )
primary_ne : neutron energy of the start particles
primary_sg : gamma energy of the start particles

fileID and entryID recoords the cal file and he entry that the vent takes palce from os you can easily go back and find a particualr event of interest in the cal files.




