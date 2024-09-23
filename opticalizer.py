# By Trevor H

# to run  
# from command line insert the name of the cal file to read and the ntp file to be created  
# eg: opticalizer.py ntpfile.root calfile.root

#typically takes 30 min to do a 80k event cal file

#the rat5 environment supports both rootpy and numpy 
# cannot run on python 3 as rootpy is deprecated  


#to set up on cedar 
#PATH=$PATH:$HOME/.local/bin:$HOME/bin 
#export PATH 
#DEAPHOME=/project/6004969 
#PATH=$PATH:$DEAPHOME/software/bin: 
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DEAPHOME/software/lib 
#MANPATH=$MANPATH:$DEAPHOME/software/share/man 
#export DEAPHOME PATH LD_LIBRARY_PATH MANPATH MACHTYPE OSTYPE 
#PATH=$PATH:$HOME/.local/bin:$HOME/bin 
#source $DEAPHOME/software/ratcage/env.sh 
#export PATH 
 


import ROOT
import rat
from os import listdir
from os.path import isfile, join
import math
import re
import math
from rootpy.tree import Tree, TreeModel, FloatCol, IntCol
from rootpy.io import root_open
from random import gauss
import sys
import time
import numpy as np
from math import sqrt
start_time = time.time()

# from command line insert the name of the cal file to read and the ntp file to be created  eg: python GeantFixer.py ntpfile.root calfile.root
ntp_file = sys.argv[1]
cal_file=sys.argv[2]


print("reading ",cal_file)
print("printing to ",ntp_file)

#get file number:
l= re.findall(r'[0-900]+',cal_file)
if len(l)>=1:
    filenumber=int(l[-1])
else:
    filenumber=0

ntp = root_open(ntp_file, "recreate")

class Event(TreeModel):
        #energy deposited from Neuton inelastics,neutron cpatures and Source gammas
        qPE=FloatCol()
        energy=FloatCol()
        energyRS=FloatCol()
        primary_ge=FloatCol()
        primary_ne=FloatCol()
        
        energy_ni=FloatCol()
        energy_nc=FloatCol()
        energy_sg=FloatCol()
        fprompt=FloatCol()
        multi_trig=FloatCol()
        DTw=FloatCol()
        

        fileID=IntCol()
        entryID=IntCol()
        
        #neutron captures
        nc_e=FloatCol()
        nc_ge=FloatCol()
        nc_ne=FloatCol()
        nc_n_theta=FloatCol()
        nc_n_phi=FloatCol()
        
        nc_Q=FloatCol()
        nc_A=IntCol()
        nc_Z=IntCol()
        
        nc_ar47=FloatCol()
        nc_fe76=FloatCol()
        nc_multi=IntCol()
        
        nc_cg1=FloatCol()
        nc_cg2=FloatCol()
        nc_cg3=FloatCol()
        nc_cg1_theta=FloatCol()
        nc_cg1_phi=FloatCol()
        nc_cg2_theta=FloatCol()
        nc_cg2_phi=FloatCol()
        nc_cg3_theta=FloatCol()
        nc_cg3_phi=FloatCol()
        
        
        #neutron inelastics 
        ni_e=FloatCol()
        ni_ge=FloatCol()
        ni_ne=FloatCol()
        ni_A=IntCol()
        ni_Z=IntCol()
        ni_multi=IntCol()
        
        
tree = Tree("data", model=Event)

# creates a list of possible elements for neutron capture
#(A,Z,Q) Q is in MEV/C^2
H1=(1,1,-2.22456)
H2=(2,1,-6.257)
ar40=(40,18,-6.0989)
ar36=(36,18,-8.787)
ar38=(38,18,-6.599)
O16=(16,8,-4.143)
Cu63=(63,29,-7.915852)
Cu65=(65,29,-7.066)
C12=(12,6,-4.946)
C13=(13,6,-8.176)
I127=(127,53,-6.826)
B10=(10,5,-11.454221)
B11=(11,5,-3.36966)
Cl35=(35,17,-8.58)
Cl37=(37,17,-6.108)
#chromium
cr50=(50,24,-9.2606)
cr52=(52,24,-7.9391)
cr53=(53,24,-9.7191)
cr54=(54,24,-6.2463)
#Iron
Fe54=(54,26,-9.2982)
Fe56=(56,26,-7.6461)
Fe57=(57,26,-10.0446)
#manganes
Mn55=(55,25,-7.2704)
#nickel
Ni58=(58,28,-8.9993)
Ni60=(60,28,-7.8201)
Ni61=(61,28,-10.5965)
Ni62=(62,28,-6.8378)
Ni64=(64,28,-6.0981)

Qvaluelist=[H1,H2,ar36,ar38,ar40,O16,Cu63,Cu65,C12,C13,I127,B10,B11,Cl35,Cl37,cr50,cr52,cr53,cr54,Fe54,Fe56,Fe57,Mn55,Ni58,Ni60,Ni61,Ni62,Ni64]

def findq(Z,A):
    #print(Z,A)
    for element in Qvaluelist:
        if element[0]==int(A) and element[1]==int(Z):
            return element[2]
    return 0

#open cal file may have to change the "Get" argument if the tree isnt called "T"
file = ROOT.TFile(cal_file)
T = file.Get("T")#_satCorr")
print("numb events ", T.GetEntries())
for entry in range(T.GetEntries()):  
    T.GetEntry(entry)
   
    if T.ds.GetMC().GetMCSummary().GetTotalScintEdep()==0 :
        continue
    mc = T.ds.GetMC()
    
    
    index_id = dict()
    parent_id = dict()

    for w in range(mc.GetMCTrackCount()):
        track = mc.GetMCTrack(w)
        index_id[track.GetTrackID()] = w
        parent_id[track.GetTrackID()] = track.GetParentID()
        
    def track_info(track):
            parent=track.GetParentID()
            track_ID = track.GetTrackID()
            pdg = track.GetPDGCode()
            process = track.GetCreatorProcessName()
            energy=track.GetMCTrackStep(0).GetKE()
            location=track.GetMCTrackStep(0).GetEndpoint().X(),track.GetMCTrackStep(0).GetEndpoint().Y(),track.GetMCTrackStep(0).GetEndpoint().Z()
            time=track.GetMCTrackStep(0).GetGlobalTime()
            momentum=step_info(track.GetMCTrackStep(0))[3]
            return (parent,track_ID,pdg,process,energy,location,momentum,time)

    def step_info(track_step):
        energy=track_step.GetKE()
        location=track_step.GetEndpoint().X(),track_step.GetEndpoint().Y(),track_step.GetEndpoint().Z()
        time=track_step.GetGlobalTime()
        volume=track_step.GetVolume()
        momentum=track_step.GetMomentum()
        return (energy,location,volume,(momentum.Mag(),momentum.Theta(),momentum.Phi()),time)


    


    def DAQ_triggers(alist): #returns a list of list with each list containing a trigger deposit
            events_sorted= sorted(alist, key=lambda tup: tup[0])  #sort the list acording to time
            trig_list=[]
            final_trigs=[]
            
            if len(events_sorted)>0:
                start_point=events_sorted[0][0]
                for event in events_sorted:
                    if event[0]<=start_point+13e6:
                        trig_list.append(event)
                    else:
                        start_point=event[0]
                        final_trigs.append(trig_list)
                        trig_list=[]
                
                if len(trig_list)>0:
                    final_trigs.append(trig_list)
            
            return final_trigs
            

    def findevents():              
        E_and_T=[]
        em=[11,-11,-13,13,22] #pdg codes of ER type particles electron muon, gamma
        for i in range(mc.GetMCTrackCount()):
            track=mc.GetMCTrack(i)
            t_info=track_info(track)
            
            for s in range(track.GetMCTrackStepCount()):
                step = track.GetMCTrackStep(s)
                if (step.GetVolume()=="cryoliquid" and s !=0 ):
                    
                    if t_info[2] not in em and step.GetTotalEdep()>0:
                        E_and_T.append((step.GetGlobalTime(),step.GetTotalEdep()*0.3))
                            
                    elif step.GetTotalEdep()>0:
                        E_and_T.append((step.GetGlobalTime(),step.GetTotalEdep()))
        
        trigs=DAQ_triggers(E_and_T)
        return trigs
                    

    def inelastic_event(ti):
        parent_id=ti[0]
        energy=ti[4]
        process=ti[3]
        times=ti[-1]
        
        #find capture element
        x=re.findall(r'\d+',process)
        Z_NI=-1
        A_NI=-1
        Q=-1
        if len(x)>1:
            Z_NI=int(x[0])
            A_NI=int(x[1])
        
        #neutron before inelastic
        neutron_track=mc.GetMCTrack(index_id[parent_id])
        s=neutron_track.GetMCTrackStepCount()
        step = neutron_track.GetMCTrackStep(s-2)
        
        nstep=step_info(step)
        neutron_energy=nstep[0]
        #neutron_momentum=nstep[3]
        
        #find other daughters
        total_energy=0
        gamma_e=0
        
        for i in range(mc.GetMCTrackCount()):
            track=mc.GetMCTrack(i)
            t_info=track_info(track)
            if t_info[0]==parent_id and process==t_info[3]:
                energy=t_info[4]
                total_energy+=energy
                if t_info[2]==22:
                    gamma_e+=energy
        
        info=(neutron_energy,total_energy,gamma_e,Z_NI,A_NI)
        return info
                                

    def ncap_event(ti):
        
        Ar47=0
        Fe76=0
        parent_id=ti[0]
        energy=ti[4]
        process=ti[3]
        
        #find capture element
        x=re.findall(r'\d+',process)
        Z_NC=-1
        A_NC=-1
        Q=-1
        if len(x)>1:
            Z_NC=int(x[0])
            A_NC=int(x[1])
            Q=findq(Z_NC,A_NC)
        #check if its the characteristic gamma
        if Z_NC==26 and A_NC==56 and 7.6<energy<7.69: Fe76=1
        if Z_NC==18 and A_NC==40 and 4.7<energy<7.79: Ar47=1
        
        #find neutron before capture
        neutron_track=mc.GetMCTrack(index_id[parent_id])
        s=neutron_track.GetMCTrackStepCount()
        step = neutron_track.GetMCTrackStep(s-2)
        
        nstep=step_info(step)
        neutron_energy=nstep[0]
        neutron_momentum=nstep[3]
        
        
        
        #find other neutron capture related information:
        charcteristic_gammas=[]
        total_energy=0
        gamma_e=0
        for i in range(mc.GetMCTrackCount()):
            track=mc.GetMCTrack(i)
            t_info=track_info(track)
            if t_info[0]==parent_id and process==t_info[3]:
                energy=t_info[4]
                total_energy+=energy
                if t_info[2]==22:
                    gamma_e+=energy
                    momentum=t_info[6]
                    charcteristic_gammas.append((energy,momentum))
        
        charcteristic_gammas.sort(key=lambda tup: tup[0])
        cg1=0
        cg2=0
        cg3=0
        cg1_phi=0
        cg1_theta=0
        cg2_phi=0
        cg2_theta=0
        cg3_theta=0
        cg3_phi=0
        cg1=charcteristic_gammas[-1][0]
        cg1_theta=charcteristic_gammas[-1][1][1]
        cg1_phi=charcteristic_gammas[-1][1][2]
        if len(charcteristic_gammas)>1:
            cg2=charcteristic_gammas[-2][0]
            cg2_theta=charcteristic_gammas[-2][1][1]
            cg2_phi=charcteristic_gammas[-2][1][2]
        if len(charcteristic_gammas)>2:
            cg3=charcteristic_gammas[-3][0]
            cg3_theta=charcteristic_gammas[-3][1][1]
            cg3_phi=charcteristic_gammas[-3][1][2]
            
        
        #done everything:
        info=(neutron_energy,neutron_momentum,Z_NC,A_NC,Q,total_energy,gamma_e,cg1,cg2,cg3,cg1_theta,cg1_phi,cg2_theta,cg2_phi,cg3_theta,cg3_phi) 
        
        return (info,Ar47,Fe76)
        

    def event_finder(track):
        ti=track_info(track)
        pid=ti[0]
        id = ti[1]
        pdg=ti[2]
        process=ti[3]
        
        found = False
        origin="not found" 
        while(not found):
            
            if parent_id[id] == 0:
                origin="Other"
                found=True
                info=0
                if pdg==22:
                    origin='Source_gamma'
                    
                
            elif pdg==22 and 'ncapture' in process.lower():
                origin="ncapture"
                info=ncap_event(ti)
                found=True
            elif pdg!=2112 and "inelastic" in process.lower():
                origin="neutron_inelastics"
                info=inelastic_event(ti)
                found=True
            else:
                track=track=mc.GetMCTrack(index_id[pid])
                ti=track_info(track)
                pid=ti[0]
                id = ti[1]
                pdg=ti[2]
                process=ti[3]
        return origin,info
        #returns source ,nc,ni, other, and then the description of what happened


    def Mev_to_qPE(E):
        C_E_Scale= 548.0 
        L_E_Scale= 6853.7 
        Q_E_Scale= -19.51 
        #L_Res= 2.53524 
        #Q_Res=4.4882e-05 
        # sigma^2 = [Linear Res.] * mu + [Quadratic Res.] * mu^2
        qPE = C_E_Scale + L_E_Scale*E + Q_E_Scale*(E**2)
        return qPE      

    def gaussres(qPE):
        #define detector respnse from carls STR 
        lin_res=2.972
        quad_res=0.00003
        sigma_sq=quad_res*qPE*qPE+(1+lin_res)*qPE
        sigma=sqrt(sigma_sq)
        mu=qPE
        s = np.random.normal(mu, sigma)
        mev=(-6853.7+sqrt(6853.7*6853.7+4*19.51*(548-s)))/(2*-19.51)
        
        return mev       

    def add_energy(e_list):
        energy=0
        for item in e_list:
            energy+=item[1]
        return energy

    def weighted_average(list1):
        num=0
        denom=0
        for item in list1:
            num+=item[0]*item[1]
            denom+=item[1]
        
        if denom!=0:
            return num/denom
        else: return "error"
            
            
            
    def deltaT(list1,list2): 
        w_1=weighted_average(list1)
        w_2=weighted_average(list2)
        if w_1=='error' or w_2=='error':
            wdelta=-1
        else:
            wdelta=w_2-w_1
        
        return wdelta

    def start_conditions():
        ge=0
        ne=0
        for i in range(mc.GetMCTrackCount()):
            track=mc.GetMCTrack(i)
            ti=track_info(track)
            if ti[0]==0:
                if ti[2]==22:
                    ge+=ti[4]
                if ti[2]==2112:
                    ne=ti[4]
        return ne,ge



    def fill_tree(Energy,Fprompt,neutron_start,gamma_start,Neutron_inelastics_E_T,ncapture_E_T,Source_E_T,ncap_info,inelastic_info,argon47,iron76):
        
        tree.qPE=Mev_to_qPE(Energy)
        tree.energyRS=gaussres(Mev_to_qPE(Energy))
        tree.energy=Energy
        tree.fprompt=Fprompt
        tree.primary_ne=neutron_start
        tree.primary_ge=gamma_start
        
            
        tree.energy_ni=add_energy(Neutron_inelastics_E_T)
        tree.energy_nc=add_energy(ncapture_E_T)
        tree.energy_sg=add_energy(Source_E_T)

        if add_energy(ncapture_E_T)>0 and add_energy(Neutron_inelastics_E_T)>0: tree.DTw = deltaT(Neutron_inelastics_E_T,ncapture_E_T)
        elif add_energy(ncapture_E_T)>0 and add_energy(Source_E_T)>0: tree.DTw = deltaT(Source_E_T,ncapture_E_T)
        else:tree.DTw=-1
        
        tree.fileID=filenumber
        tree.entryID=entry
        
        #ncapture
        if len(ncap_info)>1:tree.multievent=len(ncap_info)  
        else: tree.nc_multi=0
        
        if len(ncap_info)==1: 
            ncap_info=ncap_info[0]
            tree.nc_ne=ncap_info[0]
            tree.nc_n_theta=ncap_info[1][1]
            tree.nc_n_phi=ncap_info[1][2]
            
            tree.nc_Z=ncap_info[2]
            tree.nc_A=ncap_info[3]
            tree.nc_Q=ncap_info[4]
            
            tree.nc_e=ncap_info[5]
            tree.nc_ge=ncap_info[6]
            
            tree.nc_ar47=0
            tree.nc_fe76=0
            if argon47:tree.nc_ar47=1
            if iron76:tree.nc_fe76=1
            
            tree.nc_cg1=ncap_info[7]
            tree.nc_cg2=ncap_info[8]
            tree.nc_cg3=ncap_info[9]
            
            tree.nc_cg1_theta=ncap_info[10]
            tree.nc_cg1_phi=ncap_info[11]

            tree.nc_cg2_theta=ncap_info[12]
            tree.nc_cg2_phi=ncap_info[13]

            tree.nc_cg3_theta=ncap_info[14]
            tree.nc_cg3_phi=ncap_info[15]
                    
        else:
            
            tree.nc_ne=0
            tree.nc_n_theta=0
            tree.nc_n_phi=0
            
            tree.nc_Q=0
            tree.nc_A=0
            tree.nc_Z=0
            
            tree.nc_e=0
            tree.nc_ge=0
            
            tree.nc_cg1=0
            tree.nc_cg2=0
            tree.nc_cg3=0
            
            tree.nc_cg1_theta=0
            tree.nc_cg2_theta=0
            tree.nc_cg3_theta=0
            
            tree.nc_cg1_phi=0
            tree.nc_cg2_phi=0
            tree.nc_cg3_phi=0  
        
        #neutron inelastic 
        if len(inelastic_info)>1:tree.ni_multi=len(inelastic_info)
        else:tree.ni_multi=0
                
        if len(inelastic_info)==1:
            inelastic_info=inelastic_info[0]
            tree.ni_ne=inelastic_info[0]
            tree.ni_e=inelastic_info[1]
            tree.ni_ge=inelastic_info[2]
            tree.ni_Z=inelastic_info[3]
            tree.ni_A=inelastic_info[4]
            
                
        else:
            tree.ni_e=0
            tree.ni_ge=0
            tree.ni_ne=0
            tree.ni_A=0
            tree.ni_Z=0
            
        tree.fill()
                        
                
    def event_description(trigger):
        #print(trigger)
        Energy=add_energy(trigger)
        initial_neutron,Source_gamma=start_conditions()
        
        Neutron_inelastics_E_T=[]
        ncapture_E_T=[]
        Source_E_T=[]
        
        fprompt_e=[]
        ncap_info=[]
        #Ar47 Fe76
        argon47=False
        iron76=False
        inelastic_info=[]
        
        em=[11,-11,-13,13,22] #pdg codes of ER type particles electron muon, gamma
        for i in range(mc.GetMCTrackCount()):
            track=mc.GetMCTrack(i)
            t_info=track_info(track)
            found=False
            source=0
            info=0
            for s in range(track.GetMCTrackStepCount()):
                step = track.GetMCTrackStep(s)
                if (step.GetVolume()=="cryoliquid" and s !=0 ):
                    
                    if t_info[2] not in em and step.GetTotalEdep()>0 and (step.GetGlobalTime(),step.GetTotalEdep()*0.3) in trigger:
                        if not found:
                            source,info=event_finder(track)
                        fprompt_e.append((step.GetGlobalTime(),step.GetTotalEdep()*0.3))
                        
                        
                        if source=="Source_gamma":
                            Source_E_T.append((step.GetGlobalTime(),step.GetTotalEdep()))
                        
                        elif source=="ncapture":
                            ncapture_E_T.append((step.GetGlobalTime(),step.GetTotalEdep()))
                            if not found:
                                if info[0] not in ncap_info:
                                    ncap_info.append(info[0])
                                if info[1]!=0:argon47=True
                                if info[2]!=0:iron76=True
                                
                        elif source=="neutron_inelastics": 
                            Neutron_inelastics_E_T.append((step.GetGlobalTime(),step.GetTotalEdep()))
                            if not found:
                                if info not in inelastic_info:
                                    inelastic_info.append(info)
                        
                        found=True
                            
                    elif step.GetTotalEdep()>0 and (step.GetGlobalTime(),step.GetTotalEdep()) in trigger:
                        if not found:
                            source,info=event_finder(track)
                                        
                        if source=="Source_gamma":
                            Source_E_T.append((step.GetGlobalTime(),step.GetTotalEdep()))
                            
                        elif source=="ncapture":
                            ncapture_E_T.append((step.GetGlobalTime(),step.GetTotalEdep()))
                            if not found:
                                if info[0] not in ncap_info:
                                    ncap_info.append(info[0])
                                if info[1]!=0:argon47=True
                                if info[2]!=0:iron76=True
                            
                        elif source=="neutron_inelastics": 
                            Neutron_inelastics_E_T.append((step.GetGlobalTime(),step.GetTotalEdep()))
                            if info not in inelastic_info:
                                if not found:
                                    inelastic_info.append(info)
                        found=True
        
        
        #DEAL WITH EVERYTHING
        Fprompt=add_energy(fprompt_e)/Energy
        neutron_start, gamma_start=start_conditions()
        
        
        #print(Energy,Fprompt,neutron_start,gamma_start,ncap_info,inelastic_info,argon47,iron76)
        fill_tree(Energy,Fprompt,neutron_start,gamma_start,Neutron_inelastics_E_T,ncapture_E_T,Source_E_T,ncap_info,inelastic_info,argon47,iron76)



    #perform code
    triggers=findevents()
    #print(triggers)
    for event in triggers:
        event_description(event)

tree.write()
ntp.close()


print("My program took", time.time() - start_time, "to run")




