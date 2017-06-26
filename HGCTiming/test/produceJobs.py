#!/usr/bin/python                                                                                                                                                                                                                 

#nRUN = 245192
#nRUN = "MinBias"
## change run number

import os
currentDir = os.getcwd()
print 'currentDir = %s' %(currentDir)

##settings
#nRUN = ["PDG_130_Pt07_3fC"]
#nRUN.append("PDG_130_Pt1_3fC")
#nRUN.append("PDG_130_Pt2_3fC")
#nRUN.append("PDG_130_Pt5_3fC")
#nRUN.append("PDG_130_Pt10_3fC")
#nRUN.append("PDG_130_Pt30_3fC")
#nRUN.append("PDG_130_Pt100_3fC")
nRUN = ["PDG_22_Pt2_3fC"] 
nRUN.append("PDG_22_Pt5_3fC")
nRUN.append("PDG_22_Pt60_3fC")

print 'number of jobs = %d' %(len(nRUN))

CFDValue = 1

CELLT = 0            
FLOORV = 0.02        
LIFEA = 0            
ABSTREND = 2.        
outDir = "CSF20LBA20"

#CELLT = 0            
#FLOORV = 0.02        
#LIFEA = 1            
#ABSTREND = 2.        
#outDir = "CSF20LEA20"

######

#CELLT = 0
#FLOORV = 0.03
#LIFEA = 0
#ABSTREND = 1.
#outDir = "CSF30LBA10"

#CELLT = 0
#FLOORV = 0.03
#LIFEA = 0
#ABSTREND = 1.5 
#outDir = "CSF30LBA15"

#CELLT = 0
#FLOORV = 0.03
#LIFEA = 0
#ABSTREND = 2. 
#outDir = "CSF30LBA20"


#CELLT = 0
#FLOORV = 0.03
#LIFEA = 1
#ABSTREND = 1.
#outDir = "CSF30LEA10"

#CELLT = 0
#FLOORV = 0.03
#LIFEA = 1
#ABSTREND = 1.5 
#outDir = "CSF30LEA15"

#CELLT = 0
#FLOORV = 0.03
#LIFEA = 1
#ABSTREND = 2. 
#outDir = "CSF30LEA20"





subDir1 = currentDir+"/"+outDir;
os.system('mkdir %s' %(subDir1));


outLancia = open('%s/lancia.sh' %(currentDir),'w');

cfgTemp = "part_Time_template.py"

for num in nRUN:
    print num;
    subDir = currentDir+"/"+outDir+"/JOB_%s" %(num);
    os.system('mkdir %s' %(subDir));
    print 'subDir = %s  >>>> DID YOU CREATE THIS?' %(subDir);

    inFileList = currentDir+"/fileList/%s.txt" %(num);
    print 'inFileList = %s' %(inFileList);

    cfgFile = "partGun_NTUP_fromtemp_%s.py" %(num);
    os.system("cp %s %s/%s" %(cfgTemp, subDir, cfgFile));
    os.system("sed -i s~CFDVal~%s~g %s" %(CFDValue, subDir+"/"+cfgFile));
    os.system("sed -i s~CELLT~%s~g %s" %(CELLT, subDir+"/"+cfgFile));
    os.system("sed -i s~FLOORV~%s~g %s" %(FLOORV, subDir+"/"+cfgFile));
    os.system("sed -i s~LIFEA~%s~g %s" %(LIFEA, subDir+"/"+cfgFile));
    os.system("sed -i s~ABSTREND~%s~g %s" %(ABSTREND, subDir+"/"+cfgFile));
    os.system("sed -i s~INPUTFILELIST~%s~g %s" %(currentDir+"/fileList/"+num+".txt", subDir+"/"+cfgFile));
    os.system("sed -i s~OUTFILE~%s~g %s" %("OutTimeHGC_RecHits_"+num, subDir+"/"+cfgFile));

    outScript = open('%s/bjob.sh' %(subDir),'w');
    outScript.write('#!/bin/bash \n');
    outScript.write('cd %s \n' %(subDir));
    outScript.write('export SCRAM_ARCH=slc6_amd64_gcc530 \n');
    outScript.write('eval `scramv1 ru -sh` \n');
    #outScript.write('cd - \n');
    #outScript.write('cmsMkdir %s \n' %(outFolder));
    outScript.write('pwd \n ')
    outScript.write('cmsRun  %s \n' %(cfgFile) );
    #outScript.write('cmsStage coll_timing_%d_%d.root %s/' %(num, num+plusNum, outFolder));
    outScript.write( '\n ' );
    os.system('chmod 777 %s/bjob.sh' %(subDir));
    outLancia.write(' bsub -cwd %s -q 8nh %s/bjob.sh \n' %(subDir, subDir));

