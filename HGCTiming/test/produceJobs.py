#!/usr/bin/python                                                                                                                                                                                                                 

import os
currentDir = os.getcwd()
print 'currentDir = %s' %(currentDir)

##settings

nRUN = ["PDG_22_pt2_0PU_allEta"]
#nRUN = ["PDG_22_pt5_0PU_allEta"]
#nRUN = ["PDG_22_pt10_0PU_allEta"]
#nRUN = ["PDG_130_pt2_0PU_allEta"]
#nRUN.append("PDG_130_pt5_0PU_allEta")
#nRUN.append("PDG_130_pt10_0PU_allEta")


print 'number of jobs = %d' %(len(nRUN))

PGENPT = 2.
CFDValue = 22
TIMEOFF = 3


outLancia = open('%s/lancia.sh' %(currentDir),'w');

cfgTemp = "part_Time_template.py"

for num in nRUN:
    print num;

    outDir = "%s" %(num)
    subDir1 = currentDir+"/"+outDir;
    os.system('mkdir %s' %(subDir1));

    subDir = currentDir+"/"+outDir+"/JOB_%s" %(num);
    os.system('mkdir %s' %(subDir));
    print 'subDir = %s  >>>> DID YOU CREATE THIS?' %(subDir);

    inFileList = currentDir+"/fileListAllEta/%s.txt" %(num);
    print 'inFileList = %s' %(inFileList);

    cfgFile = "partGun_NTUP_fromtemp_%s.py" %(num);
    os.system("cp %s %s/%s" %(cfgTemp, subDir, cfgFile));
    os.system("sed -i s~PGENPT~%s~g %s" %(PGENPT, subDir+"/"+cfgFile));
    os.system("sed -i s~CFDVal~%s~g %s" %(CFDValue, subDir+"/"+cfgFile));
    os.system("sed -i s~TIMEOFF~%s~g %s" %(TIMEOFF, subDir+"/"+cfgFile));
    os.system("sed -i s~INPUTFILELIST~%s~g %s" %(currentDir+"/fileListAllEta/"+num+".txt", subDir+"/"+cfgFile));
    os.system("sed -i s~OUTFILE~%s~g %s" %("OutTimeHGC_RecHits_"+num, subDir+"/"+cfgFile));

    outScript = open('%s/bjob.sh' %(subDir),'w');
    outScript.write('#!/bin/bash \n');
    outScript.write('cd %s \n' %(subDir));
    outScript.write('export SCRAM_ARCH=slc6_amd64_gcc630 \n');
    outScript.write('eval `scramv1 ru -sh` \n');
    outScript.write('pwd \n ')
    outScript.write('cmsRun  %s \n' %(cfgFile) );
    outScript.write( '\n ' );
    os.system('chmod 777 %s/bjob.sh' %(subDir));
    outLancia.write(' bsub -cwd %s -q 8nh %s/bjob.sh \n' %(subDir, subDir));

