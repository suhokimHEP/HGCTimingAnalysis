#!/usr/bin/python                                                                                                                                                                                                                 

import os
currentDir = os.getcwd()
print 'currentDir = %s' %(currentDir)

##settings

#PDGID = 22
#nRUN = ["PDG_22_pt5_200PU_LE"]
#nRUN.append("PDG_22_pt5_200PU_HE")

#PDGID = 211
#nRUN = ["PDG_211_pt5_200PU_LE"]
#nRUN = ["PDG_211_pt5_200PU_HE"]

PDGID = 130
nRUN = ["testTimeCalib1p2Side_NuGun"]
#nRUN.append("PDG_130_pt5_200PU_LE")

#nRUN = ["PDG_22_pt10_0PU_allEta"]
#nRUN = ["PDG_130_pt2_0PU_allEta"]
#nRUN.append("PDG_130_pt5_0PU_allEta")
#nRUN.append("PDG_130_pt10_0PU_allEta")


print 'number of jobs = %d' %(len(nRUN))

PGENPT = 5.



#outLancia = open('%s/lancia.sh' %(currentDir),'w');

cfgTemp = "part_Time_template.py"

for num in nRUN:
    print num;

    outDir = "%s" %(num)
    subDir1 = currentDir+"/"+outDir;
    os.system('mkdir %s' %(subDir1));

    subDir = currentDir+"/"+outDir+"/JOB_%s" %(num);
    os.system('mkdir %s' %(subDir));
    print 'subDir = %s  >>>> DID YOU CREATE THIS?' %(subDir);

    #inFileList = currentDir+"/listTTBar_TTbar.txt" ;
    inFileList = currentDir+"/listNuGun.txt" ;
    print 'inFileList = %s' %(inFileList);

    cfgFile = "partGun_NTUP_fromtemp_%s.py" %(num);
    os.system("cp %s %s/%s" %(cfgTemp, subDir, cfgFile));
    os.system("sed -i s~PGENPT~%s~g %s" %(PGENPT, subDir+"/"+cfgFile));
    os.system("sed -i s~CALOPARTPDGID~%s~g %s" %(PDGID, subDir+"/"+cfgFile));
    os.system("sed -i s~INPUTFILELIST~%s~g %s" %(inFileList, subDir+"/"+cfgFile));
    os.system("sed -i s~OUTFILE~%s~g %s" %("OutTimeHGC_RecHitsCalib_"+num, subDir+"/"+cfgFile));

    outScript = open('%s/bjob.sh' %(subDir),'w');
    outScript.write('#!/bin/bash \n');
    outScript.write('cd %s \n' %(subDir));
    outScript.write('export SCRAM_ARCH=slc7_amd64_gcc700 \n');
    outScript.write('eval `scramv1 ru -sh` \n');
    outScript.write('pwd \n ')
    outScript.write('cmsRun  %s \n' %(cfgFile) );
    outScript.write( '\n ' );
    os.system('chmod 777 %s/bjob.sh' %(subDir));

    # outLancia.write(' bsub -cwd %s -q 8nh %s/bjob.sh \n' %(subDir, subDir));

    write_condorjob = open(currentDir+'/lanciaCondor_'+num+'.sub', 'w')
    write_condorjob.write('Universe = vanilla \n')
    write_condorjob.write('Executable  = '+subDir+'/bjob.sh \n')
    #write_condorjob.write('+JobFlavour = "'+opt.QUEUE+'" \n\n')
    #write_condorjob.write('+JobFlavour = tomorrow \n\n')
    write_condorjob.write('+MaxRuntime = 43200 \n\n')
    #write_condorjob.write('arguments   = $(ClusterID) $(ProcId) '+currentDir+' '+outDir+' '+cfgfile+' '+str(opt.LOCAL)+' '+CMSSW_VERSION+' '+CMSSW_BASE+' '+SCRAM_ARCH+' '+opt.eosArea+' '+opt.DTIER+' '+str(opt.DQM)+'\n')
    write_condorjob.write('requirements = (OpSysAndVer =?= \"CentOS7\") \n')
    write_condorjob.write('transfer_output_files = "" \n')
    write_condorjob.write('output      = '+subDir+'/output.out \n')
    write_condorjob.write('error       = '+subDir+'/error.err \n')
    write_condorjob.write('log         = '+subDir+'/log_htc.log \n\n')
    write_condorjob.write('request_cpus = 4 \n')
    write_condorjob.write('max_retries = 1\n')
    write_condorjob.write('queue JobID in ('+num+') \n')
    write_condorjob.close()
    
    
