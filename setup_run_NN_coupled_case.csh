#!/bin/csh
setenv CCSMTAG CESM2.1.2-liran
setenv CCSMROOT "/home1/00993/tg802402/repositories/$CCSMTAG"
setenv CASE        NNrealgeographyonline001
setenv CASESET     HIST_CAM%SPCAMS_CLM50%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV
setenv CASERES     f19_f19_mg17
setenv PROJECT     TG-ATM190002
setenv MACH      stampede2-knl
setenv CASEROOT  $HOME/CESM2_case/$CASE
setenv PTMP      $SCRATCH/
setenv RUNDIR    $PTMP/$CASE/run
setenv ARCHDIR   $PTMP/archive/$CASE
setenv DATADIR   $SCRATCH/SPCAM_input
cd $CCSMROOT/cime/scripts
rm -rf $CASEROOT
rm -rf $PTMP/$CASE
./create_newcase --case $CASEROOT --pecount 1024 --pesfile ./pelayout_01.xml --res f19_f19_mg17 --machine stampede2-knl --compset HIST_CAM%SPCAMS_CLM50%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV --input-dir /scratch/07088/tg863871/CESM_inputdata --output-root $SCRATCH --run-unsupported
cd  $CASEROOT
xmlchange --file env_batch.xml --id JOB_QUEUE --val normal
xmlchange --file env_workflow.xml --id JOB_WALLCLOCK_TIME --val 00:30:00
xmlchange --file env_run.xml --id RUN_STARTDATE --val 2004-01-01
xmlchange --file env_run.xml --id STOP_N --val 48
xmlchange --file env_run.xml --id run_data_archive --val "FALSE"

#Want an atmospheric IC in approx equil with actual SP so mining a 200401 IC file from the training data run...
cat <<EOF >! user_nl_cam
npr_yz = 32,4,4,32
fincl2 = 'NN2L_DSTDRY4:A','NN2L_DSTWET4:A','NN2L_DSTDRY3:A','NN2L_DSTWET3:A','NN2L_DSTDRY2:A','NN2L_DSTWET2:A','NN2L_DSTDRY1:A','NN2L_DSTWET1:A','NN2L_OCPHODRY:A','NN2L_OCPHIDRY:A','NN2L_OCPHIWET:A','NN2L_BCPHODRY:A','NN2L_BCPHIDRY:A','NN2L_BCPHIWET:A','NN2L_PSL:A','NN2L_CO2DIAG:I','NN2L_CO2PROG:I','NN2L_THBOT:I','NN2L_SOLSD:I','NN2L_SOLLD:I','NN2L_SOLS:I','NN2L_SOLL:I','NN2L_PRECL:A','NN2L_PRECC:A','NN2L_PRECSL:A','NN2L_PRECSC:A','NN2L_FLWDS:I','NN2L_NETSW:I','NN2L_TBOT:A','NN2L_ZBOT:A','NN2L_UBOT:A','NN2L_VBOT:A','NN2L_QBOT:A','NN2L_PBOT:I','NN2L_RHO:I','T:I','Q:I','QAP:I','QBP:I','TBP:I','VAP:I','PS:I','SHFLX:I','LHFLX:I','PTTEND:I','PTEQ:I','SOLIN:I','FSNT:I','FSNS:I','FLNT:I','FLNS:I','U10:I','FLDS:I','FSDS:I','PRECT:I','CLDLIQAP:I','CLDICEAP:I','PTECLDLIQ:I','PTECLDICE:I','SPLIQICESTORAGE:I','SPICESTORAGE:I','SPTKE:I','SPTKES:I','SPQTFLX:I','SPQTFLXS:I','SPQPFLX:I','SPQPFLX:I','SPMC:I','SPMCUP:I','SPMCDN:I','PRECC:I','SPQPFLX:I','QRL:I','QRS:I','QRLC:I','QRSC:I','SPQC:I','SPQI:I','CLOUD:I','TBC:A', 'QBC:A', 'CLDLIQBC', 'CLDICEBC','RELHUM:I','CLDLIQBP','CLDICEBP'
fincl3 = 'NN_BP2BC_DT','NN_BP2BC_DQ','NN_BP2BC_DQC','NN_BP2BC_DQI','SP_BP2BC_DT','SP_BP2BC_DQ','SP_BP2BC_DQC','SP_BP2BC_DQI' 
nhtfrq = 0,1,1
mfilt  = 0,16,1
ncdata = '/scratch/07088/tg863871/CESM2_f19_v13_updated_NN_pelayout01_ens_07/run/CESM2_f19_v13_updated_NN_pelayout01_ens_07.cam.i.2004-01-01-00000.nc'
/
EOF

#... unfortunately land ICs were not saved from the training data run so first tried wiring finidat to closest available item, a restart file a couple weeks prior.
#finidat = '/scratch/07088/tg863871/CESM2_f19_v13_updated_NN_pelayout01_ens_07/run/CESM2_f19_v13_updated_NN_pelayout01_ens_07.clm2.r.2003-12-18-00000.nc'
#... but this failed with errors.
#recourse = leave finidat empty, accept absurdly long init time as model auto-interpolates the following default land IC:
# /scratch/07088/tg863871/CESM_inputdata/lnd/clm2/initdata_map/clmi.BHIST.2000-01-01.0.9x1.25_gx1v7_simyr2000_c181015.nc
# But we save the answer here to avoid lengthy init for future tests:
cat <<EOF >! user_nl_clm
&clmexp
/
EOF

#! Copy in git code mods for coupled NN online and associated diagnostics:
#! diffs up to date:
cp $HOME/repositories/CESM2-neuralnet-UCI/*.F90 ./SourceMods/src.cam
./case.setup
./case.build
#./case.submit
cd ..
end
