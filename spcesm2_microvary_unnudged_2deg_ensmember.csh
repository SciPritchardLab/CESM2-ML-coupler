#!/bin/csh
#v02 -- changed microphys setting broadcast to is_first_step or is_first_restart_stop (important)
#foreach ens ( 02 03 ) #04 05 06 )
#foreach ens ( 07 08 09 10 11 12 )
foreach ens ( 07 )
set casename = "CESM2_unnudge_f19_updated_NN_v2_pelayout01_microvary_5param_v02_ens_${ens}"
./create_newcase --case $casename --pecount 1024 --pesfile ./pelayout_01.xml --res f19_f19_mg17 --machine stampede2-knl --compset HIST_CAM%SPCAMS_CLM50%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV --input-dir /scratch/07088/tg863871/CESM_inputdata --output-root /scratch/07088/tg863871 --run-unsupported

cd  $casename

xmlchange --file env_batch.xml --id JOB_QUEUE --val normal
xmlchange --file env_batch.xml --id JOB_WALLCLOCK_TIME --val 01:59:00
xmlchange --file env_run.xml --id RUN_STARTDATE --val 2003-01-02
xmlchange --file env_run.xml --id STOP_N --val 50
#xmlchange --file env_run.xml --id RESUBMIT --val 4

cat <<EOF >! user_nl_cam
npr_yz = 32,4,4,32
fincl2 = 'NN2L_DSTDRY4:I','NN2L_DSTWET4:I','NN2L_DSTDRY3:I','NN2L_DSTWET3:I','NN2L_DSTDRY2:I','NN2L_DSTWET2:I','NN2L_DSTDRY1:I','NN2L_DSTWET1:I','NN2L_OCPHODRY:I','NN2L_OCPHIDRY:I','NN2L_OCPHIWET:I','NN2L_BCPHODRY:I','NN2L_BCPHIDRY:I','NN2L_BCPHIWET:I','NN2L_PSL:I','NN2L_CO2DIAG:I','NN2L_CO2PROG:I','NN2L_THBOT:I','NN2L_SOLSD:I','NN2L_SOLLD:I','NN2L_SOLS:I','NN2L_SOLL:I','NN2L_PRECL:I','NN2L_PRECC:I','NN2L_PRECSL:I','NN2L_PRECSC:I','NN2L_FLWDS:I','NN2L_NETSW:I','NN2L_TBOT:I','NN2L_ZBOT:I','NN2L_UBOT:I','NN2L_VBOT:I','NN2L_QBOT:I','NN2L_PBOT:I','NN2L_RHO:I','T:I','Q:I','QAP:I','QBP:I','TBP:I','VAP:I','PS:I','SHFLX:I','LHFLX:I','PTTEND:I','PTEQ:I','SOLIN:I','FSNT:I','FSNS:I','FLNT:I','FLNS:I','U10:I','FLDS:I','FSDS:I','PRECT:I','CLDLIQAP:I','CLDICEAP:I','PTECLDLIQ:I','PTECLDICE:I','SPLIQICESTORAGE:I','SPICESTORAGE:I','SPTKE:I','SPTKES:I','SPQTFLX:I','SPQTFLXS:I','SPQPFLX:I','SPQPFLX:I','SPMC:I','SPMCUP:I','SPMCDN:I','PRECC:I','SPQPFLX:I','QRL:I','QRS:I','QRLC:I','QRSC:I','SPQC:I','SPQI:I','CLOUD:I','TBC:I', 'QBC:I', 'CLDLIQBC', 'CLDICEBC'
nhtfrq = 0,-24
mfilt = 0,1
/
EOF

cat <<EOF >! user_nl_clm
&clmexp
/
EOF

cp /home1/07088/tg863871/repositories/Liran_Mike_CESM2/*.F90 ./SourceMods/src.cam

./case.setup
./case.build
./case.submit 

cd ..
end
