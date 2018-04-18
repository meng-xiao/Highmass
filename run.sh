workdir=/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/Efficiency/HighMassCombInputs/prepare
#cd $CMSSW_BASE/src
#eval `scramv1 runtime -sh`
cd $workdir

#for w in 0.1 #10.0 100.0 
#do
#	for i in {0..17}  
#	do
#		m=$((130+i*10))
#			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 0 1
#			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 1 1
#			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 0 0
#			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 1 0
#	done
#done

#

#for w in 0.1 #10.0 0.1 100.0 
#do
#	for m in 1100 1250 1500 2000 2500 3000 
#	do
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 0 1
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 1 1
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 0 0
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 1 0
##
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 0 1
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 1 1
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 0 0
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 1 0
#	done
#done

#
#for w in 0.1 #10.0 100.0 
#do
#	for i in {0..13}
#	#for i in 4 
#	do
#		m=$((350+i*50))
#		if [ $m -ge 550 ]; then
#			echo "2l2q $m"
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 0 1
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 1 1
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 0 0
#			bsub -C 0 -q cmscaf1nd job_2l2q.lsf $m $w 1 0
#		fi;
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 0 1
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 1 1
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 0 0
##			bsub -C 0 -q cmscaf1nd job_4l.lsf $m $w 1 0
#	done
#done
