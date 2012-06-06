# Use the bash shell to interpret this job script
#$ -S /bin/bash

## Put the hostname, current directory, and start date
## into variables, then write them to standard output.
GSITSHOST=`/bin/hostname`
GSITSPWD=`/bin/pwd`
GSITSDATE=`/bin/date`
echo "**** JOB STARTED ON $GSITSHOST AT $GSITSDATE"
echo "**** JOB RUNNING IN $GSITSPWD"
##

cd /net/gs/vol1/home/abie/omak_gbd

/usr/local/epd_py25-4.3.0/bin/python tests/simulation_study.py $@

# unrelated notes
# for s in male female; do for y in 1990 2005; do qsub gbd_fit.sh 45 -s $s -y $y -r north_america_high_income -N -P 'confidence 500 .01, zero 0 5, smooth 100.0' -I 'confidence 500 .01, zero 0 5 smooth 100.0' -R 'confidence 500 .01, smooth 100.0' -C 'zero 0 100'; done; done

