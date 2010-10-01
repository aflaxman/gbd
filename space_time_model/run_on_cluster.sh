# Use the bash shell to interpret this job script
#$ -S /bin/bash
#

# submit this job to nodes that have
# at least 1GB of RAM free.
#$ -l mem_free=1.0G


## Put the hostname, current directory, and start date
## into variables, then write them to standard output.
cd /home/OUTPOST/abie/pymc-space-time-model/


GSITSHOST=`/bin/hostname`
GSITSPWD=`/bin/pwd`
GSITSDATE=`/bin/date`
echo "**** JOB STARTED ON $GSITSHOST AT $GSITSDATE"
echo "**** JOB RUNNING IN $GSITSPWD"
##


echo calling python -u "$@"
/usr/local/epd_py25-4.3.0/bin/python -u "$@"

