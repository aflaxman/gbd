""" script to produce plot of dismod jobs by week with Python/Pandas"""

import glob
import os
import time

import pylab as pl
import pandas

fnames = sorted(glob.glob('/home/j/Project/dismod/dismod_status/prod/dm-*'))
stats = [os.stat(f) for f in fnames]

dm_jobs = pandas.DataFrame({'year-week': [time.strftime('%Y-%W', time.localtime(s.st_ctime)) for s in stats]}, index=fnames)
dm_jobs['one'] = 1
jobs_per_week = dm_jobs.groupby('year-week')['one'].sum()

pl.clf()
pl.plot(jobs_per_week.sort_index().ix[6:].__array__())
pl.xticks(range(len(jobs_per_week.sort_index().ix[6:])), list(jobs_per_week.sort_index().ix[6:].index), rotation=90)
pl.subplots_adjust(bottom=.2)
pl.title('DisMod jobs run per week')
