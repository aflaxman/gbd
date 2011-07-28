### some notes on deleting covariates:
import glob

fnames = glob.glob('/home/j/Project/Causes of Death/Covariate Search/Covariates/DISMOD/datasets/INCORRECT COVARIATES THAT HAVE BEEN FIXED & REUPLOADED_DO NOT USE THESE/*')

slugs = [f.split('/')[-1].replace('.csv','') for f in fnames]

for s in slugs:
    try:
        ct = models.CovariateType.objects.get(slug=s)
        ct.delete()
        print s, 'deleted'
    except:
        print s, 'not found'
