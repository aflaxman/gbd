
# over-ride these in the local_settings file appropriately

DISMOD_BASE_URL = 'http://127.0.0.1:8000/'


DISMOD_USERNAME = 'cjm'
DISMOD_PASSWORD = 'kuma'

DISMOD_TWITTER_NAME = 'ihme_dismod'
DISMOD_TWITTER_PASSWORD = 's3cr3t_p4sswd'


DISMOD_LOGIN_URL = DISMOD_BASE_URL + 'accounts/login/'
DISMOD_DOWNLOAD_URL = DISMOD_BASE_URL + 'dismod/show/%s.json'
DISMOD_UPLOAD_URL = DISMOD_BASE_URL + 'dismod/upload'

DISMOD_LIST_JOBS_URL = DISMOD_BASE_URL + 'dismod/job_queue/list/?format=json'
DISMOD_REMOVE_JOB_URL = DISMOD_BASE_URL + 'dismod/job_queue/remove/'
DISMOD_INIT_LOG_URL = DISMOD_BASE_URL + 'dismod/init_log/%d/%s'
DISMOD_LOG_STATUS_URL = DISMOD_BASE_URL + 'dismod/log_status/%d/%s/%s/%s'


ON_SGE = 0
SERVER_LOAD_STATUS_HOST = '128.208.10.133'
SERVER_LOAD_STATUS_PORT = 1723
SERVER_LOAD_STATUS_SIZE = 20480

# path to job working directory
# dir = JOB_WORKING_DIR % id
#JOB_WORKING_DIR = '../../dismod_status/test/dm-%d'
JOB_WORKING_DIR = '/tmp/dismod_working/test/dm-%d'

# path to job log directory
# dir = JOB_LOG_DIR % id
JOB_LOG_DIR = '/tmp/dismod_log/test/dm-%d'

# path and name of daemon log file
DAEMON_LOG_FILE = '/tmp/daemon_test.log'

# path and name of gbd_fit lock file
GBD_FIT_LOCK_FILE = '/tmp/gbd_fit_test.lock'

# shell command string to spawn a fit process
GBD_FIT_STR = 'python gbd_fit.py %s %d >%s 2>%s'
#GBD_FIT_STR = 'qsub -cwd -o %s -e %s /home/j/dismod/gbd/gbd_fit.sh %s %d'

# time to wait (in seconds) between checking the server for new jobs
SLEEP_SECS = 2.


CSV_PATH = './'

# disease model parameters
NEARLY_ZERO = 1.e-10
MAX_AGE = 101

MISSING = -99

PRIOR_SEP_STR = ','

KEY_DELIM_CHAR = '+'

data_types = ['prevalence data',
              'incidence data',
              'remission data',
              'excess-mortality data',
              'duration data',
              'all-cause mortality data',
              ]

output_data_types = ['Prevalence',
                     'Incidence',
                     'Remission',
                     'Excess-mortality',
                     'Mortality',
                     'Relative-risk',
                     'Duration',
                     'Incidence x Duration']

stoch_var_types = output_data_types + ['bins']

gbd_regions = [u'Asia Pacific, High Income',
               u'Asia, Central',
               u'Asia, East',
               u'Asia, South',
               u'Asia, Southeast',
               u'Australasia',
               u'Caribbean',
               u'Europe, Central',
               u'Europe, Eastern',
               u'Europe, Western',
               u'Latin America, Andean',
               u'Latin America, Central',
               u'Latin America, Southern',
               u'Latin America, Tropical',
               u'North Africa/Middle East',
               u'North America, High Income',
               u'Oceania',
               u'Sub-Saharan Africa, Central',
               u'Sub-Saharan Africa, East',
               u'Sub-Saharan Africa, Southern',
               u'Sub-Saharan Africa, West']

gbd_years = ['1990', '2005']

gbd_sexes = ['Male', 'Female']

gbd_ages = [0, 1, 5, 10, 15, 20, 25, 35, 45, 55, 65, 75, 85, 100]


try:
    from local_settings import *
except:
    pass
