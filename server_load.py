import daemon
import SocketServer
import subprocess
import os
import time
import sys
import signal

SERVER_LOAD_STATUS_HOST = 'omak.ihme.washington.edu'
SERVER_LOAD_STATUS_PORT = 1723
DAEMON_LOG_FILE = '/tmp/daemon_test.log'
SERVER_LOAD_LOCK_FILE = '/tmp/server_load_test.lock'

def log(message):
    sys.stdout.write('%s  server_load  %s\n' % (time.strftime("%Y-%m-%d %H:%M:%S"), message))
    sys.stdout.flush()

def term(self, *args):
    log('server_load daemon received SIGTERM')
    sys.exit()

class RequestHandler(SocketServer.BaseRequestHandler ):
    def setup(self):
        #print self.client_address, 'connected!'
        s = subprocess.Popen(["qstat"], shell=True, stdout=subprocess.PIPE).communicate()[0]
        self.request.send(s)
        self.request.close()

    def handle(self):
        #while 1:
            #data = self.request.recv(1024)
            #s = subprocess.Popen(["qstat"], shell=True, stdout=subprocess.PIPE).communicate()[0]
            #self.request.send(s)
            #if data.strip() == 'bye':
                #return
        self.request.close()

    def finish(self):
        #print self.client_address, 'disconnected!'
        #self.request.send('bye ' + str(self.client_address) + '\n')
        self.request.close()

daemon.daemonize('/dev/null', DAEMON_LOG_FILE, DAEMON_LOG_FILE)
f = open(SERVER_LOAD_LOCK_FILE, 'w')
f.write(str(os.getpid()))
f.close()
signal.signal(signal.SIGTERM, term)
log('starting server_load daemon...')
server = SocketServer.ThreadingTCPServer((SERVER_LOAD_STATUS_HOST, SERVER_LOAD_STATUS_PORT), RequestHandler)
server.serve_forever()

