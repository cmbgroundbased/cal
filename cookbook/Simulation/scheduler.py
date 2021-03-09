import os
import numpy as np
import datetime
from apscheduler.schedulers.blocking import BlockingScheduler
from subprocess import Popen, PIPE

year   = 2022
month  = 1
day    = 1
hour   = 0
prime  = 0
second = 0

t0 = datetime.datetime(year, month, day, hour, prime, second).timestamp()

def check_job():
	p = Popen("squeue -u stefano.mandelli", stdout=PIPE, shell=True)
	out = str(p.communicate()[0]).split("\\n")
	if (np.shape(out)[0]-2 <= 6):
		print("Lancia job - N job attivi: {}".format(np.shape(out)[0]-2))
        
        p = datetime.datetime.fromtimestamp(t0)
        
        y = p.year
        m = p.month
        d = p.day
        h = p.hour
        min = p.minute
        sec = p.second
        
        data_string =str(y)+","+str(m)+","+str(d)+","+str(h)+","+str(min)+","+str(sec)
                
        # cambiare il file dei parametri
        
		sbatch = Popen("sbatch toast_batch.sl", stdout=PIPE, shell=True)
		sbatch_out = str(sbatch.communicate()[0]).split("\\n")
		so = ""
		for i in sbatch_out:
			so += i 
		print("INFO: {}".format(so))
		print("")
        t0 += 3600 # avanti di un'ora.
		
	else:
		print("Wait")
		print("Job attivi: {}".format(np.shape(out)[0]-2))
		print("")

scheduler = BlockingScheduler()
scheduler.add_executor('processpool')

scheduler.add_job(check_job, 'interval', seconds=3)

try:
	scheduler.start()
except (KeyboardInterrupt):
	pass
