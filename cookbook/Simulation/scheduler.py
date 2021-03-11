import os
import numpy as np
import datetime
from apscheduler.schedulers.blocking import BlockingScheduler
from subprocess import Popen, PIPE
from argparse import ArgumentParser, RawTextHelpFormatter

DEBUG = True

global ky
ky = 0

def check_job(N):
    global ky
    
    if DEBUG:
        p = Popen("ps aux | grep bash | wc -l", stdout=PIPE, shell=True)
        out = (str(p.communicate()[0]).split(","))
        out=str(out[0].split("\\n")[0].split("\'")[1])
        Njobs = int(out)
    else:
    	p = Popen("squeue -u stefano.mandelli", stdout=PIPE, shell=True)
    	out = (str(p.communicate()[0]).split("\\n"))
    	Njobs = int(np.shape(out)[0] - 2)
    
    if Njobs < N:
        os.system('clear')
        t0 = datetime.datetime(2022, 1, 1, 0, 0, 0).timestamp()
        t1 = 0
        t0 = datetime.datetime.fromtimestamp(t0)
        y = t0.year
        m = t0.month
        d = t0.day
        h = t0.hour
        minu = t0.minute
        seco = t0.second
        data_string_old =str(y)+","+str(m)+","+str(d)+","+str(h)+","+str(minu)+","+str(seco)
        
        print("STATUS: SUBMITTING - SUBMITTED JOBS : {}".format(Njobs))
        if DEBUG:
            sbatch = Popen("gnome-terminal", stdout=PIPE, shell=True)
            sbatch_out = str(sbatch.communicate()[0]).split("\\n")
        else:
            sbatch = Popen("sbatch toast_batch.sl", stdout=PIPE, shell=True)
            sbatch_out = str(sbatch.communicate()[0]).split("\\n")      
        
        t1 += t0.timestamp() + ky*3600 # avanti di un'ora
        t1 = datetime.datetime.fromtimestamp(t1)
        y = t1.year
        m = t1.month
        d = t1.day
        h = t1.hour
        minu = t1.minute
        seco = t1.second
        data_string =str(y)+","+str(m)+","+str(d)+","+str(h)+","+str(minu)+","+str(seco)
        Popen("/bin/sed \"s/"+data_string_old+"/"+data_string+"/\" toast_batch.sl > out.sl", stdout=PIPE, shell=True)
        
        print("DATE IN SUMISSION: {}".format(data_string))
        
        ky +=1
    else:
        os.system('clear')
        print("STATUS: WAITING FOR RESOURCE")
        print("SUBMITTED JOBS: {}".format(Njobs))
        print("")
        

if __name__ == "__main__":
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter, description="Scheduler parameters")
    parser.add_argument("-num_jobs", default=0, type=int)
    
    args = parser.parse_args()
    N = args.num_jobs
    
    scheduler = BlockingScheduler()
    scheduler.add_executor('processpool')

    scheduler.add_job(check_job, 'interval',args=[N], seconds=3)

    try:
        scheduler.start()
    except (KeyboardInterrupt):
        pass
