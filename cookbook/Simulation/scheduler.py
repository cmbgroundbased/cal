import os
import numpy as np
import datetime
from apscheduler.schedulers.blocking import BlockingScheduler
from subprocess import Popen, PIPE
from argparse import ArgumentParser, RawTextHelpFormatter
import logging
import time

DEBUG = True

def check_job(N, ky):

    p = Popen("squeue -u stefano.mandelli", stdout=PIPE, shell=True)
    out = (str(p.communicate()[0]).split("\\n"))
    Njobs = int(np.shape(out)[0] - 2)
    
    if Njobs < N:
        os.system('clear')
        print("Valore di ky: {}".format(ky))
        t0 = datetime.datetime(2022, 1, 1, 0, 0, 0).timestamp()      
        
        t1 = 0
        t2 = 0
        t3 = 0
        
        t1 += t0 + ky*3600 # Now
        t1 = datetime.datetime.fromtimestamp(t1)
        y = t1.year
        m = t1.month
        d = t1.day
        h = t1.hour
        minu = t1.minute
        seco = t1.second

        t2 += t0 + (ky+1)*3600 # avanti di un'ora
        t2 = datetime.datetime.fromtimestamp(t2)
        y2 = t2.year
        m2 = t2.month
        d2 = t2.day
        h2 = t2.hour
        minu2 = t2.minute
        seco2 = t2.second
        
        t3 += t0 + (ky+2)*3600 # avanti di due ore
        t3 = datetime.datetime.fromtimestamp(t3)
        y3 = t3.year
        m3 = t3.month
        d3 = t3.day
        h3 = t3.hour
        minu3 = t3.minute
        seco3 = t3.second

        data_string =str(y)+","+str(m)+","+str(d)+","+str(h)+","+str(minu)+","+str(seco)
        data_string_2 = str(y2)+","+str(m2)+","+str(d2)+","+str(h2)+","+str(minu2)+","+str(seco2)
        data_string_3 = str(y3)+","+str(m3)+","+str(d3)+","+str(h3)+","+str(minu3)+","+str(seco3)

        Popen("sed \"4s/2022,1,1,0,0,0/"+data_string_2+"/\" strip_file_0.par > out.par", stdout=PIPE, shell=True)
        time.sleep(2)
        Popen("sed \"6s/2022,1,1,1,0,0/"+data_string_3+"/\" out.par > par_files/strip_file_"+str(ky)+".par", stdout=PIPE, shell=True)
        time.sleep(2)
        Popen("sed \"s/strip_file_0.par/..\/par_files\/strip_file_"+str(ky)+".par/\" strip_simulation_0.sl > slurm_files/strip_simulation_"+str(ky)+".sl", stdout=PIPE, shell=True)
        # Popen("rm -rf atm_cache*", stdout=PIPE, shell=True)
        

        

        print("STATUS: SUBMITTING - SUBMITTED JOBS : {}".format(Njobs))
        if DEBUG:
            sbatch = Popen("echo ciao", stdout=PIPE, shell=True)
            sbatch_out = str(sbatch.communicate()[0]).split("\\n")
        else:
            sbatch = Popen("sbatch slurm_files\/strip_simulation_"+str(ky)+".sl; sleep 2", stdout=PIPE, shell=True)
            sbatch_out = str(sbatch.communicate()[0]).split("\\n")
        
        print("DATE IN SUMISSION: from: {} to: {}".format(data_string_2, data_string_3))
    	# scrivere qui dentro il file di parametri!!!!    
        ky += 1
    else:
        os.system('clear')
        print("STATUS: WAITING FOR RESOURCE")
        print("SUBMITTED JOBS: {}".format(Njobs))
        print("")
    
    return ky

if __name__ == "__main__":
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter, description="Scheduler parameters")
    parser.add_argument("-num_jobs", default=0, type=int)
    
    args = parser.parse_args()
    N = args.num_jobs
    
    ky = 0
    
    while(True):
        ky = check_job(N, ky)
        time.sleep(4)
    
    # scheduler = BlockingScheduler()
    # scheduler.add_executor('processpool')

    # scheduler.add_job(check_job, 'interval',args=[N], seconds=5)

    # try:
    #     p=scheduler.start()
    # except(KeyboardInterrupt, SystemExit):
    #     pass
