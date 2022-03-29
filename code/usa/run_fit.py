#!/usr/bin/python

#
#  run_fit.py - driver for running the fitting codes in parallel 
# by launching multiplt (11 in this case) jobs
#  
#
#  Created by Michal Ben-Nun on 5/19/20.
#  
import multiprocessing
from multiprocessing import Process
import os
import string
import random
import numpy as np


def info(title):
    print title
    print 'module name:',__name__
if hasattr(os,'getppid'):
    print 'parent process:',os.getpid()
    print 'process id: ',os.getpid()

##
## change which R fitting routine to use 
## usa_mcmc_fit_calendar_time.R or
## usa_mcmc_fit_pandemic_time.R
##

##
## start/end: currently set so all locations in the US are fitted using 11 threads
##

def f(start,end):
    os.system("time Rscript usa_mcmc_fit_calendar_time.R "+str(start)+" "+str(end))

nthreads = 11

start=np.arange(1,54,5)
end=np.arange(5,60,5)

if __name__ == '__main__':
    info('main line')
    p1 = Process(target=f,args=(start[0],end[0],))
    p1.start()
    p2 = Process(target=f,args=(start[1],end[1],))
    p2.start()
    p3 = Process(target=f,args=(start[2],end[2],))
    p3.start()
    p4 = Process(target=f,args=(start[3],end[3],))
    p4.start()
    p5 = Process(target=f,args=(start[4],end[4],))
    p5.start()
    p6 = Process(target=f,args=(start[5],end[5],))
    p6.start()
    p7 = Process(target=f,args=(start[6],end[6],))
    p7.start()
    p8 = Process(target=f,args=(start[7],end[7],))
    p8.start()
    p9 = Process(target=f,args=(start[8],end[8],))
    p9.start()
    p10 = Process(target=f,args=(start[9],end[9],))
    p10.start()
    p11 = Process(target=f,args=(start[10],end[10],))
    p11.start()

