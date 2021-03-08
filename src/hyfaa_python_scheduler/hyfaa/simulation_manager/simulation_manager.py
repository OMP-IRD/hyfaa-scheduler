#!/usr/bin/env python
# -*- coding: utf-8 -*-

#######################################################################
#  This code has been developped by Magellium SAS
#
#  Licensing:
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>
#######################################################################

""" 
    Author: remi.jugier@magellium.fr
"""

import os, sys, subprocess, signal
if sys.version_info.major < 3:
    raise Exception('sorry, only python 3 supported')
from datetime import datetime, timedelta
import time
import multiprocessing
from threading import Thread, Event
import traceback
from queue import Queue, Empty


        

#############################################

def get_from_queue(queue_in, timeout=None):
    try:
        if timeout is not None:
            out = queue_in.get(block=True,timeout=timeout)
        else:
            out = queue_in.get(block=False)
    except Empty:
        out = None
    return out
        
def string2bytes(str_in):
    if python_version < 3:
        return str_in
    else:
        return str_in.encode('utf-8')
        
def bytes2string(str_in):
    if python_version < 3:
        return str_in
    else:
        return str_in.decode('utf-8')


###########################################
class SimulationTasks(object):
    """Creates processes, loads tasks in a workQueue and gets result in a result queue"""
    
    def __init__(self, nprocesses, bin_path, verbose=0):
        self.verbose = verbose
        self.nprocesses = nprocesses
        self.bin_path = bin_path
        self.workQueue = multiprocessing.Queue(0)
        self.returnQueue = multiprocessing.Queue(0)
        self.processes = []
        self.exitQueues = [multiprocessing.Queue(0) for processID in range(self.nprocesses)]
        self.closed = False
        for processID in range(self.nprocesses):
            process = SimulationTask(processID, self.bin_path, self.workQueue, self.returnQueue, self.exitQueues[processID], verbose=self.verbose)
            process.start()
            self.processes.append(process)
            
    def __enter__(self):
        return self
            
    def __exit__(self, type, value, traceback):
        """Removes _lock file if it was created and exits context."""
        self.close()
        return 0
            
    def close(self):
        self.closed = True
        exitcode = 'exit'
        for ii in range(self.nprocesses):
            self.exitQueues[ii].put(exitcode)
        for process in self.processes:
            process.join()
        time.sleep(0.01)


    def run(self, jobparam_list, timeout_per_job=None):
        njobs = len(jobparam_list)
        for ijob, jobparam in enumerate(jobparam_list):
            self.workQueue.put([ijob, jobparam, timeout_per_job])

        return_dict = {ii: None for ii in range(njobs)}
        jobs_left = set(range(njobs))
        try:
            while True:
                tuple_from_queue = get_from_queue(self.returnQueue)
                if tuple_from_queue is None:
                    for process in self.processes:
                        if not process.is_alive():
                            self.close()
                            raise Exception('A multiprocessing.Process died : this should not happen')
                    time.sleep(0.1)
                    continue
                    
                ijob, jobreturndict = tuple_from_queue
                return_dict[ijob] = jobreturndict
                jobs_left.remove(ijob)
                if len(jobs_left) == 0:
                    break
        except:
            print('Terminating processes, termination signal received...')
            self.close()
            raise
        return return_dict



class SimulationTask(multiprocessing.Process):
    
    def __init__(self, processID, bin_path, workQueue, returnQueue, exitQueue, verbose=0):
        if sys.version_info.major == 3:
            super().__init__()
        else:
            super(SimulationTask, self).__init__()
        self.verbose = verbose
        self.processID = processID
        self.bin_path = bin_path
        self.workQueue = workQueue
        self.returnQueue = returnQueue
        self.exitQueue = exitQueue
        self.simulation_process = None
        self.parent_exit_command = None

        
    def update_parent_exit_command(self):
        if self.parent_exit_command == 'exit':
            return
        while(True):
            command_loc = get_from_queue(self.exitQueue)
            if command_loc is None:
                break
            if command_loc == 'exit':
                self.parent_exit_command = 'exit'
                break
                

    def run(self):
        while True:
            #check if parent issued exit command
            self.update_parent_exit_command()
            if self.parent_exit_command is not None:
                break

            #get elements from workQueue
            tuple_seq = get_from_queue(self.workQueue, timeout=0.01)
            if tuple_seq is None:
                #no element to process
                continue
            ijob, jobparams, timeout = tuple_seq
                
            #the simulation process is not started, start it
            if self.verbose > 1:
                print('Starting job: %s'%(' '.join([self.bin_path] + jobparams)))
            returndict = dict()
            returndict['start_time'] = datetime.now()
            self.simulation_process = subprocess.Popen([self.bin_path] + jobparams, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            try:
                returndict['output'], returndict['error'] = self.simulation_process.communicate(timeout=timeout)
                returndict['status'] = 'completed'
            except:
                self.simulation_process.kill()
                returndict['output'], returndict['error'] = self.simulation_process.communicate()
                returndict['status'] = 'killed'
            for el in ['output', 'error']:
                returndict[el] = returndict[el].decode('utf8')
            returndict['stop_time'] = datetime.now()
            returndict['elapsed_time'] = (returndict['stop_time'] - returndict['start_time']).total_seconds()
            returndict['exitcode'] = self.simulation_process.returncode
            self.returnQueue.put([ijob, returndict])
        if self.verbose > 1:
            print('Simulation process exited...')
        
        

