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


import os, sys, shutil
import time
from datetime import datetime

class Timer(object):
    
    def __init__(self):
        self.start_date = datetime.utcnow()
        self.start_date_str = self.start_date.strftime('%Y-%m-%dT%H:%M:%S')
        self.times = [time.time()]
        self.last_time = self.times[0]
        
    def level_up(self):
        t_old, t_new = self.step()
        self.times.append(t_new)
        
    def level_down(self):
        if len(self.times) <= 1:
            raise Exception('Cannot go to 0 level')
        self.times = self.times[0:-1]
        
    def step(self):
        t_old, t_new = self.last_time, time.time()
        self.last_time = t_new
        return t_old, t_new
        
        
    def get_full_info(self):
        t_old, t_new = self.step()
        return {'last': t_new-t_old, 'level': t_new-self.times[-1], 'start': t_new-self.times[0], 'start_date': self.start_date_str}
        
    def __str__(self):
        t_old, t_new = self.step()
        return '%.3f seconds'%(t_new-t_old)

        
        
        
