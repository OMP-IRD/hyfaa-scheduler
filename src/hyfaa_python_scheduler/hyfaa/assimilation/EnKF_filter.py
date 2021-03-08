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


from hyfaa.common.common_functions import *
import pandas as pd
import random

import hyfaa.common.yaml.yaml_parser as yaml_parser
from hyfaa.common.test.common_test_functions import make_hydrostate, success_message_utest, fail_message_utest, print_utest_message
import hyfaa.assimilation.assim_tools as assim_tools



class EnKFobj(object):
    '''EnKF entity class'''
    

    def __init__(self, assim_data, n_ensemble, nC, nU, lon, lat, model_timestep):
        '''
        This function initializes EnKF variables as defined by Evensen 2003
        Inputs :
        assim_data : DA parameters file
        n_ensemble : size of control ensemble
        nC : number of unit catchments (cells)
        nU : number of hydrological response units
        
        Outputs :
        self.Nens_i : size control ensemble
        self.list_contr_vars_s : list of control variables
        self.Nstatvars_i : number of control variables
        self.Ak_old : control ensemble
        self.Ak_new : analysis ensemble
        self.HA : matrix HA (product of observation operator H and ensemble control matrix A)
        self.Dprim : innovation matrix (difference of observation ensemble and HA)
        self.Sk : Sk=HA'=HA-median(HA)
        self.Rk : observation error covariance matrix
        '''

        if hasattr(assim_data, 'keys'):
            dico = assim_data
        elif os.path.exists(assim_data):
            dico = yaml_parser.load_yaml(assim_data, env_vars=True)
        else:
            raise Exception('There is no file for DA parameters')
        self.Nens_i = n_ensemble
        self.mesh_lon_r=lon
        self.mesh_lat_r=lat
        self.list_contr_vars_s = []
        self.dict_sizes = dict()
        self.mesh_size=nC
        
        temp = 0
        #statvars assimilation (hydrological state variables)
        for field, activated in dico['statvars'].items():
            if not activated:
                continue
            self.list_contr_vars_s.append(field) 
            if field in ['soil_moisture','water_canopy_interception_volume']:
                self.dict_sizes[field] = {'size': nC*nU}
                self.dict_sizes[field].update({'hru': nU})
                self.dict_sizes[field].update({'hru_unit_size': nC})
            else :
                self.dict_sizes[field] = {'size': nC}
                self.dict_sizes[field].update({'hru': 1})
            self.dict_sizes[field].update({'start_index': temp, 'end_index': temp+self.dict_sizes[field]['size']})
            
            temp += self.dict_sizes[field]['size']

        #static parameter assimilation
        if dico['parameters']['activate']:
            for field,activated in dico['parameters']['param_list'].items():
                if not activated:
                    continue
                self.list_contr_vars_s.append(field)
                self.dict_sizes[field] = {'size': nC}
                self.dict_sizes[field].update({'hru': 1})
                self.dict_sizes[field].update({'start_index': temp, 'end_index': temp+self.dict_sizes[field]['size']})
                temp += self.dict_sizes[field]['size']

        self.Nstatvars_i = temp
        
        self.Ak_old = None
        self.Ak_new = None
        self.HA = None
        self.Dprim = None
        self.Sk = None
        self.Rk = None

        self.loc_dico=dico['localization']
        if self.loc_dico['activate']:
            self.loc = np.loadtxt(self.loc_dico['loc_file'])
        self.model_timestep=model_timestep
        
    def clear_object_for_new_ensemble_size(self, n_ensemble):
        
        self.Nens_i = n_ensemble
        self.Ak_old = None
        self.Ak_new = None
        self.HA = None
        self.Dprim = None
        self.Sk = None
        self.Rk = None
  
    def ensemble_files_to_matrix(self, ordered_filenames):
        '''
        Build ensemble matrix from files
        Inputs : 
        oredered_filenames : list of files
        Outputs : 
        matrix_out (Nstatvars_i,nens_i): matrix containing the sensemble of control vectors Nstatvars_i: number of control variables nens_i : size of ensemble 
        '''
        if len(ordered_filenames) != self.Nens_i:
            raise Exception('length of ordered_filenames (%d) different from assumed ensemble size Nens_i (%d)'%(len(ordered_filenames), self.Nens_i))

        matrix_out = np.empty((self.Nstatvars_i, self.Nens_i), dtype=np.float64)
        for ii, filename in enumerate(ordered_filenames):
            
            with netCDF4.Dataset(filename) as dataset:
                missing_vars = set(self.list_contr_vars_s) - set(dataset.variables.keys())
                
                if len(missing_vars) > 0:
                    raise Exception('missing variables from file %s:\n%s\n'%(filename, '\n'.join([' - %s'%el for el in sorted(missing_vars)])))
                
                for field in self.list_contr_vars_s:
                   
                    if self.dict_sizes[field]['hru']>1:
                        temp_index=0
                        for hh in range(self.dict_sizes[field]['hru']):
                            matrix_out[self.dict_sizes[field]['start_index']+temp_index:self.dict_sizes[field]['start_index']+temp_index+np.shape(dataset.variables[field][:])[1],ii] = dataset.variables[field][hh,:]
                            temp_index+=np.shape(dataset.variables[field][:])[1]
                    else : 
                        matrix_out[self.dict_sizes[field]['start_index']:self.dict_sizes[field]['end_index'],ii] = dataset.variables[field][:]

        return matrix_out
        

    def build_control_ensemble(self, ordered_filenames):
        '''
        Build ensemble matrix from files
        inputs:
        ordered_filenames : list of files
        Outputs : 
        '''
        if self.Ak_old is not None:
            if np.shape(self.Ak_old) != (self.Nstatvars_i, self.Nens_i):
                raise Exception('np.shape(self.Ak_old) != (self.Nstatvars_i, self.Nens_i)')
            self.Ak_old[:] = self.ensemble_files_to_matrix(ordered_filenames)
        else:
            self.Ak_old = self.ensemble_files_to_matrix(ordered_filenames)
        
        
    
    def get_obs_ensemble(self, vector, ObsError=0.2):
        '''
        Build observation ensemble
        Inputs:
        vector : observation vector
        ObsError : standard deviation of observations
        Outputs:
        ensemble : observation ensemble matrix
        '''
        random.seed(10)
        noise_l = assim_tools.build_normal_spatial_noise(len(vector), self.Nens_i, mean=0., dev=1.)
        #noise_m = np.array([np.array(xi) for xi in noise_l])
        noise_m = vector * noise_l * ObsError
        #ensemble=np.empty((len(vector),self.Nens_i))
        #for i in range(self.Nens_i):
        ensemble=noise_m + vector 
        ensemble = np.where(ensemble<0, 0, ensemble)
        ensemble=np.transpose(ensemble)
        return ensemble    
        
   
    def get_HA(self, list_id):
       '''
       Calulates HA (projection of control vector into observation space)
       Inputs :
       list_id : list of control variables starting index
       output :
       HA : (product of observation operator H and ensemble control matrix A)
       '''
       if self.Ak_old is None:
           raise Exception('get_HA requires self.Ak_old to be initialized')
       self.HA = np.empty((len(list_id),self.Nens_i))
       index = self.list_contr_vars_s.index('streamflow_catchment')
       
       start = sum([self.dict_sizes[field]['size'] for field in self.list_contr_vars_s[0:index]])
       
       for ii,mesh_id in enumerate(list_id):
           for iens in range(self.Nens_i):
                self.HA[ii,iens] = self.Ak_old[int(start-1+mesh_id),iens]   

        
    def get_innovation_matrix(self, obs_ensemble):
        '''
        computes innovation matrix D'=Y-HA##
        inputs:
        obs_ensemble : observation ensemble matrix 
        outputs:
        Dprim : innovation matrix (difference of observation ensemble and HA)
        '''
        if self.HA is None:
           raise Exception('get_innovation_matrix requires self.HA to be initialized')
        self.Dprim = obs_ensemble - self.HA

        
        
    def get_matrix_Sk(self):
        '''
        computes Sk=HA'=HA-median(HA)
        Inputs:
        Outputs : Sk
        '''
        if self.HA is None:
           raise Exception('get_matrix_Sk requires self.HA to be initialized')
        self.Sk = np.empty((np.shape(self.HA)[0],self.Nens_i))
        for iens in range(self.Nens_i):
            for ii in range(np.shape(self.HA)[0]):
                
                self.Sk[ii,iens]=self.HA[ii,iens]-np.mean(self.HA,axis=1)[ii]
        
        

        
                
    def get_innovation_median(self, obs_ensemble):
        '''
        calculates innovation mean
        inputs :
        obs_ensemble : observation ensemble
        outputs :
        innov_median : innovation mean
        '''
        if self.HA is None:
           raise Exception('get_matrix_Sk requires self.HA to be initialized')
        self.innov_median = np.empty((np.shape(self.HA)[0],self.Nens_i))
        for iens in range(self.Nens_i):
            for ii in range(np.shape(self.HA)[0]):
                
                self.innov_median[ii,iens]=obs_ensemble[ii,iens]-np.mean(self.HA,axis=1)[ii]    
        

    def get_observation_covariance_matrix(self, obs_vector, ErrObs=0.2):
        '''computes observation error covariance matrix
           Inputs:
           obs_vector : observation vector
           ErrObs : observation error standard deviation
           Outputs:
           Rk : observation error covariance matrix
        '''
        self.Rk = np.zeros((len(obs_vector),len(obs_vector)))
        for ii in range(len(obs_vector)):
            self.Rk[ii,ii] = (obs_vector[ii]*ErrObs)**2 

    def get_HA_new(self, list_id):
       '''
       calculates projection of analysis into observation space
       Inputs : 
       list_id : list of variables starting index
       Outputs:
       HA_new : projection of analysis into observation space
       '''
       if self.Ak_new is None:
           raise Exception('get_HA requires self.Ak_new to be computed')
       self.HA_new = np.empty((len(list_id),self.Nens_i))
       index = self.list_contr_vars_s.index('streamflow_catchment')
       
       start = sum([self.dict_sizes[field]['size'] for field in self.list_contr_vars_s[0:index]])
       
       for ii,mesh_id in enumerate(list_id):
           for iens in range(self.Nens_i):
                self.HA_new[ii,iens] = self.Ak_new[int(start-1+mesh_id),iens]
                
    def get_analysis(self, obs_vec,obs_val):
        '''
        Calculates analysis using EnKF as described in Evensen 2003
        Inputs:
        Outputs:
        Ak_new : analysis ensemble matrix
        '''

        #Computes analysis
        if (self.HA is None) or (self.Dprim is None) or (self.Sk is None) or (self.Rk is None):
           raise Exception('get_matrix_Sk requires self.HA, self.Dprim,self.Sk and self.Rk  to be initialized')

        nens = self.Nens_i
        ndim = self.Nstatvars_i
        nobs = np.shape(self.Rk)[0]
        
        if self.Ak_new is None:
            self.Ak_new = np.empty((self.Nstatvars_i, self.Nens_i), dtype=np.float64)
        else:
            if np.shape(self.Ak_new) != (self.Nstatvars_i, self.Nens_i):
                raise Exception('np.shape(self.Ak_new) != (self.Nstatvars_i, self.Nens_i)')
        
        #Computes A'
        I1N = np.empty((nens,nens))
        for j in range (nens):
            for k in range (nens):
                if j == k:
                    I1N[j,k] = 1.-1./nens
                else :
                    I1N[j,k] = -1./nens


        Aprim = matmul_bigdata(self.Ak_old, I1N)
        
        if self.loc_dico['activate']:
        #matrix of localization
            rhoc=np.zeros((nobs,nobs))
            rho=np.zeros((ndim,nobs))
            rho_temp = np.zeros((ndim, nobs))
            rhoc_temp = np.zeros((nobs, nobs))
        
        
            for kk in range(nobs):
                for field in self.list_contr_vars_s:

                    if self.dict_sizes[field]['hru'] > 1:
                        temp_index = 0
                        for hh in range(self.dict_sizes[field]['hru']):
                        
                            rho[int(self.dict_sizes[field]['start_index'] + temp_index ): int(self.dict_sizes[field]['start_index'] + temp_index + self.dict_sizes[field]['hru_unit_size']) , kk ] = self.loc[int(obs_vec[kk]-1),:]
                            temp_index += self.dict_sizes[field]['hru_unit_size']
                    else:
                    
                        rho[int(self.dict_sizes[field]['start_index']):int(self.dict_sizes[field]['end_index']),kk] = self.loc[int(obs_vec[kk]-1),:]
            for kk,val in enumerate(obs_val):
                if val>=0:
                
                    rhoc[kk, :] = rho[int(obs_vec[kk] - 1), :]
        
        ##pseudo inversion of C=SS'+(N-1)R
        if nobs==1:
            nrmin=1
            W=np.full((1,1),1.0)
            eig=matmul_bigdata(self.Sk,np.transpose(self.Sk))+float(nens-1)*self.Rk
            eig=1./eig
            if self.loc_dico['activate']:
                Sst=np.full((1,1),1.0)
        
        else :
            ##R=S*S'+(nens-1)*R
            nrmin=nobs
            if self.loc_dico['activate']:
                SSt=matmul_bigdata(self.Sk,np.transpose(self.Sk))
                Cloc = np.multiply(rhoc, SSt) + (nens - 1) * self.Rk
                c1temp=np.multiply(rhoc, SSt)
                c2temp=multiply_bigdata(rhoc, SSt)
            
                ##compute eigen decomposition of R --> W*eig*W'
                (eig,W)=np.linalg.eig(Cloc) #loc
            else:
                nrmin=nobs
                Rtemp=matmul_bigdata(self.Sk,np.transpose(self.Sk))+float(nens-1)*self.Rk
                ##compute eigen decomposition of R --> W*eig*W'
                (eig,W)=np.linalg.eig(Rtemp)
        
        ##significant eigen values
        
            eig=assim_tools.get_significant_eigenvalues(eig,nobs)

        ##calculation of X3,X5
        ########
        ##Todo : coder les autres cas
        ##########
        if nobs>1:
            X3=assim_tools.genX3(nens,nobs,nrmin,eig,W,self.Dprim)
        else:
            X3=self.Dprim*eig[0]
    
        if self.loc_dico['activate']: 
            lreps=True
            Reps=matmul_bigdata(Aprim,np.transpose(self.Sk))
            Reps=np.multiply(Reps,rho)
        elif self.loc_dico['activate']==False and (2*ndim*nobs < 1*nens*(nobs+ndim)):
            lreps=True
            Reps=matmul_bigdata(Aprim,np.transpose(self.Sk))
        else:
            lreps=False
            X5=matmul_bigdata(np.transpose(self.Sk),X3)
            XX=np.einsum('ii->i',X5)
            XX+=1

        #Final ensemble update
        
        if lreps :
            self.Ak_new = self.Ak_old+matmul_bigdata(Reps,X3)
        else :
            iblkmax = min(ndim,200)
            self.Ak_new = assim_tools.multa_new(self.Ak_old,X5,ndim,nens,iblkmax)

        #constraint increment (Do not remove this line, it can be important to activate it depending on the study ... still in research)
        #self.Ak_new=np.where((abs(self.Ak_new-self.Ak_old)>2*self.Ak_old) & (abs(self.Ak_new-self.Ak_old)>0.001),self.Ak_old+2*self.Ak_old*abs(self.Ak_new-self.Ak_old)/(self.Ak_new-self.Ak_old),self.Ak_new)
        

    def perform_analysis(self, obs_db, input_files, inplace=True, verbose=0, timestep_analysis=None):

        #get observations
        start_time = time.time()
        obs_db_vector = assim_tools.get_obs_vector(obs_db)
        if verbose > 1:
            print('Observation vector %s'%(time.time()-start_time))
            # ~ print(obs_db_vector[:,0])
            start_time = time.time()
            
        self.build_control_ensemble(input_files)
        
        #this is the 2nd most time consuming task after writes at the end => its basically a read, must try with compression off to see if its faster, and try to read files in parallel
        if verbose > 1:
            print('Control variable matrix %s'%(time.time()-start_time))
            # ~ print(self.Ak_old)
            start_time = time.time()
            
        #Build observation ensemble D
        obs_ensemble = self.get_obs_ensemble(obs_db_vector[:,0])
        if verbose > 1:
            print('Observation ensemble %s'%(time.time()-start_time))
            # ~ print(obs_ensemble)
            start_time = time.time()
            
        self.get_HA(obs_db_vector[:,1])
        if verbose > 1:
            print('HA %s'%(time.time()-start_time))
            # ~ print(self.HA)
            start_time = time.time()
            
        self.get_innovation_matrix(obs_ensemble)
        if verbose > 1:
            print('innovation %s'%(time.time()-start_time))
            # ~ print(self.Dprim)
            start_time = time.time()
            
        self.get_matrix_Sk()
        if verbose > 1:
            print('S %s'%(time.time()-start_time))
            # ~ print(self.Sk)
            start_time = time.time()
            
        self.get_innovation_median(obs_ensemble)
        if verbose > 1:
            print('Innov median %s'%(time.time()-start_time))
            # ~ print(self.innov_median)
            start_time = time.time()
            
        self.get_observation_covariance_matrix(obs_db_vector[:,0])
        if verbose > 1:
            print('R %s'%(time.time()-start_time))
            # ~ print(self.Rk)
            start_time = time.time()
        
        
        self.get_analysis(obs_db_vector[:,1],obs_db_vector[:,0])
        if verbose > 1:
            print('Ak_new %s'%(time.time()-start_time))
            # ~ print(self.Ak_new)
            start_time = time.time()
        
        self.get_HA_new(obs_db_vector[:,1])
        if verbose > 1:
            print('HA_new %s'%(time.time()-start_time))
            # ~ print(self.HA_new)
            start_time = time.time()
            
        #this is the most time consuming task => easily optimized with parallel computing, and try with compression off
        corrected_files = []
        for i_ensemble, file_in in enumerate(input_files):
            if inplace:
                file_out = file_in
            else:
                file_out = file_in.replace('.nc', '_assim.nc')
                shutil.copy(file_in, file_out)
            with netCDF4.Dataset(file_out, mode='a') as ds:
                for field in self.list_contr_vars_s:
                    if self.dict_sizes[field]['hru'] > 1:
                        temp_index = 0
                        for hh in range(self.dict_sizes[field]['hru']):
                            ds.variables[field][hh,:] = self.Ak_new[self.dict_sizes[field]['start_index']+temp_index:self.dict_sizes[field]['start_index']+temp_index+self.dict_sizes[field]['hru_unit_size'],i_ensemble]
                            temp_index += self.dict_sizes[field]['hru_unit_size']
                    else:
                        ds.variables[field][:] = self.Ak_new[self.dict_sizes[field]['start_index']:self.dict_sizes[field]['end_index'],i_ensemble]

            if timestep_analysis is not None:
                with netCDF4.Dataset(file_out) as ds:
                    flood_timestep_conditions = timestep_analysis['alpha_inertial'] * timestep_analysis['river_length'] * 1000 / (9.81*ds.variables['water_depth_catchment'][:])
                    min_id = np.argmin(flood_timestep_conditions)
                    print('      Ensemble member %d: timestep=%s at mesh id %d'%(i_ensemble, flood_timestep_conditions[min_id], min_id))
            corrected_files.append(file_out)
        if verbose > 1:
            print('write assimilation results %s'%(time.time()-start_time))
        return corrected_files




def test_perform_analysis():
    
    test_dir = tempfile.mkdtemp(dir=os.path.abspath('.'), prefix='testperfanalysis_')
    from check_parameters_main import get_scheduler_basic_main_parameters_dict, get_scheduler_standard_assimilation_dict
    dico_main = get_scheduler_basic_main_parameters_dict()
    dico_assim = dico_main['assim_params_file']
    n_ensemble, n_meshes, n_hrus = 6, 11595, 10
    
    print_utest_message('Test EnKF method: success if runs without errors')
    try:
        print_utest_message('Starting EnKF_obj')
        EnKF_obj = EnKFobj(dico_assim, n_ensemble, n_meshes, n_hrus, None, None, dico_main['model_min_time_step'])
        input_files = [os.path.join(test_dir, 'hydrostate_%d.nc'%ii) for ii in range(n_ensemble)]
        print_utest_message('Making hydrostate files')
        for el in input_files:
            make_hydrostate(el, static_assimilation_params=list(dico_assim['parameters']['param_list'].keys()))
        data = [[1,1000.], [3969,1000.], [11595,1000.]]
        print_utest_message('Making assimilation data')
        obs_pd = pd.DataFrame(data, columns = ['mesh_id', 'value'])
        print_utest_message('Performing analysis')
        _ = EnKF_obj.perform_analysis(obs_pd, input_files, inplace=False, verbose=2)
        success_message_utest()
    except:
        fail_message_utest()
    shutil.rmtree(test_dir)
    
    
if __name__ == '__main__':
    test_perform_analysis()
    

