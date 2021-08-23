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
import matplotlib.pyplot as plt
from hyfaa.common.yaml.yaml_parser import load_yaml, write_dict_to_yaml_file
import math 
import pandas as pd
#from SALib.sample import saltelli


def deg2rad(deg): 
    return deg * (math.pi/180)


def getDistanceFromLatLonInKm(x,y):
  R = 6371 # Radius of the earth in km
  dLat = deg2rad(y[1]-x[1])
  dLon = deg2rad(y[0]-x[0]) 
  a = math.sin(dLat/2)*math.sin(dLat/2)+math.cos(deg2rad(x[1]))*math.cos(deg2rad(y[1]))*math.sin(dLon/2)*math.sin(dLon/2) 
  c =2*math.atan2(math.sqrt(a),math.sqrt(1-a)) 
  d = R * c # Distance in km
  return d


def get_mesh_cell_centers_from_static_data_file(static_data_file):
    with netCDF4.Dataset(static_data_file) as ds:
        longitudes = ds.variables['longitude_center'][:]
        latitudes = ds.variables['latitude_center'][:]
    return longitudes, latitudes



def cross_diff(A, B):
    return A[:,None] - B[None,:]


def build_gaussian_error_fields(n_cells, nens, square, mgb_mesh_lon, mgb_mesh_lat, dlon=1., dlat=1.):
    
    
    lonmin, lonmax = np.floor(square['lonmin'])*1., np.ceil(square['lonmax'])*1.
    latmin, latmax = np.floor(square['latmin'])*1., np.ceil(square['latmax'])*1.
    lon = np.arange(lonmin, lonmax+dlon, dlon)
    lat = np.arange(latmin, latmax+dlat, dlat)
    nlon = len(lon)
    nlat = len(lat)

    #noise0 completely random
    # ~ noise0 = np.random.normal(0,1,(nlat,nlon,nens))
    
    #spatially correlated noise
    noise = np.zeros((nlat,nlon,nens), dtype=np.float64)
    
    theta = 3.
    #sigma = 0.5
    sigma = 1.
    mean = np.zeros(nlon*nlat, dtype=np.float64)
    grid_lon, grid_lat = np.meshgrid(lon,lat)
    grid_lon = grid_lon.flatten()
    grid_lat = grid_lat.flatten()
    #cov_square = theta**2*np.exp(-0.5*(cross_diff(grid_lon, grid_lon)**2 + cross_diff(grid_lat, grid_lat)**2) / sigma **2)
    cov_square = sigma**2*np.exp(-0.5*(cross_diff(grid_lon, grid_lon)**2 + cross_diff(grid_lat, grid_lat)**2) / theta **2)
    for iens in range(nens):
        gaussfield = np.random.multivariate_normal(mean, cov_square)
        noise[:,:,iens] = np.reshape(gaussfield, (nlat,nlon))

    
    noise3_mgb=np.zeros((nens,n_cells), dtype=np.float64)
   
    for i_cell in range(n_cells):
        ilon = int(np.floor((mgb_mesh_lon[i_cell]-lonmin)/dlon))
        ilat = int(np.floor((mgb_mesh_lat[i_cell]-latmin)/dlat))
        #noise0_mgb[i_cell,:] = noise0[ilat-1,ilon-1,:]
        noise3_mgb[:,i_cell] = noise[ilat,ilon,:]

    return noise3_mgb



def build_rand_fields(nc,ne):
    '''
    build matrix of random coefficients
    Inputs: 
    nc : number of cells
    ne : size of ensemble
    Outputs:
    Error matrix of dim (ncells,nens)
    '''
    W=[]
    for ni in range(0,ne):    
        W.append(np.random.random((nc)))
    return W


def build_normal_spatial_noise(nc,ne,mean=0,dev=1):
    '''
    build matrix of coefficients following a normal distribution
    Inputs: 
    nc : number of cells
    ne : size of ensemble
    mean : center of normal distribution
    dev : standard deviation of normal distribution
    Outputs:
    Error matrix of dim (ncells,nens)
    '''
    if ne == 1:
        return np.zeros((ne,nc))
    elif ne > 1:
        return np.random.normal(np.ones((ne,nc))*mean,dev)
    raise Exception('ne < 1: %s'%ne)


def build_njissen_fields(nc,ne,mean=0,dev=1,bias=0.2, error=0.5):
    '''
    build matrix of coefficients following a njissen distribution
    Inputs: 
    nc : number of cells
    ne : size of ensemble
    mean : center of normal distribution
    dev : standard deviation of normal distribution
    bias : bias of distribution
    error : distribution error
    Outputs:
    Error matrix of dim (ncells,nens)
    '''
    W=[]
    for ni in range(0,ne):    
        dist=np.random.normal(0,1,size=nc)
        err=error**2+1
        W.append([])
        [W[ni].append(((1+bias)/np.sqrt(err))*exp(np.sqrt(np.log(err))*dist[i])) for i in range(nc)] 
    return W



def build_lognormal_fields(nc,ne,mean=0,dev=1):
    '''
    build matrix of coefficients following a lognormal distribution
    Inputs: 
    nc : number of cells
    ne : size of ensemble
    mean : center of normal distribution
    dev : standard deviation of normal distribution
    '''
    W=[]
    for ni in range(0,ne):    
        W.append(np.random.lognormal(mean,dev,size=nc))
    return W



def plot_rand_ensemble(A,mean=0.,dev=1.):
    
    fig = plt.figure()
    n=len(A)
    npoints=len(A[0])
    ny=int(sqrt(n))
    nx=ny+(n-ny**2)

    
    #for i in range (0,len(A)):
        #ax = fig.add_subplot(n,  1,i+1)
        
        #ax.imshow(A[i])
    sigma = sqrt(dev)    
    x=np.linspace(mean-3*sigma,sigma+3*dev,npoints)
    for i in range(1,ny*nx+1):            
            if i >n:
                continue
            ax = fig.add_subplot(nx,ny,i)
            ax.plot(x,A[i-1],'.')   
            
    return plt.show()


def get_obs_vector(assim_pandas_var):
        
    dim_x=len(assim_pandas_var['value'])
        
    obs_vector=np.full((dim_x,2),1e+20)
        
    obs_vector[:,0]=assim_pandas_var['value']
    obs_vector[:,1]=assim_pandas_var['mesh_id']
        
    return obs_vector



def get_perturbed_multfact_vectors_from_hydrostate(hydrostate_ordered_paths):
    n_ensemble = len(hydrostate_ordered_paths)
    with netCDF4.Dataset(hydrostate_ordered_paths[0]) as ds:
        n_meshes = ds.dimensions['n_meshes'].size
        perturbed_multfact_vectors = np.empty((n_ensemble, n_meshes), dtype=ds.variables['last_rain_data_loaded'].dtype)
        perturbed_multfact_vectors[0,:] = ds.variables['last_rain_data_loaded'][:]
    for i_ensemble in range(1, n_ensemble):
        with netCDF4.Dataset(hydrostate_ordered_paths[i_ensemble]) as ds:
            perturbed_mutlfact_vectors[i_ensemble,:] = ds.variables['last_rain_data_loaded'][:]
    return perturbed_multfact_vectors


def rain_perturbation(perturbed_previous_multfact_vectors, actual_rain_vector,error):
    """returns an array of "perturbed" rain vectors"""   
    #factnorm=np.random.normal(0,1,20)
    #normmat=np.full((20,11595),1.0)
    #for it in range(0,20):
    #    normmat[it,:]=1+factnorm[it]
    #fact = (1./np.sqrt(error**2+1))*np.exp(np.log(error**2+1)*perturbed_previous_multfact_vectors)
    fact = (1./np.sqrt(error**2+1))*np.exp(np.sqrt(np.log(error**2+1))*perturbed_previous_multfact_vectors)
    perturbed_rain = np.multiply(fact, actual_rain_vector)
    
    return perturbed_rain


def dummy_initial_rain_perturbation(actual_rain_vector, rain_vector_uncertainty, n_ensemble):
    if n_ensemble == 1:
        return np.zeros((n_ensemble, len(actual_rain_vector)), dtype=actual_rain_vector.dtype)
    return np.tile(actual_rain_vector, (n_ensemble,1)) + \
        np.multiply(np.tile(rain_vector_uncertainty, (n_ensemble,1)), \
        np.random.randn(n_ensemble, len(actual_rain_vector)))
        

def get_multfactor_ensemble(old_multfact_vector):
    return old_multfact_vector



def write_eigenvalues (eig,nobs,filein='eigenvalues.dat'):
    
    exists = os.path.isfile(filein)
    if exists:
        with  open(filein, 'a') as wfile:
            wfile.write('ZONE  F=POINT, I=%i, J=1 K=1'%nobs)
            for i in range(1,nobs):
                print(len(eig),nobs-i)
                wfile.write(i,eig[nobs-i])
    else:
    # Keep presets
        with  open(filein, 'a') as wfile:
            wfile.write('TITLE = "Eigenvalues of C"')
            wfile.write('VARIABLES = "obs" "eigenvalues"')
            wfile.write('ZONE  F=POINT, I=%i, J=1 K=1'%nobs)
            for i in range(1,nobs):
                wfile.write(i,eig[nobs-i])
                
    return
        
def get_significant_eigenvalues(eig, nobs, truncation=0.999, verbose=0):
    ########
    ##Todo : coder les autres cas
    ##########
    sigsum=np.sum(eig)
    sigsum1=0.0
    nrsigma=0
    
    for i in range(nobs,0,-1):
        if (sigsum1/sigsum < truncation) :
            nrsigma=nrsigma+1
            sigsum1=sigsum1+eig[i-1]
            eig[i-1] = 1.0/eig[i-1]
        else :
            eig[0:i]=0.0
            break
      
    if verbose > 0:
        print('analysis: Number of dominant eigenvalues: %i of %i'%(nrsigma,nobs))
        print('analysis: Share (and truncation): %f (%f)'%(sigsum1/sigsum,truncation))
          
          
    return eig
    
    
def genX3(nens,nobs,nrmin,eig,W,D):
        
    X1=np.empty((nrmin,nobs))
    for i in range(nrmin):
        for j in range(nobs):
            X1[i,j]=eig[i]*W[j,i]
                
    X2=matmul_bigdata(X1,D)
    X3=matmul_bigdata(W,X2)
        
    return X3
    
def multa_new(A,X,ndim,nens,iblkmax):
        
    v=np.empty((iblkmax,nens))
        
    for ia in range(0,ndim,iblkmax):
        ib=min(ia+iblkmax-1,ndim)
        v[1:ib-ia+1,1:nens] = A[ia:ib,1:nens]
        A[ia,1]=v[1,1]*X[1,1]
        
    return A
            
            
            
            
            
################################
#perturb static data parameters (f(cells) only)
            
def generate_perturbed_static_data_files(static_data_file_in, folder_out, dico_params, n_ensemble, perturb_type, perturb_mode):
    if not os.path.exists(folder_out):
        os.system('mkdir -p %s'%folder_out)
    if len(os.listdir(folder_out)) > 0:
        raise Exception('folder %s not empty'%folder_out)
        
    write_dict_to_yaml_file(dico_params, '%s/varying_parameters.yaml'%folder_out)
    
    with netCDF4.Dataset(static_data_file_in) as ds:
        ncells = ds.dimensions['n_cells'].size
    coef_matrix = perturb_parameters_static_data(ncells, n_ensemble, dico_params, perturb_type, perturb_mode)
    for i_ensemble in range(n_ensemble):
        static_data_file = '%s/static_data_%d.nc'%(folder_out, i_ensemble)
        shutil.copy(static_data_file_in, static_data_file)
        with netCDF4.Dataset(static_data_file, mode='a') as ds:
            for ii, dico_loc in enumerate(dico_params):
                values = ds.variables[dico_loc['name']][:]*coef_matrix[i_ensemble,ii,:]
                
                #garde fou valeurs
                if 'min' in dico_loc:
                    if dico_loc['min'] is not None:
                        values[values<dico_loc['min']] = dico_loc['min']
                if 'max' in dico_loc:
                    if dico_loc['max'] is not None:
                        values[values>dico_loc['max']] = dico_loc['max']
                        
                ds.variables[dico_loc['name']][:] = values



    
def perturb_parameters_static_data(ncells, nens, dico_parameter, perturb_type, perturb_mode):
    
    nparams=len(dico_parameter)
    error_list = [dico_loc['error'] for dico_loc in dico_parameter]
        
    if perturb_type=='normal':
        if perturb_mode=='per_variable':
            coef_matrix = build_normal_sample_per_variable(ncells, nparams, nens, error_list)
        elif perturb_mode=='per_cell':
            coef_matrix = build_normal_sample_per_cell(ncells, nparams, nens, error_list)
    elif perturb_type=='saltelli':
        if perturb_mode=='per_variable':
            coef_matrix = build_saltelli_sample_per_variable(ncells, nparams, nens, error_list)
        elif perturb_mode=='per_cell':
            raise Exception('saltelli can not be done per cell, use per variable mode')
        else:
            raise Exception('mode %s unknown'%perturb_mode)
    else:
        raise Exception('type %s unknown'%perturb_type)
        
    for ii, dico_loc in enumerate(dico_parameter):
        #garde fou valeurs coeff
        if 'coeff_min' in dico_loc:
            if dico_loc['coeff_min'] is not None:
                coef_matrix[:,ii,:][coef_matrix[:,ii,:]<dico_loc['coeff_min']] = dico_loc['coeff_min']
        if 'coeff_max' in dico_loc:
            if dico_loc['coeff_max'] is not None:
                coef_matrix[:,ii,:][coef_matrix[:,ii,:]>dico_loc['coeff_max']] = dico_loc['coeff_max']
            
    return coef_matrix
    
    
    
def gaussian_bijection_to_exp(nitems, multiplier_divider_max, sigma_max=2.):
    """ generate random values between [1/multiplier_divider_max, multiplier_divider_max], centered on 1, whose logarithm is a gaussian distribution
    """
    #generate random values within a gaussian. Hard 2-sigma limit to prevent extreme values
    mult_values = np.random.randn(nitems)
    mult_values[mult_values<-sigma_max] = -sigma_max
    mult_values[mult_values>sigma_max] = sigma_max
    #compute values for each parameter : the gaussian is dilated to a [1/a, a] interval centered on 1 using exp(mult_values)
    return np.exp(np.log(np.sqrt(multiplier_divider_max))*mult_values)

def build_saltelli_sample_per_variable(ncells,nparams,nens,error_list):
    final_sample=np.zeros((nens,nparams,ncells))
    ninput=math.ceil(nens/(2*nparams+2))
    bounds_matrix=[]
    for i,err in enumerate(error_list):
        bounds_matrix.append([-err,err])
    problem={'num_vars':nparams,
             'bounds':bounds_matrix}
    sample=saltelli.sample(problem, ninput)
    nsize=np.shape(sample)[0]
    for iens in range(nens):
        for iparam in range(nparams):
            final_sample[iens,iparam,:]=[1+sample[iens,iparam]]*ncells
    return final_sample


def build_normal_sample_per_variable(ncells,nparams,nens,error_list):
    
    sample=np.zeros((nens,nparams))
    final_sample=np.zeros((nens,nparams,ncells))
    for param in range(nparams):
        sample[:,param]=np.random.normal(1.0,error_list[param],nens)
        for iens in range(nens):
            final_sample[iens,param,:]=[sample[iens,param]]*ncells        
    return final_sample

def build_normal_sample_per_cell(ncells,nparams,nens,error_list):
    
    sample=np.zeros((nens,nparams,ncells))
    for param in range(nparams):
        for ens in range(nens):
            sample[ens,param,:]=np.random.normal(1.0,error_list[param],ncells)
    return sample


def apply_velocity_loc(mesh_size,nens,nsize,v_file,obs_vec,lon,lat,timestep):
    
    
    loc_vec=np.zeros((mesh_size))#lui donner la taille de ncells puis etendre
    loc_temp=np.zeros((mesh_size))
    
    Vmean_pd=pd.read_csv(v_file) 
    Vmean=np.array(Vmean_pd.V)
    #print(obs_vec) 
    for elem in obs_vec:
        dcrit=Vmean[int(elem)-1]*timestep*24*3600./1000.
        #print(dcrit)
        x=[lon[int(elem)-1],lat[int(elem)-1]]
        for i in range(mesh_size):
            dist=getDistanceFromLatLonInKm(x,[lon[i],lat[i]])
            
            #loc_temp=np.where(getDistanceFromLatLonInKm(x,[lon,lat])>dcrit,1.0,0.0)
            if dist>dcrit:
                loc_temp[i]=0
            else:
                #print(dist)
                loc_temp[i]=1
            
        loc_vec+=loc_temp
    
    loc_vec=np.where(loc_vec>0,1.0,0.0)
    loc_arr=np.tile(loc_vec,(nens,1))    
    num=nsize/mesh_size
    loc_arr2=np.tile(np.transpose(loc_arr),(int(num),1))
    return loc_arr2
            
            
            
    
