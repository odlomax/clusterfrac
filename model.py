import numpy as np
import pandas as pd
from clusterfrac import params
from clusterfrac.cluster import star_cluster
from maths.fields.gaussian_random_field import scalar_grf
from maths.random.probability_density_function import pdf

class cluster_model:
    
    """
    
    Class to construct fractal Brownian motion star cluster model 
    
    """
    
    def __init__(self,model_type="2d2d"):
        
        """
        
        Subroutine: initialise star cluster
        
        Arguments
        ---------
        
        model_type: string
            values should either be 2d2d, 3d2d or 3d3d
            
            2d2d: Model is completely 2 dimensional
            
            2d3d: Model is 3d, but statistics are calculated in 2d
            
            3d3d: Model is completely 3 dimensional
        
        """
        
        
        if model_type=="2d2d":   
            self.model_params=params.model_2d2d_params
            
        elif model_type=="3d2d":
            self.model_params=params.model_3d2d_params
            
        elif model_type=="3d3d":
            self.model_params=params.model_3d3d_params
            
        else:
            print("Unknown model type:",model_type,"Using default 2d2d model.")
            self.model_params=params.model_2d2d_params
            
        return
    
    def fbm_field(self,H=None,sigma=None,L=None):
        
        """
        
        Function: produce a fractional brownian motion field
        
        Arguments
        ---------
        
        H: float
            Hurst index
        
        sigma:
            geometric standard deviation (base e)
            
        L:
            gaussian kernel size as fraction of box size
            
        Result
        ------
        
        field[:,:]: float
            Fractional Brownian motion field
        
        """
        
        # if parameters are not supplied, pick values from model defaults
        if H is None:
            H=np.random.uniform(*self.model_params.H_range)
        
        if sigma is None:
            sigma=np.random.uniform(*self.model_params.sigma_range)
            
        if L is None and self.model_params.gauss_filter:
            L=np.exp(np.random.uniform(*np.log(self.model_params.L_range)))
            
        print("H:    ",H)
        print("sigma:",sigma)
        if self.model_params.gauss_filter:
            print("L:    ",L)
            
        # set power spectrum exponent
        beta=self.model_params.n_dim_field+2.*H
        field=scalar_grf((self.model_params.grid_size,)*self.model_params.n_dim_field,beta,True)
           
        # convolve field with gaussian kernel
        if self.model_params.gauss_filter:
            field.gauss_conv(L*self.model_params.grid_size)
            
        # normalise field
        field.normalise(sigma,True)
        field.com_shift()
        
        real_field=field.signal.real
        
        # integrate out one of the dimensions if 2d3d
        if self.model_params.n_dim_field>self.model_params.n_dim_cluster:
            real_field=np.trapz(real_field,axis=0)
        
        return real_field
    
    def fbm_cluster(self,N,field=None):
            
        """
        
        Function: generate fbm cluster object with N points
        
        Arguments
        ---------
        
        N: int
            number points
            
        field[:,:]: real
            fbm field
            
        Result
        ------
        
        cluster: star_cluster object
        
        
        """
        
        # generate field, if not present
        if field is None:
            field=self.fbm_field()
        
        # generate pdf object
        field_pdf=pdf(field)
        
        # return cluster
        return star_cluster(field_pdf.random(N))
    
    def make_training_data(self,N_star_min,N_star_max,N_cluster):
        
        """
        
        Function: generate a pandas table of training data
        
        Arguments
        ---------
        
        N_star_min: int
            minimum number of stars in cluster
            
        N_star_max: int
            maximum number of stars in cluster
            
        N_cluster: int
            number of training clusters
            
        Result
        ------
        
        Table: pandas table object
            training data
            
        """
        
        # set number of stars
        N_star=np.random.uniform(N_star_min,N_star_max,N_cluster).astype(np.int)
        
        # set Hurst exponent
        H=np.random.uniform(*self.model_params.H_range,N_cluster)
        
        # set geometric standard deviation
        sigma=np.random.uniform(*self.model_params.sigma_range,N_cluster)
        
        # set kernel size
        if self.model_params.gauss_filter:
            L=np.exp(np.random.uniform(*np.log(self.model_params.L_range),N_cluster))
        else:
            L=np.zeros(N_cluster)

        # initialise observables
        R_hull=np.zeros(N_cluster)
        m_mean=np.zeros(N_cluster)
        m_m2=np.zeros(N_cluster)
        m_m3=np.zeros(N_cluster)
        m_m4=np.zeros(N_cluster)
        s_mean=np.zeros(N_cluster)
        s_m2=np.zeros(N_cluster)
        s_m3=np.zeros(N_cluster)
        s_m4=np.zeros(N_cluster)
        
        # generate clusters
        for i in range(N_cluster):
            
            print("Cluster",i+1,"of",N_cluster)
            
            # generate field
            field=self.fbm_field(H[i],sigma[i],L[i])
            
            # generate cluster
            cluster=self.fbm_cluster(N_star[i],field)
            
            # record observables
            R_hull[i]=cluster.hull_length()
            m_mean[i]=cluster.m_bar()
            m_m2[i]=cluster.m_moment(2)
            m_m3[i]=cluster.m_moment(3)
            m_m4[i]=cluster.m_moment(4)
            s_mean[i]=cluster.s_bar()
            s_m2[i]=cluster.s_moment(2)
            s_m3[i]=cluster.s_moment(3)
            s_m4[i]=cluster.s_moment(4)
            
        # make table of data
        data_array=np.array((N_star,H,sigma,L,R_hull,m_mean,m_m2,m_m3,m_m4,s_mean,s_m2,s_m3,s_m4)).T
        header=("N_star","H","sigma","L","R_hull","m_mean","m_m2","m_m3","m_m4","s_mean","s_m2","s_m3","s_m4")
        
        return pd.DataFrame(data=data_array,columns=header)