import numpy as np
import pandas as pd
from maths.points.whiten import whiten_points
from maths.graphs.mst import euclidean_mst

class star_cluster:
    
    """
    Class to generate basic star cluster
    
    """
    
    def __init__(self,r,whiten=True):
        
        """
        
        Subroutine: initialise star cluster
        
        Arguments
        ---------
        
        r[:,:]: float
            array of point positions [n_point:n_dim]
        
        """
        
        # scale/rotate coordinates so that cov(r)=I
        self.r=r
        if whiten:
            self.r=whiten_points(self.r)
        
        # generate minimum spanning tree 
        self.mst=euclidean_mst(self.r)
        
        # calculate minimum spanning tree normalisation
        self.mst_norm=(self.r.shape[0]-1)/(self.r.shape[0]**((self.r.shape[1]-1)/self.r.shape[1]))
        
    def s_bar(self):
        
        """
        
        Function: return mean edge-length of complete graph
        
        """
        
        return self.mst.com_edge_lengths().mean()
    
    def m_bar(self):
        
        """
        
        Function: return mean edge-length of minimum spanning tree
        
        """
        
        return self.mst.mst_edge_lengths().mean()*self.mst_norm

    def s_moment(self,n):

        """

        Function: return nth moment of complete graph edge-lengths

        """

        return self.mst.com_moment(n)
    
    def m_moment(self,n):

        """

        Function: return nth moment of complete graph edge-lengths

        """

        return self.mst.mst_moment(n)*self.mst_norm
    
    def hull_length(self):
        
        """
        
        Function: return n_dim root of convex hull area/volume
        
        """
        
        return self.mst.volume()**(1/self.r.shape[1])
    
    def max_radius(self):
        
        """
        
        Function: return greatest distance from centre
        
        """
        
        return np.sqrt(np.ax((self.r**2).sum(axis=1)))
    
    def make_table_row(self):
        
        """
        
        Function: return a single row pandas data table of cluster properties
        
        """
        
        data_array=np.array((self.r.shape[0],self.hull_length(),self.m_bar(),self.m_moment(2),self.m_moment(3),self.m_moment(4),self.s_bar(),self.s_moment(2),self.s_moment(3),self.s_moment(4))).reshape((1,10))
        header=("N_star","R_hull","m_mean","m_m2","m_m3","m_m4","s_mean","s_m2","s_m3","s_m4")
        
        return pd.DataFrame(data=data_array,columns=header)