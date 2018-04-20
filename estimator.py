import numpy as np
from sklearn.neural_network import MLPRegressor
from clusterfrac.params import hyper_params

def get_X(data):
    
    """
    
    Function: convert pandas data table to sklearn X variable
    
    Arguments
    ---------
    
    data: panadas data table
    
    Result
    ------
    
    X[:,:]: float
        sklearn X variable
    
    """
    
    return np.array((data["R_hull"],data["m_mean"],data["m_m2"],data["m_m3"],data["m_m4"],data["s_mean"],data["s_m2"],data["s_m3"],data["s_m4"])).T
    

def get_Y(data):
    
    """
    
    Function: convert pandas data table to sklearn Y variable
    
    Arguments
    ---------
    
    data: panadas data table
    
    Result
    ------
    
    Y[:,:]: float
        sklearn Y variable
    
    """
    
    return np.array((data["H"],data["sigma"])).T
    
    
class param_estimator:
    
    """
    
    Class to construct model parameter estimator
    
    """
    
    
    def __init__(self,data_table):
        
        """
        
        Function: use data set to train model
        
        Arguments
        ---------
        
        data_table: pandas data frame
            set of training statistics from ensemble of artificial clusters
        
        """
        
        # get index to split tables
        i_split=int(hyper_params.data_split*data_table.shape[0])
        
        self.training_data=data_table.loc[:i_split-1,:]
        self.test_data=data_table.loc[i_split:,:]
        
        
        X_train=get_X(self.training_data)
        
                         
        Y_train=get_Y(self.training_data)
        
        # train predictor
        self.predictor=MLPRegressor(
                alpha=hyper_params.alpha,
                hidden_layer_sizes=hyper_params.hidden_layer_sizes,
                activation=hyper_params.activation,
                learning_rate=hyper_params.learning_rate,
                solver=hyper_params.solver)

        self.predictor.fit(X_train,Y_train)
        
        # estimate prarams for training data
        self.estimate_params(self.training_data)
        self.estimate_params(self.test_data)
        
        # calculate covariance and correlation matrices
        self.covar=np.cov(np.array((self.test_data["H_est"]-self.test_data["H"],self.test_data["sigma_est"]-self.test_data["sigma"])))
        i,j=np.meshgrid(np.arange(self.covar.shape[0]),np.arange(self.covar.shape[1]),indexing="ij")
        self.corr=self.covar/np.sqrt(self.covar[i,i]*self.covar[j,j])
        
        
        
        return
    
    def estimate_params(self,data_table):
        
        """
        
        Subroutine: estiamte parameters from data table
        
        Arguments
        ---------
        
        data_table: pandas data frame
            table which contains X values. Y values will be added to this table
        
        """
        
        X=get_X(data_table)
        Y=self.predictor.predict(X)
        
        data_table["H_est"]=Y[:,0]
        data_table["sigma_est"]=Y[:,1]
        
        return
        