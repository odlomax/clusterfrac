# Default ClusterFrac parameters

class hyper_params:
    
    # fitting hyperparameters
    alpha=0.0001
    hidden_layer_sizes=(40,)
    activation="tanh"
    learning_rate="adaptive"
    solver="lbfgs"
    data_split=0.7
    
class model_2d2d_params:
    
    # parameters for 2d model
    n_dim_field=2
    n_dim_cluster=2
    grid_size=1024
    gauss_filter=True
    L_range=(0.0005,0.005)
    H_range=(0.,1.)
    sigma_range=(0.2,2.5)
    
class model_3d2d_params:
    
    # parameters for 3d model projected to 2d
    n_dim_field=3
    n_dim_cluster=2
    grid_size=128
    gauss_filter=False
    L_range=(0.0005,0.005)
    H_range=(0.,1.)
    sigma_range=(0.5,3.5)
    
class model_3d3d_params:
    
    # parameters for fully 3d model
    n_dim_field=3
    n_dim_cluster=3
    grid_size=128
    gauss_filter=False
    L_range=(0.0005,0.005)
    H_range=(0.,1.)
    sigma_range=(0.5,3.5)
    