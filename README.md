# Variable-Time-Normalization

A python implementation of the analysis done in [this](https://github.com/hnagib/Variable-Time-Normalization/blob/master/papers/Bur-s-2016-Angewandte_Chemie_International_Edition.pdf) paper. The data used for testing model were obtained from [here](https://github.com/hnagib/Variable-Time-Normalization/blob/master/papers/anie201609757-sup-0001-misc_information.pdf).

While the VTN method is simple enough that it can be easily implemented in Excel, the Python implementation provides convenient methods for plotting and algorithmic determination of reaction orders. 

The author suggests using visual inspection of the plots to determine reaction orders. This can be a posed as a total variation minimization problem. The module developed here provides methods for computing, plotting and minimizing this objective function. 
