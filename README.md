# Variable-Time-Normalization

A python implementation of the analysis done in [this](https://github.com/hnagib/Variable-Time-Normalization/blob/master/papers/Bur-s-2016-Angewandte_Chemie_International_Edition.pdf) paper. The data used for testing model were obtained from the paper's [supplementary materials](https://github.com/hnagib/Variable-Time-Normalization/blob/master/papers/anie201609757-sup-0001-misc_information.pdf).

Check out the example usage demo notebook [here](https://nbviewer.jupyter.org/github/hnagib/Variable-Time-Normalization/blob/master/notebooks/01-hnagib-oop.ipynb?flush_cache=True). 

# Python vs. Excel implementation
While the VTN method is simple enough that it can be easily implemented in Excel, the Python implementation provides convenient methods for plotting and algorithmic determination of reaction orders. The author of the paper suggests using visual inspection of the plots to determine reaction orders. I have posed this as a total variation minimization problem. This allows us to use mathematical optimization to determine reaction order. The module developed here provides methods for:

- Computing: `tv()`, 
- Minimizing: `min_tv()`
- Plotting: `plot_tv()`, `plot_min_tv()` 
