# LLNL_ToFi

LLNL\_ToFi is a small python code for tomographic filtering of hypothetical seismic mantle structure (vs or vp) using the resolution matrix of the LLNL-G3D-JPS model by Simmons et al. (2015). The routine LLNL_ToFi.py performs the matrix-vector multiplication R\*m=m' to obtain the filtered version m' of the given seismic model m. To be able to perform this operation, m needs to be given in the parametrization of the LLNL-G3D-JPS model.

Main author: Bernhard Schuberth (Geophysics, LMU Munich, Germany, bernhard.schuberth@lmu.de)

Co-author: Nathan Simmons (Lawrence Livermore National Laboratory, CA, USA, simmons27@llnl.gov)

## License

[![License: GLPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

## Development

Development is hosted on the GitLab server of the Leibniz Supercomputing Centre (LRZ) in the
[bschuberth/LLNL_ToFi repository](https://gitlab.lrz.de/bschuberth/LLNL_ToFi).

## Data
Input data required by LLNL_ToFi.py are located on the LLNL server:

  * [Data](https://www-gs.llnl.gov/nuclear-threat-reduction/nuclear-explosion-monitoring/global-3d-seismic-tomography)

or by email request to simmons27@llnl.gov

## Documentation

  * Information on the resolution matrix and the parametrization of the LLNL-G3D-JPS tomographic model can be found in _LLNL-G3D-JPS\_R\_Matrix\_TomoFilter\_README.pdf_

  * The LLNL-G3D-JPS model is described in 
    _Simmons, N. A., S. C. Myers, G. Johannesson, E. Matzel, and S. P. Grand (2015), Evidence for long-lived subduction of an ancient tectonic plate beneath the southern Indian Ocean, Geophys. Res. Lett., 42, doi:10.1002/2015GL066237._

  * An example of applying the resolution matrix R to a geodynamic model is described in 
  _N. A. Simmons, B. S. A. Schuberth, S. C. Myers, D. R. Knapp (2019), Resolution and Covariance of the LLNL-G3D-JPS Global Seismic Tomography Model: Applications to Travel Time Uncertainty and Tomographic Filtering of Geodynamic Models, Geophys. J. Int., in press_

  * Running the code:
    1. To run the code, please first get the necessary input data (i.e., the resolution matrix files) from the source given above and put them into the directory _./DATA_.

    2. Edit the file _model.py_: Modify the code such that the dummy method _project\_model\_3D_ returns the values of your specific seismic velocity model at the given coordinates (radius, lat, lon).

    3. Run the code in a python 2.7 environment with: _python LLNL\_ToFi.py_ 

       Specify the option _-r|--reparam_ if you run the code the first time. This will perform the reparametrization of your seismic model onto the parametrization of the LLNL-G3D-JPS tomographic model.
       
       Output files will be stored in the directory _./OUTPUT\_FiLES_. Files containing the reparametrized model will be named according to the variable _OUTFILE\_PARM\_PREFIX_ [default:_LLNL\_G3D\_JPS\_Parm\_layer_] and the tomographically filtered model will be stored according to the variable _OUTFILE\_FILT\_PREFIX_ [default:_LLNL\_G3D\_JPS\_ToFi\_layer_].


