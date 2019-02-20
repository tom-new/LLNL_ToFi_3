# LLNL_ToFi

LLNL\_ToFi is a small python program for tomographic filtering of hypothetical seismic mantle structure (vs or vp) using the resolution matrix of the LLNL-G3D-JPS model by Simmons et al. (2015). The routine LLNL_ToFi.py performs the matrix-vector multiplication R\*m=m' to obtain the filtered version m' of the given seismic model m. To be able to perform this operation, m needs to be given in the parametrization of the LLNL-G3D-JPS model.

Main author: Bernhard Schuberth (Geophysics, LMU Munich, Germany, bernhard.schuberth@lmu.de)

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)

## Development

Development is hosted on the GitLab server of the Leibniz Supercomputing Centre (LRZ) in the
[bschuberth/LLNL_ToFi repository](https://gitlab.lrz.de/bschuberth/LLNL_ToFi).

## Data
Input data required by LLNL_ToFi.py are located on the LLNL server:

  * [Data](https://www-gs.llnl.gov/nuclear-threat-reduction/nuclear-explosion-monitoring/global-3d-seismic-tomography)

or by email request to simmons27@llnl.gov

## Documentation

  * Information on the resolution matrix and the parametrization of the LLNL-G3D-JPS tomographic model can be found in *LLNL-G3D-JPS\_R\_Matrix\_TomoFilter\_README.pdf*

  * The LLNL-G3D-JPS model is described in  
    *Simmons, N. A., S. C. Myers, G. Johannesson, E. Matzel, and S. P. Grand (2015), Evidence for long-lived subduction of an ancient tectonic plate beneath the southern Indian Ocean, Geophys. Res. Lett., 42, doi:10.1002/2015GL066237.*

  * An example of applying the resolution matrix R to a geodynamic model is described in  
  *N. A. Simmons, B. S. A. Schuberth, S. C. Myers, D. R. Knapp (2019), Resolution and Covariance of the LLNL-G3D-JPS Global Seismic Tomography Model: Applications to Travel Time Uncertainty and Tomographic Filtering of Geodynamic Models, Geophys. J. Int., in press*

  * Running the code:
    1. To run the code, please first get the necessary input data (i.e., the resolution matrix files) from the source given above and put them into the directory *./DATA*.

    2. Edit the file *model.py*: Modify the code such that the dummy routine *project\_model\_3D* returns the values of your specific seismic velocity model at the given coordinates (radius, lat, lon).

    3. Run the code in a python 2.7 environment using one of the following commands: 
        * *python LLNL\_ToFi.py*
        * *./LLNL\_ToFi.py*  

        for serial processing, or:

        * *mpirun -np {number of processes} LLNL\_Tofi.py*  
        
        for parallel processing.

       Specify the option *-n|--no-reparam* if you run the code the several times and you do not want to perform the reparametrization again. This assumes that the reparametrized version of your seismic model (i.e., on the parametrization of the LLNL-G3D-JPS tomographic model) is already stored in the directory *OUTPUT\_FILES*.
       
       Output files will be stored in the directory *./OUTPUT\_FiLES*. Files containing the reparametrized model will be named according to the variable *OUTFILE\_PARM\_PREFIX* [default: '*LLNL\_G3D\_JPS\_Parm\_layer*'] and the tomographically filtered model will be stored according to the variable *OUTFILE\_FILT\_PREFIX* [default: '*LLNL\_G3D\_JPS\_ToFi\_layer*'].


