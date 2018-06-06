# grain-properties
This Software estimates granular properties such as tortuosity, porosity and co-ordination number from granular particles position data from discrete element method softwares or tomography data.
Tortuosity_Porosity.py does the following things:

    1) creates cross sectional jpeg files from the positional data of particles.
    
    2) tiff files are created from the above jpeg files and trimmed to the required dimension.
    
    3) porosity of each tiff file is calculated by calculating the number of dark pixels and the total number of pixels.
    
    4) Call Perl and R scripts for estimating tortuosity using programs developed by Professor Nakashima.
    
    5) Renames the result files.
    
    6) Cleans up the unwanted files.
    
The python script is in itself well documented, if one reads the python script one can understand the logic used.
The same cross sectioning can be used to estimate tortuosity from X-Ray tomography data also. 

This software uses Tortuosity estimator developed by Professor Nakashima Yoshito Home Page: <https://staff.aist.go.jp/nakashima.yoshito/index-e.htm>. Link to the Program download <https://staff.aist.go.jp/nakashima.yoshito/progeng.htm>. Tortuosity estimator uses random walk method to estimate effective diffusivity from which tortuosity is estinated. The citation for the original paper is followed.
1) Nakashima, Y. and Kamiya, S. (2007) Mathematica Programs for the Analysis of Three-Dimensional Pore Connectivity and Anisotropic Tortuosity of Porous Rocks using X-ray Computed Tomography Image Data. Journal of Nuclear Science and Technology, 44, 1233-1247. 
2) Nakashima, Y. and Nakano, T. (2012) Steady-State Local Diffusive Fluxes in Porous Geo-Materials Obtained by Pore-Scale Simulations. Transport in Porous Media 93, 657-673.

Features yet to be implemented:
Co-ordination number calculation.
