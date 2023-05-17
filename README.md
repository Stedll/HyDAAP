# **Hy**perspectral to **D**em **A**utomatic **A**lignment and **P**rojection
A general overview of the project is available here [here](HyDAAP.pdf)

## Related projects
### georef_webcam
The main structure of json_writer.py an the general interface with the PRACTISE package has been adapted from the code provided in [georef_webcam](https://github.com/SebBuchelt/georef_webcam)
### PRACTISE
The PRACTISE code in src/PRACTISE to compute the viewshed and project the DEM points in SRF is adapted from the one provided by [georef_webcam](https://github.com/SebBuchelt/georef_webcam) and [PRACTISE](https://github.com/shaerer/PRACTISE)
### SOFT
The SOFT package has been adapted starting from this [unofficial repo](https://gitlab.com/ccpem/soft), the original code was discussed in:  
- Kostelec, P.J., Rockmore, D.N. FFTs on the Rotation Group. J Fourier Anal Appl 14, 145–179 (2008). https://doi.org/10.1007/s00041-008-9013-5

## How to run
The src/main.py file is the main interface for the entire code, in order to run you need to modify the parameters in main.py according to your specific task, such as filenames and the various control switches

### SOFT package
You need to build the hydaap routine in the SOFT submodule using 'make hydaap' inside SOFT/  
To build the package fftw3 >3.3.10 is required, modification can be made removing the multithreading support in order to build with an older version of fftw3  

### Conda environment
The python environment is provided as a conda .yml file that can be easily imported using "conda env create -f environment.yml"  
All necessary packages should be then installed  

### MATLAB
PRACTISE has been adapted from octave to MATLAB for various problems regarding mainly the performances, to correctly run the code you need MATLAB (tested on 2022b) and the [mapping toolbox](https://it.mathworks.com/products/mapping.html)

## Folder structure
```
.
├── data
│   ├── dem
│   │   ├── *merged
│   │   │   └── demfile.tif
│   │   ├── dem1_path
│   │   │   └── dem1.tif
│   │   └── dem2_path
│   │       └── dem2.tif
│   ├── img
│   │   └── hyp
│   │       ├── acquisition.azg
│   │       ├── acquisition.hdr
│   │       └── acquisition.img
│   ├── out
│   │   └── *run_name
│   │       ├── dem_skyline.png
│   │       ├── pattern.dat
│   │       ├── signal.dat
│   │       ├── output.tif
│   │       └── spectral_skyline.png
│   └── tmp
└── src
    ├── PRACTISE
    │   ├── Cam_PRACTISE.m
    │   ├── Input_PRACTISE.m
    │   ├── IsOctave_PRACTISE.m
    │   ├── PRACTISE_spher.m
    │   ├── Proj_PRACTISE_spher.m
    │   ├── Unix_Slash_PRACTISE.m
    │   └── template.m
    ├── dem_functions.py
    ├── dem_reprojector.py
    ├── geotiff_writer.py
    ├── image_functions.py
    ├── json_writer.py
    ├── main.py
    ├── position_optimizer.py
    ├── spherical_functions.py
    └── utils.py
```

Path preceded by '*' are created at runtime if needed
