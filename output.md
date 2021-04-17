
# Output files

Once calling the subroutuine [saving](https://github.com/PrincetonUniversity/FOCUS/blob/master/sources/saving.f90), present calculating results will be saved to local files.
The frequency of writing outputs is controlled by the variable *save_freq*.
The basic output file is written in hdf5 format.

* hdf5 file
  
  All the results can be seen in *focus_example.h5*. 
  Any hdf5 file reading functions can be used for reading the results. 
  
If you want to view other output files, please turn on the flags in the namelist.

* save_coils

  Optimized coils will be stored in *example.focus* and *example.coils*. 
  They are FOCUS format coils and the standard coils file mentioned above, respectively. 
  So they can be used as the initial coils of next optimization.
  
* save_harmonics

  The Bn spectrum produced by current coils are stored in *example.harmonics*.
  
* save_filaments

Intermediate coil data (XYZ points in space) are stored in binary file *.example.filaments.000001*. 
  So you can use them to plot coil evolution movie. (Better way is to use the data in hdf5 file.)

Various post-processing options are available and they are controlled by `case_postproc`.
Here is a list of available functions.

| post_proc | actions                                                              |
|:---------:|----------------------------------------------------------------------|
|     0     | None                                                                 |
|     1     | Basic post-processing                                                |
|     2     | Basic post-processing + SPEC input preparation                       |
|     3     | Basic post-processing + Poincare plots                               |
|     4     | Basic post-processing + Poincare plots + Boozer spectrum calculation |
|     5     | Basic post-processing + Write MGRID file                             |

&nbsp;

# Plotting
There are several tools for processing the data. If you need to use one of them, please contact the author(s).

  - **Python:** There is an [python package](https://zhucaoxiang.github.io/CoilPy/) written by Dr. Caoxiang Zhu, using [Matplotlib](https://matplotlib.org/) and [Mayavi](http://docs.enthought.com/mayavi/mayavi/).
  - **OMFIT:** In the [OMFIT](https://gafusion.github.io/OMFIT-source/) frame, there is a powerful module [*focus*](https://docs.google.com/document/d/1aGpRUMpYxBJmQXfkOK2OFMZ4P0Mn5H1_OSAjJVvxpHU/edit?ts=5ad10d28#) managed by Dr. Nikolas Logan.
  - **IDL:** There is a GUI interface Echidna in IDL written by Dr. Stuart Hudson.
  - **Matlab:** There are also some MATLAB scripts by Dr. Samuel Lazerson.  
