# miglocplot

This package plots various types of seimograms and DE location results obtained by using
migloc3d package (Gharti et al 2010). Data sources 
for this program can be of several types: 1) Synthetic data (binary SAC files) computed using 3d viscoelastic 
finite difference code E3D (Larsen and Grieger, 1998), 2) Synthetic data computed using 
SPECFEM3D (Komatitsch and Tromp, 1998), 3) MAT file processed with microseismic monitoring 
package MIMO (Oye and Roth, 2003), 4) SEG2 files, 5) SEGY files, and 6) MAT file containing 
DE result of migloc3d . If you want to process other types of data, you can look at the file 
load_data.m add/modify as necessary.

Please see doc/miglocplot_manual for more detail.
