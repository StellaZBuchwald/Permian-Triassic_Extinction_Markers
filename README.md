This repository contains the code and all raw- and metadata to perform the calculations and reproduce the figures
from the manuscript:
<b> Buchwald SZ, Birgel D, Senger K, Mosociova T, Pei Y, Zuchuat V, Tarhan LG, Frank AB, Galasso F, Gómez Correa MA, Koşun E, Karapunar B, Wang X, Kustatscher E, Prinoth H, Lahajnar N, Steinkrauss R, Peckmann J, and Foster WJ
<br> "Phytoplankton blooms on the Barents Shelf, Svalbard, associated with the Permian–Triassic mass extinction"

</B> To properly run the script, the folder structure must mirror the structure in this repository:
The R scripts `Extinction_marker.R` and `Extinction_marker_functions.R` need to be in the working directory.
Create a sub-folder titled `Raw_data` in the working directory, which should contain all raw- and metadata, including
`Extinction_marker_metadata.xlsx`, `MPI_raw.xlsx`, `Nabbefeld.xlsx` and the folders `area_mol_sieve` and `area_n_alk`, in which the integrated
peak areas of the compounds of interest are saved in an individual `.txt`-files per sample.

When initiating the project, a folder `Output` is created in the working directory,
that will contain all output produced.

Additionally, the output of quantified n-ACHs and phytanyl toluene (in ug/g TOC) is provided in the Excel file `FID_mol_quantified_ug_TOC_Chol.xlsx`, and quantified n-alkanes in the Excel file `FID_nalk_quantified_ug_TOC_Chol.xlsx`. 
