This repository contains the code and all raw- and metadata to perform the calculations and reporduce the figures
from the manuscript:
<b> Buchwald SZ, Birgel D, Senger K, Mosociova T, Pei Y, Frank AB, Galasso F, Gómez Correa MA, Koşun E, Karapunar B, Wang X, Kustatscher E, Prinoth H, Steinkrauss R, Peckmann J, and Foster WJ
<br> "Primary productivity blooms on the Barents Shelf, Svalbard, associated with the Permian–Triassic mass extinction"

</B> To properly run the script, the folder structure must mirror the structure in this repositioy:
The R scripts `Extinction_marker.R` and `Extinction_marker_functions.R` need to be in the working directory.
Create a sub-folder titled `Raw_data` in the working directory, which should contain all raw and metadata, including
`Extinction_marker_functions.xlsx`, `MPI_raw.xlsx` and the folder `area_molar_sieve`, in which the integrated
peak areas of the compounds of interest are saved in an individual `.txt`-files per sample.

When initiating the project, a folder `Output` is created in the working directory,
that will contain all data produced.

Additionally, the output of quantified n-ACHs and phytanyltoluene (in ug/g TOC) is provided in the Excel file `FID_mol_quantified_ug_TOC_Chol.xlsx`. 
