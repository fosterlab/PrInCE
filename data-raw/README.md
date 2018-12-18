## Raw data

This directory contains raw data that is subsequently processed into datasets bundled with the `PrInCE` R package, following Hadley Wickham's [guide](http://r-pkgs.had.co.nz/data.html) to including data in packages. 

The example dataset that is included with PrInCE by default is taken from Scott et al., _Mol. Syst. Biol._ 2017 (doi: 10.15252/msb.20167067). The replicate included with the package is the heavy isotope channel from replicate #1 of the cytosolic interactome (i.e., the SEC-PCP-SILAC condition); it was chosen because it is the smallest of the 12 co-elution matrices generated in this study after compression.

To ensure complete reproducibility, the code included in this directory generates the data matrix that is bundled with the R package from an Excel spreadsheet that is included in the supplementary materials associated with the paper, and which is available for download from the [journal website](http://msb.embopress.org/content/13/1/906). 

The package also includes a set of 477 "gold standard" complexes obtained from the EBI [Complex Portal](ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/homo_sapiens.tsv) (downloaded December 17, 2018). These are included to demonstrate the prediction functionality. 

