# SpectraAlign
Combine Spectra from different ranges by aligning them and finding a smooth overlap. Call it with two 2-column ASCII files containing the spectra

    spectraaligner spectrum_1.asc spectrum_2.asc > combined_spectrum.asc
It does not matter, if spectrum_1 oder spectrum_2 is the spectrum with higher wavenumbers/energies/wavelengths. The spectra just need to go from high to low values and have to overlap by at least one grid point.

