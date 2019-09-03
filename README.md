# Introduction

This MATLAB toolbox contains code to perform localization of metal objects in MRI using fast simulation and template matching. For given scan parameters and a model of the metal object, a library of templates is simulated for many different orientations of the object. Then, to find this type of object in an MR image, the library is loaded and matched using Phase-Only Cross Correlation.

The base methods and brachytherapy seed application are described in:
- Zijlstra, F., Bouwman, J. G., Braškute, I. , Viergever, M. A. and Seevinck, P. R. (2017), Fast Fourier-based simulation of off-resonance artifacts in steady-state gradient echo MRI applied to metal object localization. Magn. Reson. Med, 78: 2035-2041. doi:10.1002/mrm.26556
- Zijlstra, F., Moerland, M. A., Voort van Zyp, J. R., Noteboom, J. L., Viergever, M. A. and Seevinck, P. R. (2017), Challenges in MR-only seed localization for postimplant dosimetry in permanent prostate brachytherapy. Med. Phys., 44: 5051-5060. doi:10.1002/mp.12505

The gold fiducial application is described in:
- Maspero M., van den Berg, C.A.T., Zijlstra, F., Sikkes, G.G., de Boer H.C.J., Meijer G.J., Kerkmeijer L.G.W., Viergever M.A., Lagendijk J.J.W. and Seevinck, P.R. (2017), Evaluation of an automatic MR-based gold fiducial marker localisation method for MR-only prostate radiotherapy, Phys. Med. Biol., 62. doi:10.1088/1361-6560/aa875f 


# Usage

Three examples are included, which show the basic usage of the simulation and template matching code. All examples used a similar multi-echo gradient echo MRI scan. The examples are located in the following directories:
- `cylinder`: A phantom containing a titanium cylinder, as shown in the FORECAST paper.
- `brachyseeds`: A phantom containing 3 brachytherapy seeds in three orientations. Although only one brachytherapy seed type is used, the example includes simulation of a second seed type which was used in the brachytherapy seed
- `goldfiducials`: A phantom containing 3 gold fiducials in three orientations.

Run `download_data.m` to download the example data. Each example contains a `run_createLibrary` and `run_localization` script. Run both to see a visualization of the localization of the object(s).

The object model is defined in scanner space (Z dimension points along the main magnetic field), and is transformed to the image space before simulation. The object model is also responsible for calculating the magnetic field shift caused by the object. We use the calculateFieldShift function by Job Bouwman (included in this repository).

In the main directory there is a `run_test.m` script, that validates that the detected AP/RL/FH dimensions are correct by simulating a template that overlays the letters on an image (`testObjectFunction.m`). Note that dimensions may be flipped, this is a known issue.

# Included packages

- FORECAST: https://www.mathworks.com/matlabcentral/fileexchange/56680-mri-simulation-using-forecast-fourier-based-off-resonance-artifact-simulation-in-the-steady-state
- Dicomseries reader: https://github.com/FrankZijlstra/dicomseries-matlab
- Forward field-shift calculation for MRI: https://www.mathworks.com/matlabcentral/fileexchange/37278-forward-field-shift-calculation-for-mri

# Known issues
- Patient position is not used, it is assumed to be HFS.
- Scan angulation is not supported, scans must be axis-aligned.
- Flips in image directions are not correctly processed.
- Scan type is always assumed to be gradient echo by `getImageParametersFromDicomPhilips.m`.