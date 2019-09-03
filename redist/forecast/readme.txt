FORECAST: Fourier-based Off-REsonanCe Artifact simulation in the STeady-state
Version: 1.05
Date: 19/01/2017


FORECAST is a fast alternative to Bloch simulation for simulating
off-resonance effects in steady-state MRI. FORECAST accelerates simulation
of steady-state pulse sequences by using multiple Fast Fourier Transforms
to evaluate the signal equation, which can include proton density, T2, and
off-resonance effects. Currently the simulation is limited to Cartesian
pulse sequences, but we plan to add support for non-Cartesian pulse
sequences as well.

Please note that FORECAST uses a lefthanded notation for phase evolution.
If you would like to compare your simulations with scans in a righthanded
notation (for example, scans from Philips systems), then you need to take
the complex conjugate of either the scan or the simulated image.


Getting started

Place the contents of the .zip file in a directory of your choice (for
example 'c:\Matlab\Forecast'). You can run the demo_gradientecho.m file
to see if our example simulations run properly. If you would like to use
FORECAST from other MATLAB scripts, add the FORECAST directory to your
MATLAB path.


Usage examples

See demo_gradientecho.m, demo_spinecho.m, and demo_brainweb.m for examples
on how to use FORECAST.


Referring to FORECAST

If you use FORECAST in your research, please include a reference to our
MRM paper and a link to the most recent code:
F. Zijlstra, J.G. Bouwman, I. Braskute, M.A. Viergever, and P.R. Seevinck,
"Fast Fourier-based simulation of off-resonance artifacts in steady-state
gradient echo MRI applied to metal object localization", Magn. Reson.
Med., 2016


Changelog

1.01:
- Minor fix in brainweb example.

1.02:
- Fixed a small bug in calculateCartesianSamplingTimes.

1.03:
- Fixed a bug in the brainweb example that caused it to not use the T2 values

1.04:
- Removed dependencies on the Image Processing Toolbox

1.05:
- Included support for spin echo sequences
- Updated reference to our recently published paper


Contact

If you have any questions, suggestions, or find any bugs, feel free to
contact us.

Frank Zijlstra (f.zijlstra@gmail.com) and Job Bouwman (jgbouwman@hotmail.com), 2017
