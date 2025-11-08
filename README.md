# pcMRI-MATLAB-pipeline
This MATLAB pipeline processes phase-contrast MRI (pcMRI) data to analyze brain tissue motion.
It performs the following steps:

Velocity inhomogeneity correction – removes background bias in the velocity field.

Semi-automatic ROI segmentation – lets the user select a region of interest interactively.

Velocity waveform extraction – computes the mean velocity within the ROI across cardiac phases.

Displacement computation – integrates the velocity waveform to obtain displacement.

Quantitative analysis – calculates

the temporal maximum displacement per voxel, and

the spatially maximum displacement within the ROI.

📄 Citation

If you use or adapt this code, please cite:

Karamzadeh, M., et al. (2025).
Cardiac-induced brain tissue motion in Chiari Malformation type 1 subjects and its relationship to symptomatology, morphometrics, and surgical outcomes.
Magnetic Resonance Imaging, 2025.
https://www.sciencedirect.com/science/article/pii/S0730725X25002140
