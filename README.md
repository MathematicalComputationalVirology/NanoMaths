# NanoMaths
The three files should be read into MATLAB.
Files called "cage_10_analysis_final_with_theta.m", "cage_14_analysis_final_with_theta.m" and "cage_16_analysis_final_with_theta" are for KWOCAs cages:
In each code the range of S is currently specified according to the experimentaly observed range.
The input is the range of theta, and the output is the corresponding range of delta. In the paper (Fig. 4d) we indicate the value of
delta closest to the experimentally observed range.

.dat files (coordinates of cages) are being generated using the Carbon Generator (CaGe) programme: Brinkmann, G., Delgado-Friedrichs, O., Lisken, S., Peeters, A. & Van Cleemput, N. CaGe: a virtual environment for studying some special classes of plane graphs : an update. MATCH-COMMUNICATIONS IN MATHEMATICAL AND IN COMPUTER CHEMISTRY 63, 533–552 (2010).

Files called "trimer_angle_distribution.m", "finding_minimal_arm_length.m", "plotting_k_range.m" and "trimer_distribution_range_check.m" are for I32-10 cages:

plotting_k_range.m finds minimal frequency (mf) for various scaling factors (k) and plots the results.

trimer_angle_distribution.m: The imput is k and the output is the minimum and maximum arm length, and mf.

finding_minimal_arm_length.m finds the ideal arm length range. In this code we first identify the smallest range of arm lengths containing all the experimentally vable ranges: 0.8394 – 1.0039. In this range the lowest and highest value might be outliers as they would squash, respectively stretch, the trimeric arms the most. If we exclude these outliers we get a smaller ranger, i.e. 0.8627 – 0.9693, but this new range is small as it is not possible to reconstruct some of the observed cages using this range. Thus, the ideal range should be between these two ranges. To find the new range, we interpolate the space between these two ranges  in steps of 0.0058, resulting in 35 intervals. For all these, we then check if it is possible to build the observed cages with positive values of mf. This identifies the unique interval 0.845225 – 0.975067.

trimer_distribution_range_check.m: The input is the minimum and maximum arm lenght and the output are the mf values.


