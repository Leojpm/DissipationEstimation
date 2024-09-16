The latest code for applying this method is DDCdissipationV4. Please feel free to let me know at leo.middleton@whoi.edu 
with any bugs etc. or if you just want to chat more about double diffusion. I'm always happy to explain this method in more detail.

The data as input to this code should be gridded vertically in 
along-track profiles i.e. P should be a single vector for every profile in the matrices T and S.
The profiles do not need to be a constant distance,
apart and the method may give more accurate values if you use the true measured
horizontal distances between profiles as opposed to an artificial 
horizontal grid. Our final output will be on a half grid, as it uses both
vertical and horizontal gradients i.e. if your temperature field is an
(n,m) matrix, your dissipation prediction will be an (n-1,m-1) matrix.

Inputs:
- S : Absolute Salinity on a pressure-time grid, both time and pressure
spacing may vary.
- T : Conservative Temperature on a pressure-time grid
- p : Pressure (dbar)
- lat : Latitude
- lon : Longitude
- kscale : The exponent for the power spectrum of spice. In Middleton et
al. (2021) we used k^(-1) (i.e. kscale = -1), following observations from Mackinnon et al.
(2016) and others. In Fine et al. (2022), we also used k^(-2), as k^(-1)
overestimated observed dissipation rates.

Outputs:
- epsilon : Turbulent Kinetic Energy Dissipation Rate estimate for averaged
epsilon between observations based on the Middleton et al. (2021) method.
- N : Buoyancy frequency, N^2 = -dbdz

Notes:
 - Dissipation rate has only performed well for low buoyancy Reynolds number
(Re_b = epsilon/nu/N^2) data. Testing shows values less that 100 work
well in comparison to measured epsilon: values of 20 have also been used
as thresholds. I recommend varying this threshold in your data if you are
making a comparison to microstructure.

- The spacing between observations is important for this method. As shown
in Middleton et al. (2021), as the spacing increases between
observations, the skill of the parameterization decreases. Recommended
spacing is somewhere between 1 km and 10 km. Larger spacing will still
give an answer, but using a k^(-1) spectrum only will give overinflated
dissipation rates (using k^(-2) or k^(-3) may give something more
accurate).

- kscale is limited to values >-3 due to the form of the integral in the
iterative loop. For values of -3, you should be able to use kscale =
-2.99.

- The method is best applied to the original sampled data. e.g. if you
sampled at a variable grid spacing of [1km 1.5km 0.5km 3km] for four
profiles, then do not interpolate this onto a consistent grid spacing,
rather use the natural sampling resolution. The same applies to glider
data. If you can use the spacing between observations derived from the
original see-saw pattern, you will get better results. For a version of
this code that works for glider data, feel free to send me an email at
leo.middleton@whoi.edu.

- This method is not a substitute for microstructure observations! It may
provide insight when no observations are available, but it will not give
the full picture of turbulent dissipation rate. It can also be used to
compare to microstructure to understand whether given observations could
feasibly be explained by double-diffusive processes. However, this
parameterization does not involve understanding the feedbacks between
shear-driven turbulence and double-diffusive convection, so it will
likely underestimate the integrated effects of double diffusion.
