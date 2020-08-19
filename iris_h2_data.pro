
PRO iris_h2_data, file = file

;*******************************************************************
;
; NAME:
;     IRIS_H2_DATA
;
; PURPOSE:
;     Measument of Doppler and Nonthermal velocities for IRIS H2 line at 1333.79 Angstrom
;     
; CATEGORY:
;     IRIS; measurement of physical parameters from H2 line widths.
;
; INPUTS:
;     Restore .save file and use stored variables as input
;
; VARIABLE NAMEs:
;     The following parameters were derived from the H2 spectral line (single Gaussian) fitting 
;
;     h2_slit_no: IRIS slit numbers
;
;     h2_ut_time: IRIS slit time in UT
;
;     h2_pix: IRIS pixels along the slit which are used to create H2 spectra 
;
;     label: H2 line label
;
;     wavelength: wavelength of H2 line
;
;     shift: line shift determined from the calculation - O I line was used for wavelength calibration
;     more details can be found at IRIS ITN 20 - https://iris.lmsal.com/documents.html
;
;     h2_fwhm: FWHM measured for H2 line (Gaussian fitting)
;
;     h2_fwhm_err: Errors on measured FWHM
;
;     h2_cent: centroid of H2 line
;
;     h2_cent_err: erroe on the centroid positions
;
; OUTPUTS:
;     Returns Doppler and Nonthermal velocities
;
; EXAMPLE:         
;
;     IDL> file = 'iris_flare_h2_parameters.save'
;     IDL> file = file
;
;     IDL> .r iris_h2_data
;     IDL> iris_h2_data, file = file
;
; MODIFICATION HISTORY:
;     Ver.1, 19-August-2020, Dr. Sargam Madhav Mulay
; 
;*************************************************************************

restore, file, /verbose

; RESTORE: Restored variable: H2_SLIT_NO.
; RESTORE: Restored variable: H2_UT_TIME.
; RESTORE: Restored variable: H2_PIX.
; RESTORE: Restored variable: LABEL.
; RESTORE: Restored variable: WAVELENGTH.
; RESTORE: Restored variable: SHIFT.
; RESTORE: Restored variable: H2_FWHM.
; RESTORE: Restored variable: H2_FWHM_ERR.
; RESTORE: Restored variable: H2_CENT.
; RESTORE: Restored variable: H2_CENT_ERR.


;**********************************
; H2 Doppler velocity calculation
;**********************************

ref_wave = wavelength - shift

pos_plus = h2_cent + h2_cent_err
pos_minus = h2_cent - h2_cent_err

h2_dop_vel = (h2_cent - ref_wave) / ref_wave*3d5 
h2_dop_vel_plus = (pos_plus - ref_wave) / ref_wave*3d5
h2_dop_vel_minus = (pos_minus - ref_wave) / ref_wave*3d5

n = n_elements(h2_dop_vel)
h2_dop_vel_avg = fltarr(n)
h2_dop_vel_err = fltarr(n)

for i = 0, n-1 do begin 
aa = [h2_dop_vel_plus[i], h2_dop_vel[i], h2_dop_vel_minus[i]]
h2_dop_vel_avg[i] = average(aa)
h2_dop_vel_err[i] = stddev(aa)
endfor

print, 'Average Doppler velocity = ', h2_dop_vel_avg
print, 'Doppler velocity error = ', h2_dop_vel_err


;**********************************
; H2 thermal width measurement
;**********************************

; Thermal width of molecular hydrogen in the unit of km/s

; element name
element = 'H' 

; mass of H2 atom and H2 molecule
mass = eis_element2mass(element,/kg) 
mol_mass = 2*mass
print, mol_mass ; = 3.34757e-27 in kg

; Temperature of formation for Hydrogen molecule in Kelvin (from Innes, D. 2008, A&A, 481, L41)
ti_use = 4200

; wavelength of molecular hydrogen in Angstrom
wavelength = 1333.79 

; Boltzmann constant k = 1.380 649. 10^-23 J / K
k_const = 1.3806488E-23

; speed of light in km/s
vel_light = 2.9979E+5 

; Thermal width of Molecular Hydrogen in km/s
h2_thermal_width = sqrt(2*(ti_use)*k_const/mol_mass)/1.0E+3
print, 'Thermal width of Molcular Hydrogen in km/s = ', h2_thermal_width ; = 5.88595 km/s

;**********************************
; H2 Nonthermal width measurement
;**********************************

; wavelength of molecular hydrogen in Angstrom
wavelength = 1333.79 

; Speed of light in km/s
vel_light = 2.9979E+5 

; According to the IRIS paper, the spectral resolution (FWHM) is 26 mÅ in the FUV.
instr_fwhm = 0.026 ; in Angstrom for FUV

; instrumental width in the unit of km/s
instr_v = instr_fwhm/sqrt(alog(2.))/2./wavelength*vel_light 

print, 'IRIS instrumental width = ', instr_v

; (delta lambda / lambda) = (v/c)
; observed line width in km/s
h2_obs_width = (h2_fwhm * vel_light / wavelength) 
print, 'Observed width of Molcular Hydrogen km/s = ', h2_obs_width

; Nonthermal FWHM in km/s
h2_nth_width_sq = (h2_obs_width^2 - instr_v^2 - h2_thermal_width^2)
h2_nth_width = sqrt(h2_nth_width_sq)

; width error - plus
width_plus = h2_fwhm + h2_fwhm_err
h2_obs_width_plus = (width_plus * vel_light / wavelength) 
h2_nth_width_sq_plus = (h2_obs_width_plus^2 - instr_v^2 - h2_thermal_width^2)
h2_nth_width_plus = sqrt(h2_nth_width_sq_plus)

; width error - minus
width_minus = h2_fwhm - h2_fwhm_err
h2_obs_width_minus = (width_minus * vel_light / wavelength)
h2_nth_width_sq_minus = (h2_obs_width_minus^2 - instr_v^2 - h2_thermal_width^2)
h2_nth_width_minus = sqrt(abs(h2_nth_width_sq_minus))

n = n_elements(h2_nth_width)
h2_nth_vel_avg = fltarr(n)
h2_nth_vel_err = fltarr(n)

for i = 0, n-1 do begin 
bb = [h2_nth_width_plus[i], h2_nth_width[i], h2_nth_width_minus[i]]
h2_nth_vel_avg[i] = average(bb)
h2_nth_vel_err[i] = stddev(bb)
endfor

print, 'Average Nonthermal velocity = ', h2_nth_vel_avg
print, 'Nonthermal velocity error = ', h2_nth_vel_err

save, filename = 'iris_flare_h2_velocity.save', h2_slit_no, h2_ut_time, h2_pix, label, wavelength, shift, h2_fwhm, h2_fwhm_err, h2_cent, h2_cent_err, h2_dop_vel_avg, h2_dop_vel_err, h2_thermal_width, h2_nth_vel_avg, h2_nth_vel_err


;*****************************************************************


end
