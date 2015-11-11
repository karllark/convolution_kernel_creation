;+
; NAME:
;       CREATE_GAUSS_PSF
;
; PURPOSE:
;       Procedure to create a FITS file with a gaussian PSF.  Useful
;       if the FWHM of a gaussian is all that is know for a PSF.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;       CREATE_GAUSS_PSF,fwhm,out_file
;
;
; INPUTS:
;      fwhm : full width, half maximum (FWHM) of gaussian
;      out_file : name of output FITS file
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Written by        : Karl Gordon (18 Feb 2004)
;       21 Sep 2006 (KDG) : fixed gaussian equation (2 was inside ^2)
;-

pro create_gauss_psf,fwhm,out_file,subsamp=subsamp

; check we got at least 2 parameters
if (N_params() LT 2) then begin
    print, 'Syntax - create_gauss_psf,gauss_fwhm,out_file'
    return
endif

; make the size of the psf much larger than the fwhm while 
; sampling the psf well
if (not keyword_set(subsamp)) then subsamp = 10
sub_sample_factor = subsamp
image_scale = fwhm/float(sub_sample_factor)
x_npts = 10*sub_sample_factor + 1
y_npts = x_npts

; make a radius image
x_vals = (findgen(x_npts) + 0.5 - x_npts/2.0)#replicate(1.0,y_npts)
y_vals = replicate(1.0,x_npts)#(findgen(y_npts) + 0.5 - y_npts/2.0)
rad2_image = (image_scale*x_vals)^2 + (image_scale*y_vals)^2

; make psf
psf_image = exp(-rad2_image/(2.*(fwhm/2.35)^2))

; ensure the total flux in PSF is equal to 1
psf_image /= total(psf_image)

; make header
fxhmake,out_header,psf_image
sxaddpar,out_header,'PIXSCALE',image_scale
sxaddpar,out_header,'CD1_1',image_scale/3600.0
sxaddpar,out_header,'CD1_2',0.0
sxaddpar,out_header,'CD2_1',0.0
sxaddpar,out_header,'CD2_2',image_scale/3600.0
sxaddpar,out_header,'CRPIX1',x_npts/2
sxaddpar,out_header,'CRPIX2',y_npts/2
sxaddpar,out_header,'CRVAL1',0.0
sxaddpar,out_header,'CRVAL2',0.0
sxaddpar,out_header,'CTYPE1','RA---TAN'
sxaddpar,out_header,'CTYPE2','DEC--TAN'
sxaddpar,out_header,'HISTORY','Gaussian PSF with fwhm = ' + strtrim(string(fwhm),2)
sxaddpar,out_header,'HISTORY','Created by create_gauss_psf.pro (IDL program)'
sxaddpar,out_header,'HISTORY','written by Karl D. Gordon (kgordon@as.arizona.edu)'

; output file
fits_write,out_file,psf_image,out_header

end
