;+
; NAME:
;       CREATE_KERNEL
;
; PURPOSE:
;       This procedure takes two PSFs and creates the convolution
;       kernel which coverts the first PSF to the second PSF.  This
;       program was created to create kernels to take a higher
;       resolution PSF and turn it accurately into a lower resolution
;       PSF.
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;       CREATE_KERNEL, from_file, to_file, out_file
;
; INPUTS:
;       from_file : FITS file with from PSF (must have )
;       to_file : FITS file with to PSF (must have )
;       out_file : FITS filename to write convolution kernal to
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;       This program uses Fast Fourier Transforms (FFTs) to create the
;       new convolution kernel.  The key to getting this method to
;       work is to properly surpress the high frequencies in the FFT
;       of the from PSF.  These high frequencies only contain noise
;       and dominate the convolution kernal if not surpressed.  For
;       details, see "The Fast Fourier Transform and Its Applications"
;       section 14.5 (page 345) by E. Oran Brigham.
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
; 	Written by:	Karl Gordon w/ help from JD Smith (17 Feb 2004)
;           2 Jun 2005 (KDG) : updated to properly center kernel in image
;
;-
pro create_kernel,from_file,to_file,out_file,div_val=div_val,To=To, $
                  show_kernel=show_kernel,silent=silent

; check we got at least 3 parameters
if (N_params() LT 3) then begin
    print, 'Syntax - create_kernel,from_file,to_file,out_file'
    return
endif

psf_path = ''
if (not keyword_set(div_val)) then div_val = 10.

; print current output filename
if (not keyword_set(silent)) then print,'working on ' + out_file + '...'

; get from psf
if (file_test(psf_path+from_file)) then begin
    fits_read,psf_path+from_file,from_psf_image,from_psf_header
    from_psf_scale = fxpar(from_psf_header,'PIXSCALE',count=count)
    if (count EQ 0) then begin ; for IRAC
        from_psf_scale = fxpar(from_psf_header,'SECPIX',count=count)
    endif
    if (count EQ 0) then begin ; for IRAF produced?
        from_psf_scale = fxpar(from_psf_header,'SCALE',count=count)
    endif
    if (count EQ 0) then begin
        getrot,from_psf_header,from_rot,from_cdelt
        from_psf_scale = abs(from_cdelt[0])*3600.
        if (from_psf_scale EQ 0) then count = 0 else count = 1
    endif
    if (count EQ 0) then begin
        from_psf_scale = abs(fxpar(from_psf_header,'CD1_1',count=count))*3600.
        print,from_psf_scale
    endif
    if (from_psf_scale EQ 0) then begin
        print,from_psf_header
        print,'from_psf_scale = 0'
        return
    endif
endif else begin
    print,'From PSF file (' + psf_path+from_file + ') does not exist.'
    return
endelse

; get to psf
if (file_test(psf_path+to_file)) then begin
    fits_read,psf_path+to_file,to_psf_image,to_psf_header
    to_psf_scale = fxpar(to_psf_header,'PIXSCALE',count=count)
    if (count EQ 0) then begin ; for IRAC
        to_psf_scale = fxpar(to_psf_header,'SECPIX',count=count)
    endif
    if (count EQ 0) then begin
        getrot,to_psf_header,to_rot,to_cdelt
        to_psf_scale = abs(to_cdelt[0])*3600.
        if (to_psf_scale EQ 0) then count = 0 else count = 1
    endif
    if (to_psf_scale EQ 0) then begin
        print,'to_psf_scale = 0'
        return
    endif
endif else begin
    print,'To PSF file (' + psf_path+to_file + ') does not exist.'
    return
endelse

size_orig_from_psf_image = size(from_psf_image)

; put the to psf on the same scale as the from psf
to_psf_image = change_image_scale(to_psf_image,to_psf_header,to_psf_scale,from_psf_scale)

; make sure both PSFs are centered
ensure_psf_centered,from_psf_image
ensure_psf_centered,to_psf_image

; now make the two psf images the same size
size_from_psf_image = size(from_psf_image)
size_to_psf_image = size(to_psf_image)

if (size_from_psf_image[1] GT size_to_psf_image[1]) then begin
    new_psf_image = replicate(0.0,size_from_psf_image[1],size_from_psf_image[2])
    i1 = 0 + (size_from_psf_image[1] - size_to_psf_image[1])/2
    i2 = i1 + size_to_psf_image[1] - 1
    j1 = 0 + (size_from_psf_image[2] - size_to_psf_image[2])/2
    j2 = i1 + size_to_psf_image[2] - 1
    new_psf_image[i1:i2,j1:j2] = to_psf_image
    ensure_psf_centered,new_psf_image
    to_psf_image = new_psf_image
endif else begin
    new_psf_image = replicate(0.0,size_to_psf_image[1],size_to_psf_image[2])
    i1 = 0 + (size_to_psf_image[1] - size_from_psf_image[1])/2
    i2 = i1 + size_from_psf_image[1] - 1
    j1 = 0 + (size_to_psf_image[2] - size_from_psf_image[2])/2
    j2 = i1 + size_from_psf_image[2] - 1
    new_psf_image[i1:i2,j1:j2] = from_psf_image
    ensure_psf_centered,new_psf_image
    from_psf_image = new_psf_image
endelse


; create radius image
size_psf_image = size(from_psf_image)
x_npts = size_psf_image[1]
y_npts = size_psf_image[2]
x_vals = (findgen(x_npts) + 0.5 - x_npts/2)#replicate(1.0,y_npts)
y_vals = replicate(1.0,x_npts)#(findgen(y_npts) + 0.5 - y_npts/2)
rad_image = sqrt(x_vals^2 + y_vals^2)

; using FFTs, make the convolution kernel
size_psf = size(from_psf_image)
;print,size_psf
shift_x = size_psf[1]/2 + (size_psf[1] mod 2)
shift_y = size_psf[2]/2 + (size_psf[2] mod 2)
npix = size_psf[1]*size_psf[2]

;print,'getting ffts of psfs'
; get PSF FFTs
from_psf_fft = fft(from_psf_image,-1)
to_psf_fft = fft(to_psf_image,-1)

; determine 1/from_FFT making sure to cut all high frequencies
abs_from_psf_fft = abs(from_psf_fft)

inverse_from_psf_fft = 1.0/from_psf_fft

; clip if asked
if (keyword_set(To)) then begin
    if (not keyword_set(silent)) then print,'clipping from psf fft'
;    To = size_orig_from_psf_image[1]/2
    if (not keyword_set(silent)) then print,'To = ',To
; determine cutting function (Hanning)
    Wh = 0.5*(1.0 + cos(2*!PI*rad_image/To))
;    Wh = fltarr(x_npts,y_npts)
    indxs = where(rad_image GE To/2.0,n_indxs)
    if (n_indxs GT 0) then Wh[indxs] = 0.0
;    Wh = (To/2.)*(sin(!PI*freq*To)/(!PI*freq*To(1.- (freq*To)^2)))

;    fits_write,'wh.fits',Wh

;    cut_val = max(abs_from_psf_fft)/div_val
;    indxs = where(abs_from_psf_fft LT cut_val,n_indxs)

;help,inverse_from_psf_fft
;abs_from_psf_fft[indxs] = 0.0
    inverse_from_psf_fft = inverse_from_psf_fft*shift(Wh,shift_x,shift_y)
;    from_psf_image *= Wh

;    fits_write,'from_psf_wh.fits',from_psf_image
endif

;fits_write,'onefft_abs.fits',shift(abs_from_psf_fft,shift_x,shift_y)
;fits_write,'one.fits',float(fft(from_psf_fft,1))
;from_psf_fft[indxs] = 0.0
;fits_write,'one_clamp.fits',float(fft(from_psf_fft,1))

;fits_write,'oneoverfft_abs.fits',abs(inverse_from_psf_fft)
;fits_write,'oneoverfft_fl.fits',float(inverse_from_psf_fft)
;fits_write,'oneoverfft_im.fits',imaginary(inverse_from_psf_fft)
;if (n_indxs GT 0) then inverse_from_psf_fft[indxs] = 0.0
;fits_write,'oneoverfft_clamp_abs.fits',abs(inverse_from_psf_fft)
;fits_write,'oneoverfft_clamp_fl.fits',float(inverse_from_psf_fft)
;fits_write,'oneoverfft_clamp_im.fits',imaginary(inverse_from_psf_fft)

; make conv kernel
if (not keyword_set(silent)) then print,'making kernel'
;fits_write,'test.fits',float(to_psf_fft)
;fits_write,'test2.fits',float(from_psf_fft)
conv_kernel_fft = to_psf_fft*inverse_from_psf_fft
conv_kernel = npix*shift(float(fft(conv_kernel_fft,1)),shift_x,shift_y)

; ensure psf is centered
ensure_psf_centered,conv_kernel

; check if the kernel is too large and rebin smaller if so
size_conv_kernel = size(conv_kernel)
max_kernel_size = 1501
if (size_conv_kernel[1] GT max_kernel_size) then begin
    old_from_psf_scale = from_psf_scale
    from_psf_scale = from_psf_scale*(float(size_conv_kernel[1])/max_kernel_size)
    sxaddpar,from_psf_header,'NAXIS1',size_conv_kernel[1]
    sxaddpar,from_psf_header,'NAXIS2',size_conv_kernel[1]
    conv_kernel = change_image_scale(conv_kernel,from_psf_header,old_from_psf_scale,from_psf_scale)
    if (not keyword_set(silent)) then print,'changing final kernel from ', old_from_psf_scale, ' to ', from_psf_scale
endif

; normalize kernel to have a total of 1
conv_kernel = conv_kernel/total(conv_kernel)

;print,conv_kernel[88:92,88:92]/max(conv_kernel)

; make output header
fxhmake,out_header,conv_kernel
sxaddpar,out_header,'PIXSCALE',from_psf_scale
sxaddpar,out_header,'CD1_1',from_psf_scale/3600.0
sxaddpar,out_header,'CD1_2',0.0
sxaddpar,out_header,'CD2_1',0.0
sxaddpar,out_header,'CD2_2',from_psf_scale/3600.0
sxaddpar,out_header,'HISTORY','Custom Convolution Kernel'
sxaddpar,out_header,'HISTORY','Created by create_kernel.pro (IDL program)'
sxaddpar,out_header,'HISTORY','written by Karl D. Gordon (kgordon@stsci.edu)'
sxaddpar,out_header,'HISTORY','From PSF file = ' + from_file
sxaddpar,out_header,'HISTORY','To PSF file = ' + to_file
if (keyword_set(To)) then begin
    sxaddpar,out_header,'T_o',To,' Hanning function T_o value used'
endif

; write conv kernel to the specified FITS files
fits_write,psf_path+out_file,conv_kernel,out_header

if (keyword_set(show_kernel)) then begin
; display kernel
    create_png,psf_path+out_file,'tmp.png',/ss_scale,scaled_image=scaled_image

    size_image = size(scaled_image)
    xwindow = 600.0
    scale_fac = xwindow/size_image[1]
    if (scale_fac LT 1) then begin
        scale_fac = fix(1.0/scale_fac)
        scaled_image = rebin(scaled_image,size_image[1]/scale_fac,size_image[2]/scale_fac)
    endif else begin
        scale_fac = fix(scale_fac)
        scaled_image = rebin(scaled_image,size_image[1]*scale_fac,size_image[2]*scale_fac)
    endelse
    size_image = size(scaled_image)
    
    window,0,xsize=size_image[1],ysize=size_image[2]
    loadct,3,/silent
    tv,scaled_image
endif

end
