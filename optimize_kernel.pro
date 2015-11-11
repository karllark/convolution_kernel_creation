;+
; NAME:
;       OPTIMIZE_KERNEL
;
; PURPOSE:
;       This procedure optimizes the To paramter used in creating a
;       kernel using create_kernel.pro.  It is iterative and can take
;       quite a long time!
;
; CATEGORY:
;       Fun.
;
; CALLING SEQUENCE:
;       OPTIMIZE_KERNEL, from_file, to_file, out_file
;
; INPUTS:
;       from_file : FITS file with from PSF (must have )
;       to_file : FITS file with to PSF (must have )
;       out_file : FITS filename to write convolution kernal to (must have)
;       to_fwhm : FWHM of the PSF of the to_file (must have)
;
; KEYWORD PARAMETERS:
;       range_To : range of T_o (2 element vector, default = [50.,500.])
;                  T_o is the clamping period
;       chisqr_tol : tolerance for chisqr [default = 0.01]
;
; OUTPUTS:
;
; RESTRICTIONS:
;       Use at your own risk.  Good luck.
;
;
; PROCEDURE:
;       The code iterates using create_kernel and check_kernel to
;       optimize To.
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
; 	Written by: Karl Gordon (19 Feb 2008)
;       updated: 8 Apr 2008 (KDG)
;-

pro optimize_kernel,from_file,to_file,out_file,to_fwhm,range_To=range_To, $
                    chisqr_tol=chisqr_tol,norm_range=norm_range,use_radial=use_radial

if (n_params() LT 4) then begin
    printf,'optimize_kernel,from_file,to_file,out_file,to_fwhm'
retall & end

if (not keyword_set(range_To)) then range_To = [50.,500.]
if (not keyword_set(chisqr_tol)) then chisqr_tol = 0.01
silent = 1

size_norm = 3.*to_fwhm

; get the to psf
fits_read,to_file,to_psf_image,to_psf_header
getrot,to_psf_header,rot,cdelt
to_psf_cdelt = cdelt

; determine the kernels at the ends and middle of range_To
min_To = range_To[0]
max_To = range_To[1]

print,'creating kernel with T_o = ', min_To
create_kernel,from_file,to_file,'To_min.fits',To=min_To,silent=silent
check_kernel,'To_min.fits',from_file,to_file,min_chisqr,silent=silent,norm_out=norm_out,norm_range=norm_range
print,'radial chisqr = ', min_chisqr
; get image chisqr
fits_read,'test.fits',To_image,To_header
getrot,To_header,rot,cdelt

to_psf_image = change_image_scale(to_psf_image,to_psf_header,abs(to_psf_cdelt[0])*3600.,abs(cdelt[0])*3600.)
size_image = size(to_psf_image)
getrot,to_psf_header,rot,cdelt
hsize = min([size_norm/(abs(cdelt[0])*3600.),(size_image[1]/2)-1])
to_x1 = size_image[1]/2-hsize
to_x2 = size_image[1]/2+hsize
to_y1 = size_image[2]/2-hsize
to_y2 = size_image[2]/2+hsize
to_psf_image /= total(to_psf_image[to_x1:to_x2,to_y1:to_y2])

;To_image = change_image_scale(To_image,To_header,abs(cdelt[0])*3600.,abs(to_psf_cdelt[0])*3600.)
indxs = where(finite(To_image),n_indxs)
size_image = size(To_image)
getrot,To_header,rot,cdelt
hsize = min([size_norm/(abs(cdelt[0])*3600.),(size_image[1]/2)-1])
x1 = size_image[1]/2-hsize
x2 = size_image[1]/2+hsize
y1 = size_image[2]/2-hsize
y2 = size_image[2]/2+hsize

To_image /= total(To_image[x1:x2,y1:y2])
if (not keyword_set(use_radial)) then begin
    min_chisqr = total((To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2])^2)
    print,'image chisqr = ', min_chisqr
endif
fits_write,'To_min_diff.fits',To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2]
; get structure chisqr
To_image[x1:x2,y1:y2] = !values.f_nan
findxs = where(finite(To_image))
struct_chisqr = total((To_image[findxs] - mean(To_image[findxs]))^2)
print,'struct chisqr = ', struct_chisqr
min_chisqr = struct_chisqr

print,'creating kernel with T_o = ', max_To
create_kernel,from_file,to_file,'To_max.fits',To=max_To,silent=silent
check_kernel,'To_max.fits',from_file,to_file,max_chisqr,silent=silent,norm_out=norm_out,norm_range=norm_range
print,'radial chisqr = ', max_chisqr
; get image chisqr
fits_read,'test.fits',To_image,To_header
;To_image = change_image_scale(To_image,To_header,abs(cdelt[0])*3600.,abs(to_psf_cdelt[0])*3600.)
indxs = where(finite(To_image),n_indxs)
To_image /= total(To_image[x1:x2,y1:y2])
if (not keyword_set(use_radial)) then begin
    max_chisqr = total((To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2])^2)
    print,'image chisqr = ', max_chisqr
endif
fits_write,'To_max_diff.fits',To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2]
; get structure chisqr
To_image[x1:x2,y1:y2] = !values.f_nan
findxs = where(finite(To_image))
struct_chisqr = total((To_image[findxs] - mean(To_image[findxs]))^2)
print,'struct chisqr = ', struct_chisqr
max_chisqr = struct_chisqr

mid_To = range_To[0] + (range_To[1] - range_To[0])/2.
print,'creating kernel with T_o = ', mid_To
create_kernel,from_file,to_file,'To_mid.fits',To=mid_To,silent=silent
check_kernel,'To_mid.fits',from_file,to_file,mid_chisqr,silent=silent,norm_out=norm_out,norm_range=norm_range
print,'radial chisqr = ', mid_chisqr
; get image chisqr
fits_read,'test.fits',To_image,To_header
;To_image = change_image_scale(To_image,To_header,abs(cdelt[0])*3600.,abs(to_psf_cdelt[0])*3600.)
indxs = where(finite(To_image),n_indxs)
To_image /= total(To_image[x1:x2,y1:y2])
if (not keyword_set(use_radial)) then begin
    mid_chisqr = total((To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2])^2)
    print,'image chisqr = ', mid_chisqr
endif
fits_write,'To_mid_diff.fits',To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2]
; get structure chisqr
To_image[x1:x2,y1:y2] = !values.f_nan
findxs = where(finite(To_image))
struct_chisqr = total((To_image[findxs] - mean(To_image[findxs]))^2)
print,'struct chisqr = ', struct_chisqr
mid_chisqr = struct_chisqr

done = 0
while (not done) do begin 
    chisqrs = [min_chisqr,mid_chisqr,max_chisqr]
    max_val = max(chisqrs,mindx)
    case mindx of
        0: begin
            min_To = mid_To
            min_chisqr = mid_chisqr
        end
        1: begin
            print,'whoops'
            print,min_To,mid_To,max_To
            print,min_chisqr,mid_chisqr,max_chisqr
            done = 1
        end
        2: begin
            max_To = mid_To
            max_chisqr = mid_chisqr
        end
    endcase

    delta_chisqr = 0.5*abs(max_chisqr - min_chisqr)/(max_chisqr + min_chisqr)
    if (delta_chisqr GT chisqr_tol) then begin
        mid_To = min_To + (max_To - min_To)/2.
        print,'creating kernel with T_o = ', mid_To
        create_kernel,from_file,to_file,'To_mid.fits',To=mid_To,silent=silent
        check_kernel,'To_mid.fits',from_file,to_file,mid_chisqr,silent=silent,norm_out=norm_out,norm_range=norm_range
        print,'radial chisqr = ', mid_chisqr
        ; get image chisqr
        fits_read,'test.fits',To_image,To_header
;        To_image = change_image_scale(To_image,To_header,abs(cdelt[0])*3600.,abs(to_psf_cdelt[0])*3600.)
        indxs = where(finite(To_image),n_indxs)
        To_image /= total(To_image[x1:x2,y1:y2])
        if (not keyword_set(use_radial)) then begin
            mid_chisqr = total((To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2])^2)
        print,'image chisqr = ', mid_chisqr
        endif
        fits_write,'To_mid_diff.fits',To_image[x1:x2,y1:y2] - to_psf_image[to_x1:to_x2,to_y1:to_y2]
; get structure chisqr
        To_image[x1:x2,y1:y2] = !values.f_nan
        findxs = where(finite(To_image))
        struct_chisqr = total((To_image[findxs] - mean(To_image[findxs]))^2)
        print,'struct chisqr = ', struct_chisqr
        mid_chisqr = struct_chisqr
    endif else begin
        done = 1
    endelse
endwhile
chisqrs = [min_chisqr,mid_chisqr,max_chisqr]
print,chisqrs
T_o = [min_To,mid_To,max_To]
print,T_o
min_val = min(chisqrs,mindx)
best_To = T_o[mindx]

print,'final T_o = ', best_To
create_kernel,from_file,to_file,out_file,To=best_To,silent=silent
check_kernel,out_file,from_file,to_file,mid_chisqr,silent=silent

file_delete,'To_min.fits'
file_delete,'To_mid.fits'
file_delete,'To_max.fits'
file_delete,'test.fits'

end
