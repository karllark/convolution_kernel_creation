; program to check how good a kernel is by
; convolving the from PSF with the kernel and comparing
; it to the to PSF using radial profiles (to start with)

pro check_kernel,kernel_file,from_file,to_file,chisqr,ylog=ylog,paper=paper,eps=eps, $
                 show_plot=show_plot,silent=silent,norm_out=norm_out,norm_range=norm_range

norm_psf = 1
inner_range = [0,2]
rad_npts = 500

fits_read,kernel_file,kern_image,kern_header
fits_read,from_file,from_image,from_header
fits_read,to_file,to_image,to_header

size_kern = size(kern_image)
radave_psf,kernel_file,size_kern[1:2]/2.,kern_rad,kern_psf,kern_psf_unc,kern_psf_med,kern_psf_npts, $
           norm_psf=norm_psf,inner_range=inner_range,rad_npts=rad_npts,silent=silent

size_from = size(from_image)
radave_psf,from_file,size_from[1:2]/2.,from_rad,from_psf,from_psf_unc,from_psf_med,from_psf_npts, $
           norm_psf=norm_psf,inner_range=inner_range,rad_npts=rad_npts,silent=silent

size_to = size(to_image)
radave_psf,to_file,size_to[1:2]/2.,to_rad,to_psf,to_psf_unc,to_psf_med,to_psf_npts, $
           norm_psf=norm_psf,inner_range=inner_range,rad_npts=rad_npts,silent=silent

; now do the real check
test_file = 'test.fits'
conv_image,from_file,kernel_file,test_file,silent=silent

fits_read,test_file,test_image,test_header
size_test = size(test_image)
radave_psf,test_file,size_test[1:2]/2.,test_rad,test_psf,test_psf_unc,test_psf_med,test_psf_npts, $
           norm_psf=norm_psf,inner_range=inner_range,rad_npts=rad_npts,silent=silent

; renormalize the PSFs
norm_pix = 2.
if (not keyword_set(norm_range)) then norm_range = [5.,20.]

;indxs = where(kern_psf_npts GT 0)
;sindxs = sort(abs(kern_rad[indxs] - norm_pix))
;kern_psf /= kern_psf[sindxs[indxs[0]]]
indxs = where((kern_psf_npts GT 0) AND ((kern_rad GT norm_range[0]) AND (kern_rad LT norm_range[1])))
kern_psf /= mean(kern_psf[indxs])

;indxs = where(from_psf_npts GT 0)
;sindxs = sort(abs(from_rad[indxs] - norm_pix))
;from_psf /= from_psf[sindxs[indxs[0]]]
indxs = where((from_psf_npts GT 0) AND ((from_rad GT norm_range[0]) AND (from_rad LT norm_range[1])))
from_psf /= mean(from_psf[indxs])

;indxs = where(to_psf_npts GT 0)
;sindxs = sort(abs(to_rad[indxs] - norm_pix))
;to_psf /= to_psf[sindxs[indxs[0]]]
indxs = where((to_psf_npts GT 0) AND ((to_rad GT norm_range[0]) AND (to_rad LT norm_range[1])))
norm_to = mean(to_psf[indxs])
to_psf /= norm_to

;indxs = where(test_psf_npts GT 0)
;sindxs = sort(abs(test_rad[indxs] - norm_pix))
;test_psf /= test_psf[sindxs[indxs[0]]]
indxs = where((test_psf_npts GT 0) AND ((test_rad GT norm_range[0]) AND (test_rad LT norm_range[1])))
norm_test = mean(test_psf[indxs])
test_psf /= norm_test

norm_out = [norm_to,norm_test]

; determine the chisqr

; interpolate the test_psf to the to_psf radial scale
indxs = where(test_psf_npts GT 0)
test_psf_chisqr = interpol(test_psf[indxs],test_rad[indxs],to_rad)
indxs = where((to_psf_npts GT 0) AND $
              ((to_rad GT norm_range[0]) AND (to_rad LT norm_range[1])),n_indxs)
save_indxs = indxs
if (n_indxs GT 0) then begin
    chisqr = total((abs(test_psf_chisqr[indxs] - to_psf[indxs])^2))/(n_indxs-1)
endif else begin
    print,'no points in common between the two psfs'
    stop
endelse

if (keyword_set(show_plot)) then begin
    ps_file = 'check_kernel'
    setup_ps_output,ps_file,eps=eps

    setup_colors,base_color,back_color,blue_color,red_color,green_color, $
                 yel_color,purple_color,light_blue_color, $
                 line_color=line_color,red=red,green=green,blue=blue, $
                 bw=bw
    

    if (keyword_set(ylog)) then kyplot_type = 'o' else kyplot_type = 'i'
    kxrange = krange([kern_rad,from_rad,to_rad],kplot_type='o')
    kxrange = [1.,100.]
    kyrange = krange([kern_psf,from_psf,to_psf],kplot_type=kyplot_type)
    off_val = 0.0               ;1e-5
;    if (keyword_set(ylog)) then kyrange = [5e-4,1.2] else kyrange = [0.0,1.2]
    if (keyword_set(ylog)) then kyrange = [5e-4,10.] else kyrange = [0.0,2.0]
    kplot,[1],[1],/no_data,xrange=kxrange,yrange=kyrange,color=base_color,kplot_type='o' + kyplot_type, $
          xtitle='radius [arcsec]',ytitle='normalized amplitude'

    koplot,to_rad[save_indxs],to_psf[save_indxs],psym=1,color=base_color
    
    indxs = where(kern_psf_npts GT 0,n_indxs)
    min_kern = min(kern_psf[indxs])
    if (min_kern GT 0) then min_kern = 0.0 else min_kern = abs(min_kern) + off_val
    if (keyword_set(ylog)) then kern_psf[indxs] += min_kern
    koplot,kern_rad[indxs],kern_psf[indxs],psym=100,color=green_color,linestyle=2
    
    indxs = where(from_psf_npts GT 0,n_indxs)
    if (keyword_set(ylog)) then from_psf[indxs] += min_kern
    koplot,from_rad[indxs],from_psf[indxs],psym=100,color=blue_color,linestyle=1

    indxs = where(to_psf_npts GT 0,n_indxs)
    if (keyword_set(ylog)) then to_psf[indxs] += min_kern
    koplot,to_rad[indxs],to_psf[indxs],psym=100,color=base_color,linestyle=0

    indxs = where(test_psf_npts GT 0,n_indxs)
    if (keyword_set(ylog)) then test_psf[indxs] += min_kern
    koplot,test_rad[indxs],test_psf[indxs],psym=100,color=red_color,linestyle=3
    
    ktype = [[0,100],[1,100],[2,100],[3,100]]
    if (keyword_set(paper)) then begin
        kstr = ['MIPS 70','MIPS 24','24->70 kernel','MIPS 24 with kernel']
    endif else begin
        kstr = ['to PSF','from PSF','from->to kernel','from PSF with kernel']
    endelse
    kcolor = [base_color,blue_color,green_color,red_color]
    
    klegend,[0.15,0.4],ktype,kstr,line_color=kcolor,color=base_color,kplot_type='o' + kyplot_type,/box, $
            charsize=1.2
    close_ps_output,eps=eps
endif

end
