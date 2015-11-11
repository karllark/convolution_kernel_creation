; program to make radial average via median
;
; updated to be more general - KDG 29 May 2006

pro make_ave_radial,file,coords,rad,rad_vals,rad_vals_unc,xrange=xrange,yrange=yrange, $
                    xtitle=xtitle,ytitle=ytitle,mult_fac=mult_fac, $
                    max_rad=max_rad,rad_npts=rad_npts,silent=silent, $
                    ps=ps,eps=eps,bw=bw,show_plot=show_plot,save_rad=save_rad, $
                    exclude_regions_file=exclude_regions_file

if (not keyword_set(max_rad)) then max_rad = 150.
if (not keyword_set(rad_npts)) then rad_npts = 100

if (not keyword_set(silent)) then print,file
fits_read,file,image,header

if (keyword_set(mult_fac)) then begin
    image = image*mult_fac
endif

if (keyword_set(exclude_regions_file)) then begin
    get_exclude_regions,exclude_regions_file,exclude_regions
    nan_dce_body,image,header,exclude_regions=exclude_regions
endif

extast,header,ast_info
ad2xy,coords[0],coords[1],ast_info,pix_x,pix_y
getrot,ast_info,rot,cdelt
pixscale = abs(cdelt[0])*60.0
if (not keyword_set(silent)) then print,'pix_scale [arcmin/pixel] = ',pixscale

size_image = size(image)
x_npts = size_image[1]
y_npts = size_image[2]
x_vals = (findgen(x_npts) + 0.5 - pix_x)#replicate(1.0,y_npts)
y_vals = replicate(1.0,x_npts)#(findgen(y_npts) + 0.5 - pix_y)
rad_image = sqrt(x_vals^2 + y_vals^2)*pixscale

;fits_write,'test.fits',rad_image

delta_rad = max_rad/rad_npts
rad_min = 0.0 + findgen(rad_npts)*delta_rad
rad_max = rad_min + delta_rad
rad = (rad_min + rad_max)/2.0
rad_vals = fltarr(rad_npts)
rad_vals_unc = fltarr(rad_npts)
for i = 0,(rad_npts-1) do begin
;    print,i
    indxs = where((rad_image GE rad_min[i]) AND (rad_image LT rad_max[i]),n_indxs)
;    print,rad_min[i],rad_max[i]
;    print,image[indxs]
;    print,i,n_indxs
    if (n_indxs GT 0) then begin
        indxs2 = where((image[indxs] NE 0.0) AND (finite(image[indxs]) EQ 1),n_indxs)
    endif
;    print,i,n_indxs
    if (n_indxs GT 0) then begin
        indxs = indxs[indxs2]

        iter_sigma_clip,image,indxs,sigma_vals=sigma_vals,/silent
;        rad_vals[i] = median(image[indxs])
        rad_vals[i] = sigma_vals[0]
        rad_vals_unc[i] = sigma_vals[1]
    endif
endfor

if (keyword_set(show_plot)) then begin
    output_filename = strmid(file,0,strlen(file)-5) + '_rad'
    setup_ps_output,output_filename,ps=ps,eps=eps,bw=bw

    setup_colors,base_color,back_color,blue_color,red_color,green_color, $
                 yel_color,purple_color,light_blue_color, $
                 line_color=line_color,red=red,green=green,blue=blue, $
                 bw=bw

    xplottype = 'o'
    yplottype = 'o'
    if (not keyword_set(xrange)) then xrange = krange(rad,kplot_type=xplottype)
    if (not keyword_set(yrange)) then yrange = krange(rad_vals,kplot_type=yplottype)

    kplot,[1],[1],/nodata,kplot_type=xplottype+yplottype, $
          xrange=xrange,yrange=yrange, $
          xtitle=xtitle,ytitle=ytitle, $
          color=base_color,background=back_color
    
    koplot,rad,rad_vals,psym=1,color=red_color
;koplot,rad_image[indxs],image[indxs],psym=21,color=blue_color

    close_ps_output,ps=ps,eps=eps
endif

; write radial average to a ASCII file

if (keyword_set(save_rad)) then begin
    filebase = save_rad
;    filebase = remove_fits_str(file)
    openw,unit1,filebase + '.rad.dat',/get_lun
    printf,unit1,'# radial average'
    printf,unit1,'# FITS file = ' + file
    printf,unit1,'# radius[arcmin], median val'
    for i = 0,(rad_npts-1) do begin
        printf,unit1,rad[i],rad_vals[i]
    endfor
    free_lun,unit1
endif

end
