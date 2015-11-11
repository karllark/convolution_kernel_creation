; program to automatically create the irac/mips to irac/mips convolution kernels

pro create_many_kernels,irac=irac,mips=mips,galex=galex,uit=uit,halpha=halpha, $
  temp_str=temp_str,twomass=twomass,scuba=scuba

if (not keyword_set(temp_str)) then temp_str = '100'

To_c1 = 750.
To_c2 = 750.
To_c3 = 750.
To_c4 = 750.

; both nuv and fuv uit have same gaussian?
uit_nuv = 'UIT/uit_nuv_a1.fits'
To_uit = 200.

; galex
galex_nuv = 'GALEX/galex_nuv.fits'
galex_fuv = 'GALEX/galex_fuv.fits'
To_galex = 150.

; ha for M81
ha_m81 = 'Ha/ha_m81.fits'
To_ha = 200.

; 2mass
twomass_j = '2MASS/2mass_j.fits'
twomass_h = '2MASS/2mass_h.fits'
twomass_k = '2MASS/2mass_k.fits'
To_2mass = 200.

;radio_m81 = 'radio_m81.fits'

irac_path = 'IRAC/'
irac_c1 = irac_path+'irac_ch1_flight.fits'
irac_c2 = irac_path+'irac_ch2_flight.fits'
irac_c3 = irac_path+'irac_ch3_flight.fits'
irac_c4 = irac_path+'irac_ch4_flight.fits'
mips_path = 'MIPS/'
mips_24 = mips_path + 'mips_24_'+temp_str+'K.fits'
mips_70 = mips_path + 'mips_70_'+temp_str+'K.fits'
mips_160 = mips_path + 'mips_160_'+temp_str+'K.fits'
;mips_path = 'stinytim/'
;mips_24 = mips_path + 'mips24_'+temp_str+'K_sub10_smooth16.fits'
;mips_70 = mips_path + 'mips70_'+temp_str+'K_sub10_smooth13.fits'
;mips_160 = mips_path + 'mips160_'+temp_str+'K_sub10_smooth18.fits'

; scuba
To_scuba_450 = 100.
To_scuba_850 = 75.
scuba_path = 'SCUBA/'
scuba_450 = scuba_path + 'scuba_450.fits'
scuba_850 = scuba_path + 'scuba_850.fits'

if (keyword_set(uit)) then begin
    create_kernel,uit_nuv,mips_24,'uit_nuv_to_mips_24_'+temp_str+'K.fits',To=To_uit
    create_kernel,uit_nuv,mips_70,'uit_nuv_to_mips_70_'+temp_str+'K.fits',To=To_uit
    create_kernel,uit_nuv,mips_160,'uit_nuv_to_mips_160_'+temp_str+'K.fits',To=To_uit
endif

if (keyword_set(galex)) then begin
;    create_kernel,galex_nuv,mips_24,'galex_nuv_to_mips_24_75K.fits',To=To_galex
    create_kernel,galex_nuv,mips_70,'galex_nuv_to_mips_70_'+temp_str+'K.fits',To=To_galex
    create_kernel,galex_nuv,mips_160,'galex_nuv_to_mips_160_'+temp_str+'K.fits',To=To_galex
;    create_kernel,galex_fuv,mips_24,'galex_fuv_to_mips_24_75K.fits',To=To_galex
    create_kernel,galex_fuv,mips_70,'galex_fuv_to_mips_70_'+temp_str+'K.fits',To=To_galex
    create_kernel,galex_fuv,mips_160,'galex_fuv_to_mips_160_'+temp_str+'K.fits',To=To_galex
endif

if (keyword_set(halpha)) then begin
    create_kernel,ha_m81,mips_24,'ha_m81_to_mips_24_'+temp_str+'K.fits',To=To_ha
    create_kernel,ha_m81,mips_70,'ha_m81_to_mips_70_'+temp_str+'K.fits',To=To_ha
    create_kernel,ha_m81,mips_160,'ha_m81_to_mips_160_'+temp_str+'K.fits',To=To_ha
endif

;create_kernel,radio_m81,mips_160,'radio_m81_to_mips_160_50K.fits'

if (keyword_set(twomass)) then begin
    create_kernel,twomass_j,mips_24,'2mass_j_to_mips_24_'+temp_str+'K.fits',To=To_2mass
    create_kernel,twomass_j,mips_70,'2mass_j_to_mips_70_'+temp_str+'K.fits',To=To_2mass
    create_kernel,twomass_j,mips_160,'2mass_j_to_mips_160_'+temp_str+'K.fits',To=To_2mass

    create_kernel,twomass_h,mips_24,'2mass_h_to_mips_24_'+temp_str+'K.fits',To=To_2mass
    create_kernel,twomass_h,mips_70,'2mass_h_to_mips_70_'+temp_str+'K.fits',To=To_2mass
    create_kernel,twomass_h,mips_160,'2mass_h_to_mips_160_'+temp_str+'K.fits',To=To_2mass

    create_kernel,twomass_k,mips_24,'2mass_k_to_mips_24_'+temp_str+'K.fits',To=To_2mass
    create_kernel,twomass_k,mips_70,'2mass_k_to_mips_70_'+temp_str+'K.fits',To=To_2mass
    create_kernel,twomass_k,mips_160,'2mass_k_to_mips_160_'+temp_str+'K.fits',To=To_2mass
endif

if (keyword_set(irac)) then begin
    create_kernel,irac_c1,irac_c2,'irac_c1_to_c2.fits',To=30.
    create_kernel,irac_c1,irac_c3,'irac_c1_to_c3.fits',To=35.
    create_kernel,irac_c1,irac_c4,'irac_c1_to_c4.fits',To=35.
    create_kernel,irac_c2,irac_c3,'irac_c2_to_c3.fits',To=25.
    create_kernel,irac_c2,irac_c4,'irac_c2_to_c4.fits',To=35.
    create_kernel,irac_c3,irac_c4,'irac_c3_to_c4.fits',To=35.

    if (keyword_set(mips)) then begin
        create_kernel,irac_c1,mips_24,'irac_c1_to_mips_24_'+temp_str+'K.fits',To=To_c1
        create_kernel,irac_c1,mips_70,'irac_c1_to_mips_70_'+temp_str+'K.fits',To=To_c1
        create_kernel,irac_c1,mips_160,'irac_c1_to_mips_160_'+temp_str+'K.fits',To=To_c1
        create_kernel,irac_c2,mips_24,'irac_c2_to_mips_24_'+temp_str+'K.fits',To=To_c2
        create_kernel,irac_c2,mips_70,'irac_c2_to_mips_70_'+temp_str+'K.fits',To=To_c2
        create_kernel,irac_c2,mips_160,'irac_c2_to_mips_160_'+temp_str+'K.fits',To=To_c2
        create_kernel,irac_c3,mips_24,'irac_c3_to_mips_24_'+temp_str+'K.fits',To=To_c3
        create_kernel,irac_c3,mips_70,'irac_c3_to_mips_70_'+temp_str+'K.fits',To=To_c3
        create_kernel,irac_c3,mips_160,'irac_c3_to_mips_160_'+temp_str+'K.fits',To=To_c3
        create_kernel,irac_c4,mips_24,'irac_c4_to_mips_24_'+temp_str+'K.fits',To=To_c4
        create_kernel,irac_c4,mips_70,'irac_c4_to_mips_70_'+temp_str+'K.fits',To=To_c4
        create_kernel,irac_c4,mips_160,'irac_c4_to_mips_160_'+temp_str+'K.fits',To=To_c4
    endif 
endif

if (keyword_set(mips)) then begin
    create_kernel,mips_24,mips_70,'mips_24_to_70_'+temp_str+'K.fits',To=225.
    create_kernel,mips_24,mips_160,'mips_24_to_160_'+temp_str+'K.fits',To=200.
    create_kernel,mips_70,mips_160,'mips_70_to_160_'+temp_str+'K.fits',To=85.
endif

if (keyword_set(scuba)) then begin
    create_kernel,scuba_450,mips_160,'scuba_450_to_mips_160_'+temp_str+'K.fits',To=To_scuba_450
    create_kernel,scuba_850,mips_160,'scuba_850_to_mips_160_'+temp_str+'K.fits',To=To_scuba_850
endif

end
