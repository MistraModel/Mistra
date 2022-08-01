;*************************************************
; Joyce_2014_Fig4.ncl
;************************************************
;
; Plot O3, N2O5, ClNO2, HNO3 (gaz)
;  and NO3-, NH4+ (aqueous)
;  in a similar way as Joyce et al 2014 (Fig 4)
;
; Author: J. Bock
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

  outdir = "~/Mistra_2019/out/Joyce2014bc/"

  GridFile = addfile(outdir+"grid.nc","r")
  MeteoFile = addfile(outdir+"meteo.nc","r")
  GasFile = addfile(outdir+"gas.nc","r")
  AqFile = addfile(outdir+"aq.nc","r")

; define bounds
; -------------
; plot up to ~600m: end at grid index 62 (remember index 1 is skipped (zero-height layer), index 2 middle of layer = 5m, index 62 = 605m)
  gidx = 62
; start at t=0
  tidx = 0

; set the axis variables: height and time
; ---------------------------------------
  height = GridFile->eta(0,1:gidx,0,0) ; file is (rec=1, n=150, y=1, x=1) ; units m
  height!0 = "height"
  height@positive = "up"

; create a time array, that will be used to redefine time dimension of all arrays
  nday = GasFile->lday(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units day
  nhour = GasFile->lst(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units hour
  nmin = GasFile->lmin(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units min
  time = nday*24. + nhour + nmin/60.


; import data
; -----------
  ; rho, to convert units
  rhoT = MeteoFile->rho(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units kg/m3
  rho = transpose(rhoT)
; convert to m3/mol
  M_air=28.96546e-3 ; kg/mol
  xair = M_air / rho

  ; ozone
  O3T = GasFile->O3(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units mol/m3
  O3 = transpose(O3T)
; convert to nmol/mol
  O3 = O3 * xair * 1e9

  ; N2O5
  N2O5T = GasFile->N2O5(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units mol/m3
  N2O5 = transpose(N2O5T)
; convert to pmol/mol
  N2O5 = N2O5 * xair * 1e12

  ; ClNO2
  ClNO2T = GasFile->ClNO2(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units mol/m3
  ClNO2 = transpose(ClNO2T)
; convert to pmol/mol
  ClNO2 = ClNO2 * xair * 1e12

  ; NO3-
  NO3mT = AqFile->NO3m(tidx:,1:gidx,0,0) \\ ; file is (rec=, n=150, y=1, nkc_l=2) ; units mol/m3
        + AqFile->NO3m(tidx:,1:gidx,0,1)
  NO3m = transpose(NO3mT)
; convert to mu g/mol
  M_no3 = 62e-3 ; kg/mol
  NO3m = NO3m * M_no3 * 1e9

  ; HNO3
  HNO3T = GasFile->HNO3(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units mol/m3
  HNO3 = transpose(HNO3T)
; convert to nmol/mol
  HNO3 = HNO3 * xair * 1e9

  ; NH4+
  NH4pT = AqFile->NH4p(tidx:,1:gidx,0,0) \\ ; file is (rec=, n=150, y=1, nkc_l=2) ; units mol/m3
        + AqFile->NH4p(tidx:,1:gidx,0,1)
  NH4p = transpose(NH4pT)
; convert to mu g/mol
  M_nh4 = 20e-3 ; kg/mol
  NH4p = NH4p * M_nh4 * 1e9


; set axis to data variables
; --------------------------
  O3!0 = "height"
  O3!1 = "time"
  O3&height = height
  O3&time = time

  N2O5!0 = "height"
  N2O5!1 = "time"
  N2O5&height = height
  N2O5&time = time

  ClNO2!0 = "height"
  ClNO2!1 = "time"
  ClNO2&height = height
  ClNO2&time = time

  NO3m!0 = "height"
  NO3m!1 = "time"
  NO3m&height = height
  NO3m&time = time

  HNO3!0 = "height"
  HNO3!1 = "time"
  HNO3&height = height
  HNO3&time = time

  NH4p!0 = "height"
  NH4p!1 = "time"
  NH4p&height = height
  NH4p&time = time



;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("pdf","Joyce_2014_Fig4")          ; send graphics to PDF file
  plot = new(6,graphic)

  res                              = True
  res@gsnDraw                      = False
  res@gsnFrame                     = False
  res@gsnMaximize                  = True            ; use full page
  res@gsnPaperOrientation          = "portrait"      ; prevent automatic rotation
  res@gsnLeftStringFontHeightF     = 0.03
  res@gsnRightStringFontHeightF    = 0.03

  ; adjust plot to get same height and width
  res@vpHeightF           = 0.6
  res@vpWidthF            = 0.6

  ; Contour levels settings
  res@cnFillOn             = True            ; turn on color to fill areas between contour levels
  res@cnFillPalette        = "WhiteBlueGreenYellowRed"
  res@cnLinesOn            = False
  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnMinLevelValF       =  2.
  res@cnMaxLevelValF       =  38.
  res@cnLevelSpacingF      =  2.
  res@lbOrientation        = "Vertical"  ; Default is horizontal

  ; Axis settings
  res@tmXTOn               = False
  res@tmXBLabelFontHeightF = 0.030
  res@tmYLLabelFontHeightF = 0.030
  res@tmXBMinorOn          = False
  res@tmYLMode             = "Explicit"
  res@tmYLValues           = (/5,50,100,150,200,250,300,350,400,450,500,550/)
  res@tmYLLabels           = "" + res@tmYLValues
  res@trYMinF              = 5.
  res@trYMaxF              = 600.
  res@tiXAxisString     = "Time (h)"
  res@tiYAxisString     = "Height (m)"
  res@tiXAxisFontHeightF = 0.03
  res@tiYAxisFontHeightF = 0.03
  ; Informations displayed on top
  res@gsnLeftString = "O~B~3~N~ (nmol mol~S~-1~N~)"
  res@gsnRightString = "max="+max(O3)

  plot(0) = gsn_csm_contour(wks,O3,res)

  ; Contour levels settings
  res@cnMinLevelValF       =  10.
  res@cnMaxLevelValF       =  565.
  res@cnLevelSpacingF      =  40.
  res@gsnLeftString = "N~B~2~N~O~B~5~N~ (pmol mol~S~-1~N~)"
  res@gsnRightString = "max="+sprintf("%5.1f",max(N2O5))
  plot(1) = gsn_csm_contour(wks,N2O5,res)

  res@cnMinLevelValF       =  10.
  res@cnMaxLevelValF       =  50.
  res@cnLevelSpacingF      =  4.
  res@gsnLeftString = "ClNO~B~2~N~ (pmol mol~S~-1~N~)"
  res@gsnRightString = "max="+sprintf("%5.1f",max(ClNO2))
  plot(2) = gsn_csm_contour(wks,ClNO2,res)

  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnMinLevelValF       =  0.4
  res@cnMaxLevelValF       =  6.0
  res@cnLevelSpacingF      =  0.4
  res@gsnLeftString = "NO~B~3~N~~S~-~N~ (~F33~m~F~g m~S~-3~N~)"
  res@gsnRightString = "max="+sprintf("%5.1f",max(NO3m))
  plot(3) = gsn_csm_contour(wks,NO3m,res)

  res@cnMinLevelValF       =  0.1
  res@cnMaxLevelValF       =  2.0
  res@cnLevelSpacingF      =  0.2
  res@gsnLeftString = "HNO~B~3~N~ (nmol mol~S~-1~N~)"
  res@gsnRightString = "max="+sprintf("%5.1f",max(HNO3))
  plot(4) = gsn_csm_contour(wks,HNO3,res)

  res@cnMinLevelValF       =  0.04
  res@cnMaxLevelValF       =  2.0
  res@cnLevelSpacingF      =  0.2
  res@gsnLeftString = "NH~B~4~N~~S~+~N~ (~F33~m~F~g m~S~-3~N~)"
  res@gsnRightString = "max="+sprintf("%5.1f",max(NH4p))
  plot(5) = gsn_csm_contour(wks,NH4p,res)

;************************************************
; create panel
;************************************************
  resP = True
  resP@gsnPanelXWhiteSpacePercent = 3
  resP@gsnPanelYWhiteSpacePercent = 4
  gsn_panel(wks,plot,(/3,2/),resP)             ; now draw as one plot

end
