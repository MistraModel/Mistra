;*************************************************
; Joyce_2014_Fig3.ncl
;************************************************
;
; Plot NOx (gaz)
;  and PM2.5 (aqueous)
;  in a similar way as Joyce et al 2014 (Fig 3)
;
; Author: Josué Bock
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

  ntim=dimsizes(time) ; timesteps

; import data
; -----------
  ; rho, to convert units
  rhoT = MeteoFile->rho(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units kg/m3
  rho = transpose(rhoT)
; convert to m3/mol
  M_air=28.96546e-3 ; kg/mol
  xair = M_air / rho

  ; NOx
  NOT = GasFile->NO(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units mol/m3
  NO2T = GasFile->NO2(tidx:,1:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units mol/m3
  NOxT = NO2T + NOT
  NOx = transpose(NOxT)
; convert to nmol/mol
  NOx = NOx * xair * 1e9

  ; PM2.5: approximated as chemical bin 1 content.
; the water content is NOT accounted for in PM2.5 mass, only major ions
  ;wat = AqFile->cw(:,1:gidx,0,0) * 1e12 ; file is (rec=289, n=150, y=1, nkc_l=2) ; units m3 m-3 converted ug m-3
  NO3m = AqFile->NO3m(:,1:gidx,0,0) * 62e6 ; file is (rec=289, n=150, y=1, nkc_l=2) ; units mol m-3 converted ug m-3
  NH4p = AqFile->NH4p(:,1:gidx,0,0) * 18e6 ; ...
  NH3 = AqFile->NH3(:,1:gidx,0,0) * 17e6 ; ...
  SO42m = AqFile->SO42m(:,1:gidx,0,0) * 96e6 ; ...
  HSO4m = AqFile->HSO4m(:,1:gidx,0,0) * 97e6 ; ...
  Clm = AqFile->Clm(:,1:gidx,0,0) * 35.45e6 ; ...
  DOM = AqFile->DOM(:,1:gidx,0,0) * 122.21e6 ; ... mass is octene, cf. initm values of PJ

  sumliqT = NO3m+NH4p+SO42m+HSO4m+Clm+DOM
  PM2p5 = new((/gidx,ntim/),float)
  PM2p5!0 = "height"
  PM2p5!1 = "time"
  PM2p5&height = height
  PM2p5&time = time

  PM2p5 = (/transpose(sumliqT)/)


; set axis to data variables
; --------------------------
  NOx!0 = "height"
  NOx!1 = "time"
  NOx&height = height
  NOx&time = time


;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("pdf","Joyce_2014_Fig3")          ; send graphics to PDF file
  plot = new(2,graphic)

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
  res@cnMinLevelValF       =  1.
  res@cnMaxLevelValF       =  59.
  res@cnLevelSpacingF      =  4.
  res@lbOrientation        = "Vertical"  ; Default is horizontal

  ; Axis settings
  res@trYMinF = 5.
  res@trYMaxF = 600.
  res@tmXTOn               = False
  res@tmXBLabelFontHeightF = 0.03
  res@tmXBMinorOn          = False
  res@tmYLMode             = "Explicit"
  res@tmYLValues           = (/5,50,100,150,200,250,300,350,400,450,500,550/)
  res@tmYLLabels           = "" + res@tmYLValues
  res@tmYLLabelFontHeightF = 0.03
  res@tiXAxisString      = "Time (h)"
  res@tiYAxisString      = "Height (m)"
  res@tiXAxisFontHeightF = 0.03
  res@tiYAxisFontHeightF = 0.03
  ; Informations displayed on top
  res@gsnLeftString = "NO~B~x~N~ (nmol mol~S~-1~N~)"
  res@gsnRightString = "max = "+sprintf("%5.1f",max(NOx))

  plot(0) = gsn_csm_contour(wks,NOx,res)

  ; Contour levels settings
  res@cnMinLevelValF       =  1.
  res@cnMaxLevelValF       =  19.
  res@cnLevelSpacingF      =  1.
  res@gsnLeftString = "Total PM~B~2.5~N~ (~F33~m~F~g m~S~-3~N~)"
  res@gsnRightString = "max = "+sprintf("%5.1f",max(PM2p5))
  plot(1) = gsn_csm_contour(wks,PM2p5,res)


;************************************************
; create panel
;************************************************
  resP = True
  resP@gsnPanelYWhiteSpacePercent = 2
  gsn_panel(wks,plot,(/2,1/),resP)             ; now draw as one plot

end
