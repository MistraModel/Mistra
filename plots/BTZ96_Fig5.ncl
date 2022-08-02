;*************************************************
; Fig5_BTZ96.ncl
;************************************************
;
; Plot cloud water content in a similar way as BTZ96
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

  outdir = "~/Mistra_2019/out/BTZ96/"
  GridFile = addfile(outdir+"grid.nc","r")
  MeteoFile = addfile(outdir+"meteo.nc","r")

; plot up to ~900m: end at grid index 92 (remember index 1 is skipped (zero-height layer), index 2 middle of layer = 5m, index 92 = 905m)
  gidx = 92

  height = GridFile->eta(0,:gidx,0,0) ; file is (rec=1, n=150, y=1, x=1) ; units m
  height!0 = "height"
  height@positive = "up"

; skip first day of simulation, spinup: start with time index 96
  tidx = 96

  LWCt = MeteoFile->xm2(tidx:,:gidx,0,0) ; file is (rec=1, n=150, y=1, x=1) ; units kg/m3
  LWC = transpose(LWCt)
; convert to g/m3
  LWC = LWC * 1e3

; create a time array, that will be used to redefine time dimension of all arrays
  nday = MeteoFile->lday(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units day
  nhour = MeteoFile->lst(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units hour
  nmin = MeteoFile->lmin(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units min
  time = (nday-1)*24. + nhour + nmin/60.

  LWC!0 = "height"
  LWC!1 = "time"
  LWC&height = height
  LWC&time = time


;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("pdf","Fig5_BTZ96")          ; send graphics to PDF file

  res                              = True
  res@gsnMaximize                  = True            ; use full page
  res@gsnPaperOrientation          = "portrait"      ; prevent automatic rotation

  ; adjust plot to get same height and width
  res@vpHeightF           = 0.6
  res@vpWidthF            = 0.6

  ; Contour levels settings
  res@cnFillOn            = True            ; turn on color to fill areas between contour levels
  res@cnFillPalette       = "WhiteBlue"
  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnMinLevelValF       =  0.0
  res@cnMaxLevelValF       =  0.7
  res@cnLevelSpacingF      =  0.1

  res@lbLabelAutoStride = True

  ; Axis settings
  res@trYMinF              = 0.
  res@trYMaxF              = 900.
  res@tiXAxisString     = "time (hours)"
  res@tiXAxisOffsetYF = 0.21
  res@tiYAxisString     = "height (m)"
  ; Informations displayed on top
  res@gsnLeftString = "LWC (g m~S~-3~N~)"
  res@gsnRightString = "max = "+ sprintf("%5.2f",max(LWC))

  plot = gsn_csm_contour(wks,LWC,res)


end
