;*************************************************
; Joyce_2014_Fig2.ncl
;************************************************
;
; Plot specific humidity, T and theta in a similar way as Joyce et al 2014 (Fig 2)
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

; define bounds
; -------------
; plot up to ~600m: end at grid index 62 (remember index 1 is skipped (zero-height layer), index 2 middle of layer = 5m, index 62 = 605m)
  gidx = 101;62
; start at t=0
  tidx = 0

; set the axis variables: height and time
; ---------------------------------------
  height = GridFile->eta(0,:gidx,0,0) ; file is (rec=1, n=150, y=1, x=1) ; units m
  height!0 = "height"
  height@positive = "up"

; create a time array, that will be used to redefine time dimension of all arrays
  nday = MeteoFile->lday(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units day
  nhour = MeteoFile->lst(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units hour
  nmin = MeteoFile->lmin(tidx:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units min
  time = nday*24. + nhour + nmin/60.


; import data
; -----------
  ; specific humidity is variable xm1
  qt = MeteoFile->xm1(tidx:,:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units kg/kg
  q = transpose(qt)
; convert to g/kg
  q = q * 1e3

  ; temperature and potential temperature
  tempt = MeteoFile->temp(tidx:,:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units kg/kg
  temp = transpose(tempt)
  thetat = MeteoFile->theta(tidx:,:gidx,0,0) ; file is (rec=, n=150, y=1, x=1) ; units kg/kg
  theta = transpose(thetat)

; set axis to data variables
; --------------------------
  q!0 = "height"
  q!1 = "time"
  q&height = height
  q&time = time

  temp!0 = "height"
  temp!1 = "time"
  temp&height = height
  temp&time = time

  theta!0 = "height"
  theta!1 = "time"
  theta&height = height
  theta&time = time


;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("pdf","Joyce_2014_Fig2")          ; send graphics to PDF file

  res                              = True
  res@gsnDraw                      = False
  res@gsnFrame                     = False
  res@gsnMaximize                  = True            ; use full page
  res@gsnPaperOrientation          = "portrait"      ; prevent automatic rotation

  ; adjust plot to get same height and width
  res@vpHeightF           = 0.6
  res@vpWidthF            = 0.6

  ; Contour levels settings
  res@cnFillOn             = True            ; turn on color to fill areas between contour levels
  res@cnFillPalette        = "WhiteBlue"
  res@cnLinesOn            = False           ; no lines for humidity contour levels
  res@cnLevelSelectionMode =  "ManualLevels"
  ;res@cnLineLabelBackgroundColor = -1 ; make background transparent (default=0=background)
  res@cnMinLevelValF       =  0.35
  res@cnMaxLevelValF       =  0.8
  res@cnLevelSpacingF      =  0.05
  res@lbOrientation        = "Vertical"  ; Default is horizontal

  ; Axis settings
  res@trYMinF              = 5.
  res@trYMaxF              = 600.
  res@tiXAxisString = "Time (h)"
  res@tiYAxisString = "Height (m)"
  res@tiXAxisFontHeightF = 0.02
  res@tmXBMinorOn = False
  res@tmXTOn      = False
  res@tmYLMode    = "Explicit"
  res@tmYLValues  = (/5,50,100,150,200,250,300,350,400,450,500,550/)
  res@tmYLLabels  = "" + res@tmYLValues
  ; Informations displayed on top: nothing in this plot
  res@gsnLeftString = ""
  res@gsnRightString = ""

  plot1 = gsn_csm_contour(wks,q,res)

  ; Contour levels settings
  res@cnFillOn            = False
  res@cnInfoLabelOn       = False
  res@cnLinesOn           = True           ; turn on the lines for t and theta contour levels
  res@cnLineDashPattern   = 16             ; dashed lines for temperature
  res@cnLineLabelInterval = 1              ; display a label for each line (default=2)
  res@cnLineLabelDensityF = 0.55           ; decrease the number of labels for each line (adjusted to be nicely presented)
  res@cnMinLevelValF      =  250.
  res@cnMaxLevelValF      =  260.
  res@cnLevelSpacingF     =  1.
  plot2 = gsn_csm_contour(wks,temp,res)

  res@cnLineDashPattern    = 0             ; full lines for potential temperature
  res@cnLevelSelectionMode =  "ExplicitLevels"
  res@cnLevels             =  (/255.7,256.4,256.7/)
  plot3 = gsn_csm_contour(wks,theta,res)

  overlay (plot1,plot2)
  overlay (plot1,plot3)
  draw (plot1)
  frame (wks)

end