;*************************************************
; Buys_2013_Fig3.ncl
;************************************************
;
; Plot Br2 and BrCl
;  in a similar way as Buys et al 2014 (Fig 3)
;
; Author: Josué Bock
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

  outdir = "~/Mistra_2019/out/Buys13_0D/"

  GridFile = addfile(outdir+"grid.nc","r")
  MeteoFile = addfile(outdir+"meteo.nc","r")
  GasFile = addfile(outdir+"gas.nc","r")


; define bounds
; -------------
  gidx = 1 ; BOX: output is in layer 2 (index 1 in NCL)

; create a time array, that will be used to redefine time dimension of all arrays
  nday  = GasFile->lday(:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units day
  nhour = GasFile->lst(:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units hour
  nmin  = GasFile->lmin(:,0,0,0) ; file is (rec=, noz=1, y=1, x=1) ; units min

  time = nday*24. + nhour + nmin/60.

  ntstep = dimsizes(nday)

; import data
; -----------
  ; rho, to convert units
  rho = MeteoFile->rho(:,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units kg/m3
; convert to m3/mol
  M_air=28.96546e-3 ; kg/mol
  xair = M_air / rho

  ; Br2
  Br2 = GasFile->Br2(:,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3

  ; BrCl
  BrCl = GasFile->BrCl(:,gidx,0,0) ; file is (rec=, n=2, y=1, x=1) ; units mol/m3

; convert to pmol/mol
  Br2  = Br2 * xair * 1e12
  BrCl = BrCl * xair * 1e12


;************************************************
; create plots
;************************************************
  wks = gsn_open_wks("pdf","Buys_2013_Fig3")          ; send graphics to PDF file
  plot = new(2,graphic)

  res                              = True
  res@gsnDraw                      = False
  res@gsnFrame                     = False
  res@gsnMaximize                  = True            ; use full page
  res@gsnPaperOrientation          = "portrait"      ; prevent automatic rotation
  ; adjust plot to get same height and width
  res@vpHeightF           = 0.44
  res@vpWidthF            = 0.7

  ; X-Axis settings
  res@trXMinF           = -1.5
  res@trXMaxF           = 27.
  res@tmXBMode          = "Manual"
  res@tmXBTickStartF    = 0.
  res@tmXBTickEndF      = 24.
  res@tmXBTickSpacingF  = 6.
  res@tmXBMinorPerMajor = 5
  res@tmXTOn            = False
  res@tiYAxisString     = "Time (h)"
  res@tiXAxisFontHeightF = 0.02
  ; Y-Axis settings
  ;res@trYMinF           = 0.
  res@trYMaxF           = 15.
  res@tmYLMode          = "Manual"
  res@tmYLTickStartF    = 0.
  res@tmYLTickEndF      = 14.
  res@tmYLTickSpacingF  = 2.
  res@tmYLMinorPerMajor = 1
  res@tmYROn            = False
  res@tiYAxisString     = "Br~B~2~N~ (pmol mol~S~-1~N~)"
  res@tiYAxisFontHeightF = 0.02
  ; Informations displayed on top
  res@gsnLeftString = ""
  res@gsnRightString = ""

  plot(0) = gsn_csm_xy(wks,time,Br2,res)

  ; Adjust plot settings for BrCl
  ; X-Axis settings
  res@trXMinF           = -2.5
  res@trXMaxF           = 26.5
  ; Y-Axis settings
  res@trYMinF           = -0.05
  res@trYMaxF           = 3.
  res@tmYLMode          = "Manual"
  res@tmYLTickStartF    = 0.
  res@tmYLTickEndF      = 3.
  res@tmYLTickSpacingF  = 0.5
  res@tmYLMinorPerMajor = 4
  res@tiYAxisString     = "BrCl (pmol mol~S~-1~N~)"
  plot(1) = gsn_csm_xy(wks,time,BrCl,res)

;************************************************
; create panel
;************************************************
  gsn_panel(wks,plot,(/1,2/),False)             ; now draw as one plot

end
