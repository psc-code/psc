; NAME:
;   COLORBAR
;
; PURPOSE:
;       The purpose of this routine is to add a color bar to the current
;       graphics window.
;
; CATEGORY:
;       Graphics, Widgets.
;
; CALLING SEQUENCE:
;       COLORBAR
;
; INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;
;       BOTTOM: The lowest color index of the colors to be loaded in
;                 the bar.
;
;       CHARSIZE: The character size of the color bar annotations. Default is 1.0.
;
;       COLOR:    The color index of the bar outline and characters. Default
;                 is ncolors - 1 + bottom.
;
;       DIVISIONS: The number of divisions to divide the bar into. There will
;                 be (divisions + 1) annotations. The default is 2.
;
;       FORMAT:   The format of the bar annotations. Default is '(F6.2)'.
;
;       MAX:      The maximum data value for the bar annotation. Default is
;                 NCOLORS-1.
;
;       MIN:      The minimum data value for the bar annotation. Default is 0.
;
;       NCOLORS:  This is the number of colors in the color bar.
;
;       POSITION: A four-element array of normalized coordinates in the same
;                 form as the POSITION keyword on a plot. Default is
;                 [0.88, 0.15, 0.95, 0.95] for a vertical bar and
;                 [0.15, 0.88, 0.95, 0.95] for a horizontal bar.
;
;       PSCOLOR:  This keyword is only applied if the output is being sent to
;                 a PostScript file. It indicates that the PostScript device
;                 is configured for color output. If this keyword is set, then
;                 the annotation is drawn in the color specified by the COLOR
;                 keyword. If the keyword is not set, the annotation is drawn
;                 in the color specified by the !P.COLOR system variable
;                 (usually this will be the color black). In general, this
;                 gives better looking output on non-color or gray-scale
;                 printers. If you are not specifically setting the annotation
;                 color (with the COLOR keyword), it will probably
;                 be better NOT to set this keyword either, even if you
;                 are outputting to a color PostScript printer.
;
;       RIGHT:    This puts the labels on the right-hand side of a vertical
;                 color bar. It applies only to vertical color bars.
;
;       TITLE:    This is title for the color bar. The default is to have
;                 no title.
;
;       TOP:      This puts the labels on top of the bar rather than under it.
;                 The keyword only applies if a horizontal color bar is rendered.
;
;       VERTICAL: Setting this keyword give a vertical color bar. The default
;                 is a horizontal color bar.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       Color bar is drawn in the current graphics window.
;
; RESTRICTIONS:
;       The number of colors available on the display device (not the
;       PostScript device) is used unless the NCOLORS keyword is used.
;
; EXAMPLE:
;       To display a horizontal color bar above a contour plot, type:
;
;       LOADCT, 5, NCOLORS=100
;       CONTOUR, DIST(31,41), POSITION=[0.15, 0.15, 0.95, 0.75], $
;          C_COLORS=INDGEN(25)*4, NLEVELS=25
;       COLORBAR, NCOLORS=100
;
; MODIFICATION HISTORY:
;       Written by: David Fanning, 10 JUNE 96.
;       10/27/96: Added the ability to send output to PostScript. DWF
;       11/4/96: Substantially rewritten to go to screen or PostScript
;           file without having to know much about the PostScript device
;           or even what the current graphics device is. DWF
;       1/27/97: Added the RIGHT and TOP keywords. Also modified the
;            way the TITLE keyword works. DWF
;       7/15/97: Fixed a problem some machines have with plots that have
;            no valid data range in them. DWF
;-

PRO COLORBAR, BOTTOM=bottom, CHARSIZE=charsize, COLOR=color, DIVISIONS=divisions, $
   FORMAT=format, POSITION=position, MAX=max, MIN=min, NCOLORS=ncolors, $
   PSCOLOR=pscolor, TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right

   ; Is the PostScript device selected?

postScriptDevice = (!D.NAME EQ 'PS')

  ; Check and define keywords.

IF N_ELEMENTS(ncolors) EQ 0 THEN BEGIN

   ; Most display devices do not use the 256 colors available to
   ; the PostScript device. This presents a problem when writing
   ; general-purpose programs that can be output to the display or
   ; to the PostScript device. This problem is especially bothersome
   ; if you don't specify the number of colors you are using in the
   ; program. One way to work around this problem is to make the
   ; default number of colors the same for the display device and for
   ; the PostScript device. Then, the colors you see in PostScript are
   ; identical to the colors you see on your display. Here is one way to
   ; do it.

   IF postScriptDevice THEN BEGIN
      oldDevice = !D.NAME

         ; What kind of computer are we using? SET_PLOT to appropriate
         ; display device.

      thisOS = !VERSION.OS_FAMILY
      thisOS = STRMID(thisOS, 0, 3)
      thisOS = STRUPCASE(thisOS)
      CASE thisOS of
         'MAC': SET_PLOT, thisOS
         'WIN': SET_PLOT, thisOS
         ELSE: SET_PLOT, 'X'
      ENDCASE

         ; Open a window (to make sure !D.N_COLORS is accurate).

      WINDOW, /FREE, /PIXMAP, XSIZE=10, YSIZE=10
      WDELETE, !D.WINDOW

         ; Here is how many colors we should use.

      ncolors = !D.N_COLORS
      SET_PLOT, oldDevice
    ENDIF ELSE ncolors = !D.N_COLORS
ENDIF
IF N_ELEMENTS(bottom) EQ 0 THEN bottom = 0B
IF N_ELEMENTS(charsize) EQ 0 THEN charsize = 1.0
IF N_ELEMENTS(format) EQ 0 THEN format = '(F6.2)'
IF N_ELEMENTS(color) EQ 0 THEN color = ncolors - 1 + bottom
IF N_ELEMENTS(min) EQ 0 THEN min = 0.0
IF N_ELEMENTS(max) EQ 0 THEN max = FLOAT(ncolors) - 1
IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 2
IF N_ELEMENTS(title) EQ 0 THEN title = ''
pscolor = KEYWORD_SET(pscolor)

IF KEYWORD_SET(vertical) THEN BEGIN
   bar = REPLICATE(1B,10) # BINDGEN(256)
   IF N_ELEMENTS(position) EQ 0 THEN position = [0.88, 0.15, 0.95, 0.95]
ENDIF ELSE BEGIN
   bar = BINDGEN(256) # REPLICATE(1B, 10)
   IF N_ELEMENTS(position) EQ 0 THEN position = [0.15, 0.88, 0.95, 0.95]
ENDELSE

   ; Scale the color bar.

bar = BYTSCL(bar, TOP=ncolors-1) + bottom

   ; Get starting locations in DEVICE coordinates.

xstart = position(0) * !D.X_VSIZE
ystart = position(1) * !D.Y_VSIZE

   ; Get the size of the bar in DEVICE coordinates.

xsize = (position(2) - position(0)) * !D.X_VSIZE
ysize = (position(3) - position(1)) * !D.Y_VSIZE

   ; For PostScript output only, draw the annotation in !P.COLOR
   ; unless "pscolor" is set. This makes better output on grayscale
   ; printers.

IF postScriptDevice AND (pscolor NE 1) THEN BEGIN
   oldcolor = color
   color = !P.COLOR
ENDIF

   ; Display the color bar in the window. Sizing is
   ; different for PostScript and regular display.

IF postScriptDevice THEN BEGIN

   TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize

ENDIF ELSE BEGIN

   bar = CONGRID(bar, CEIL(xsize), CEIL(ysize), /INTERP)
   TV, bar, xstart, ystart

ENDELSE

   ; Annotate the color bar.

IF KEYWORD_SET(vertical) THEN BEGIN

   IF KEYWORD_SET(right) THEN BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', YTICKLEN=0.1 , $
         YRANGE=[min, max], YTITLE=title
         
      AXIS, YAXIS=1, YRANGE=[min, max], YTICKFORMAT=format, YTICKS=divisions, $
         YTICKLEN=0.1, YSTYLE=1, COLOR=color, CHARSIZE=charsize

   ENDIF ELSE BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=1, $
         YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT=format, XTICKFORMAT='(A1)', YTICKLEN=0.1 , $
         YRANGE=[min, max]

      AXIS, YAXIS=1, YRANGE=[min, max], YTICKFORMAT='(A1)', YTICKS=divisions, $
         YTICKLEN=0.1, YTITLE=title, YSTYLE=1, COLOR=color, CHARSIZE=charsize

   ENDELSE

ENDIF ELSE BEGIN

   IF KEYWORD_SET(top) THEN BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=9, YSTYLE=1, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', XTICKLEN=0.1, $
         XRANGE=[min, max], XTITLE=title

      AXIS, XTICKS=divisions, XSTYLE=1, COLOR=color, CHARSIZE=charsize, $
         XTICKFORMAT=format, XTICKLEN=0.1, XRANGE=[min, max], XAXIS=1

   ENDIF ELSE BEGIN

      PLOT, [min,max], [min,max], /NODATA, XTICKS=divisions, $
         YTICKS=1, XSTYLE=1, YSTYLE=1, $
         POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
         YTICKFORMAT='(A1)', XTICKFORMAT=format, XTICKLEN=0.1, $
         XRANGE=[min, max], TITLE=title

    ENDELSE

ENDELSE

   ; Restore color variable if changed for PostScript.

IF postScriptDevice AND (pscolor NE 1) THEN color = oldcolor

END
