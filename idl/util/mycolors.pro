function mycolors,load=load,table=table,goodgreen=goodgreen,badgreen=badgreen
; loads color tables and defines indices for colors

if not keyword_set(table) then table=39
if keyword_set(load) then loadct,table ;load rainbow + white color table
ntab=!d.table_size
black=0
brick=250*ntab/256.
red = 220*ntab/256.
blue = 46.*ntab/256.
purple = 26.*ntab/256.
cyan = 80.*ntab/256.
turquoise = 100.*ntab/256.
bluegreen = 130.*ntab/256.
green = 150.*ntab/256.  ;140
yellow = 180.*ntab/256.
orange = 200.*ntab/256.
white=ntab-1
!p.background=white
!p.color=black


; Make green not terrible!
if keyword_set(goodgreen) then begin
    tvlct,r,g,b,/get
    r(green) = 0
    g(green) = 130
    b(green) = 0

    tvlct,r,g,b
 endif

; Make green terrible again!
if keyword_set(badgreen) then begin
    tvlct,r,g,b,/get
    r(green) = 0
    g(green) = 255
    b(green) = 38

    tvlct,r,g,b
endif


a={black:black, red:red, blue:blue, white:white, cyan:cyan, green:green, $
   yellow:yellow, orange:orange, purple:purple, turquoise:turquoise, $
   bluegreen:bluegreen, brick:brick}
return,a
end
