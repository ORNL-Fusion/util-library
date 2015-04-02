pro closeall
while !d.window ne -1 do wdelete, !d.window
end


function sym_flipper

flip = [4,6,7,1,5]
flip = [flip,flip,flip]
flip = [flip,flip,flip]

return,flip
end
function color_flipper,blackfirst=blackfirst

z=mycolors()
;flip = [z.black,z.brick,z.blue,z.green,z.purple,z.cyan,z.orange]
flip = [z.black,z.brick,z.blue,z.cyan,z.orange,z.purple,z.green]

flip=[flip,flip,flip,flip,flip,flip,flip,flip]

return,flip
end

function alphabet
return,['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
end

function plot_positions,nrow=nrow,ncol=ncol,bottom=bottom,left=left,top=top,right=right,$
                        horzspace=horzspace,vertspace=vertspace

if not keyword_set(nrow) then nrow=1
if not keyword_set(ncol) then ncol=1
if not keyword_set(bottom) then bottom = 0.15
if not keyword_set(left) then left = 0.18
if not keyword_set(top) then top=0.99
if not keyword_set(right) then right=0.99
if not keyword_set(horzspace) then horzspace=0.1
if not keyword_set(vertspace) then vertspace=0.01


width = (right-left-horzspace*(ncol-1))/ncol
height = (top-bottom-vertspace*(nrow-1))/nrow
pos = fltarr(4,nrow*ncol)

for j=0,ncol-1 do begin
    for i=0,nrow-1 do begin
        pos(0,j*nrow+i) = left + j*(width+horzspace)
        pos(1,j*nrow+i) = top - height - i*(height+vertspace)
        pos(2,j*nrow+i) = left + width + j*(width+horzspace)
        pos(3,j*nrow+i) = top - i*(height+vertspace)
    endfor
endfor

return,pos
end


pro ylabel,text,pos,offset=offset,charthick=charthick,charsize=charsize

if not keyword_set(offset) then offset=0.

xyouts,pos(0)-0.1+offset,0.5*(pos(1)+pos(3)),textoidl(text),alignment=0.5,orientation=90,$
  charthick=charthick,charsize=charsize,/normal

end


pro xlabel,text,pos,offset=offset,charthick=charthick,charsize=charsize

if not keyword_set(offset) then offset=0.

xyouts,0.5*(pos(0)+pos(2)),pos(1)-0.05+offset,textoidl(text),alignment=0.5, $
  charthick=charthick,charsize=charsize,/normal 

end

pro legend,text,pos,offsetx=offsetx,offsety=offsety,corner=corner,spacing=spacing, $
           charthick=charthick,charsize=charsize,colors=colors

if not keyword_set(offsetx) then offsetx=0.
if not keyword_set(offsety) then offsety=0.
if not keyword_set(spacing) then spacing=0.03
if not keyword_set(corner) then corner = 'ul'
if not keyword_set(colors) then colors = color_flipper()

case corner of
    'ul': begin
        xx=pos(0)+.02+offsetx 
        yy=pos(3)-0.03+offsety 
        align=0.0 
        sign=-1.
    end
    'bl': begin
        xx=pos(0)+.02+offsetx 
        yy=pos(1)+0.01+offsety 
        align=0.0 
        sign=-1.
    end
    'ur': begin
        xx=pos(2)-.02+offsetx 
        yy=pos(3)-0.03+offsety 
        align=1.0 
        sign=1.
    end
    'br': begin
        xx=pos(2)-.02+offsetx 
        yy=pos(1)+0.01+offsety 
        align=1.0 
        sign=1.
    end
    else: print,'Bad choice of corner, pick one of ul,bl,ur,br'
endcase

for i =0,n_elements(text)-1 do $
  xyouts,xx,yy+i*sign*spacing,textoidl(text(i)),alignment=align, $
  charthick=charthick,charsize=charsize,/normal,color=colors(i)

end


pro make_nice_plots,ps=ps

if not keyword_set(ps) then begin
  !p.charsize=1.5
  !p.charthick=1.6
  !x.thick=1.5
  !y.thick=1.5
  !p.thick=2.0
endif else begin
  !p.charsize=1.2
  !p.charthick=2.0
  !x.thick=2.0
  !y.thick=2.0
  !p.thick=3.0
end

end

;Function to generate a linearly spaced double precision vector
function dlinspace, xstart, xend, numel
  if ( numel eq 1 ) then myvec = xstart else begin
  myvec = xstart + dindgen(numel)/(numel-1.D)*(xend-xstart)
  endelse
  return, myvec
end

;Function to generate a logarithmically spaced double precision vector
function dlogspace, xstart, xend, numel
  myvec = 10.D^dlinspace(alog10(xstart),alog10(xend),numel)
  return, myvec
end


