;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Collection of routines to perform magnetic field-line
; tracing, similar to TRIP3D (T. Evans, GA)
;
; Started 1/21/09  JC
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_bicub_mat

a = dblarr(16,16)

;;;;;; function values at corners
a(0,0) = 1.0
a(1,[0,1,2,3]) = [1.,1.,1.,1.]
a(2,[0,4,8,12]) = [1.,1.,1.,1.]
a(3,*) = 1.0

;;;;;;; 1st derivatives at corners: x direction
a(4,1)=1.0
a(5,[1,2,3]) = [1.0,2.0,3.0]
a(6,[1,5,9,13]) = [1.0,1.0,1.0,1.0]
a(7,[1,2,3]) = [1.,2.,3.]
a(7,[5,6,7]) = [1.,2.,3.]
a(7,[9,10,11]) = [1.,2.,3.]
a(7,[13,14,15]) = [1.,2.,3.]

;;;;;;; 1st derivatives at corners: y direction
a(8,4) = 1.0
a(9,4:7) = [1.,1.,1.,1.]
a(10,[4,8,12]) = [1.,2.,3.]
a(11,[4,8,12]) = [1.,2.,3.]
a(11,[5,9,13]) = [1.,2.,3.]
a(11,[6,10,14]) = [1.,2.,3.]
a(11,[7,11,15]) = [1.,2.,3.]

;;;;;;; cross derivatives at corners
a(12,5) = 1.
a(13,[5,6,7]) = [1.,2.,3.]
a(14,[5,9,13]) = [1.,2.,3.]
a(15,[5,9,13]) = [1.,2.,3.]
a(15,[6,10,14]) = [2.,4.,6.]
a(15,[7,11,15]) =[3.,6.,9.]

return,a
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_psi_bicub_coeffs,g

nr = g.mw
nz = g.mh

dr = g.r(1)-g.r(0)
dz = g.z(1)-g.z(0)

ip_sign=-g.cpasma/abs(g.cpasma)
psi2d=ip_sign*g.psirz

dsdr = (shift(psi2d,-1,0)-shift(psi2d,1,0))/(2*dr)
dsdz = (shift(psi2d,0,-1)-shift(psi2d,0,1))/(2*dz)
d2sdrdz = (shift(psi2d,-1,-1) - shift(psi2d,1,-1) - shift(psi2d,-1,1) + shift(psi2d,1,1))/(4*dr*dz)

for ir = 1,nr-2 do begin
    dsdz(ir,0) = (psi2d(ir,1) - psi2d(ir,0))/dz
    dsdz(ir,nz-1) = (psi2d(ir,nz-1) - psi2d(ir,nz-2))/dz
    d2sdrdz(ir,0) = (psi2d(ir+1,1)-psi2d(ir+1,0)-psi2d(ir-1,1)+psi2d(ir-1,0))/(2*dr*dz)
    d2sdrdz(ir,nz-1) = (psi2d(ir+1,nz-1)-psi2d(ir+1,nz-2)-psi2d(ir-1,nz-1)+psi2d(ir-1,nz-2))/(2*dr*dz)
endfor

for iz = 1,nz-2 do begin
    dsdr(0,iz) = (psi2d(1,iz)-psi2d(0,iz))/dr
    dsdr(nr-1,iz) = (psi2d(nr-1,iz)-psi2d(nr-2,iz))/dr
    d2sdrdz(0,iz) = (psi2d(1,iz+1)-psi2d(0,iz+1)-psi2d(1,iz-1)+psi2d(0,iz-1))/(2*dr*dz)
    d2sdrdz(nr-1,iz) = (psi2d(nr-1,iz+1)-psi2d(nr-2,iz+1)-psi2d(nr-1,iz-1)+psi2d(nr-2,iz-1))/(2*dr*dz)
endfor

a = get_bicub_mat()

c=dblarr(nr*nz,4,4)
for ir=0,nr-2 do begin
    for iz=0,nz-2 do begin
        index = iz + nz*ir
        b = [ psi2d(ir,iz),        psi2d(ir+1,iz),        psi2d(ir,iz+1),        psi2d(ir+1,iz+1), $
              dsdr(ir,iz)*dr,      dsdr(ir+1,iz)*dr,      dsdr(ir,iz+1)*dr,      dsdr(ir+1,iz+1)*dr, $
              dsdz(ir,iz)*dz,      dsdz(ir+1,iz)*dz,      dsdz(ir,iz+1)*dz,      dsdz(ir+1,iz+1)*dz, $
              d2sdrdz(ir,iz)*dr*dz,d2sdrdz(ir+1,iz)*dr*dz,d2sdrdz(ir,iz+1)*dr*dz,d2sdrdz(ir+1,iz+1)*dr*dz]
        coeff = transpose(invert(a))##b
        c(index,*,0) = coeff(0:3)
        c(index,*,1) = coeff(4:7)
        c(index,*,2) = coeff(8:11)
        c(index,*,3) = coeff(12:15)
    endfor
endfor

return,c
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_bfield_bicub_coeffs,f

tstart = systime(1)

nr = f.nr
nz = f.nz
nt = f.nt

dr = f.dr
dz = f.dz

br2d=f.br
bz2d=f.bz
bphi2d=f.bphi

a = get_bicub_mat()
cbr=dblarr(nr*nz*nt,4,4)
cbz=dblarr(nr*nz*nt,4,4)
cbt=dblarr(nr*nz*nt,4,4)

for it=0,nt-1 do begin
   psi2d = reform(br2d(*,it,*))

   dsdr = (shift(psi2d,-1,0)-shift(psi2d,1,0))/(2*dr)
   dsdz = (shift(psi2d,0,-1)-shift(psi2d,0,1))/(2*dz)
   d2sdrdz = (shift(psi2d,-1,-1) - shift(psi2d,1,-1) - shift(psi2d,-1,1) + shift(psi2d,1,1))/(4*dr*dz)

   for ir = 1,nr-2 do begin
      dsdz(ir,0) = (psi2d(ir,1) - psi2d(ir,0))/dz
      dsdz(ir,nz-1) = (psi2d(ir,nz-1) - psi2d(ir,nz-2))/dz
      d2sdrdz(ir,0) = (psi2d(ir+1,1)-psi2d(ir+1,0)-psi2d(ir-1,1)+psi2d(ir-1,0))/(2*dr*dz)
      d2sdrdz(ir,nz-1) = (psi2d(ir+1,nz-1)-psi2d(ir+1,nz-2)-psi2d(ir-1,nz-1)+psi2d(ir-1,nz-2))/(2*dr*dz)
   endfor

   for iz = 1,nz-2 do begin
      dsdr(0,iz) = (psi2d(1,iz)-psi2d(0,iz))/dr
      dsdr(nr-1,iz) = (psi2d(nr-1,iz)-psi2d(nr-2,iz))/dr
      d2sdrdz(0,iz) = (psi2d(1,iz+1)-psi2d(0,iz+1)-psi2d(1,iz-1)+psi2d(0,iz-1))/(2*dr*dz)
      d2sdrdz(nr-1,iz) = (psi2d(nr-1,iz+1)-psi2d(nr-2,iz+1)-psi2d(nr-1,iz-1)+psi2d(nr-2,iz-1))/(2*dr*dz)
   endfor

   for ir=0,nr-2 do begin
      for iz=0,nz-2 do begin
         index = iz + nz*ir + nr*nz*it
         b = [ psi2d(ir,iz),        psi2d(ir+1,iz),        psi2d(ir,iz+1),        psi2d(ir+1,iz+1), $
               dsdr(ir,iz)*dr,      dsdr(ir+1,iz)*dr,      dsdr(ir,iz+1)*dr,      dsdr(ir+1,iz+1)*dr, $
               dsdz(ir,iz)*dz,      dsdz(ir+1,iz)*dz,      dsdz(ir,iz+1)*dz,      dsdz(ir+1,iz+1)*dz, $
               d2sdrdz(ir,iz)*dr*dz,d2sdrdz(ir+1,iz)*dr*dz,d2sdrdz(ir,iz+1)*dr*dz,d2sdrdz(ir+1,iz+1)*dr*dz]
         coeff = transpose(invert(a))##b
         cbr(index,*,0) = coeff(0:3)
         cbr(index,*,1) = coeff(4:7)
         cbr(index,*,2) = coeff(8:11)
         cbr(index,*,3) = coeff(12:15)
      endfor
   endfor

   psi2d = reform(bz2d(*,it,*))

   dsdr = (shift(psi2d,-1,0)-shift(psi2d,1,0))/(2*dr)
   dsdz = (shift(psi2d,0,-1)-shift(psi2d,0,1))/(2*dz)
   d2sdrdz = (shift(psi2d,-1,-1) - shift(psi2d,1,-1) - shift(psi2d,-1,1) + shift(psi2d,1,1))/(4*dr*dz)

   for ir = 1,nr-2 do begin
      dsdz(ir,0) = (psi2d(ir,1) - psi2d(ir,0))/dz
      dsdz(ir,nz-1) = (psi2d(ir,nz-1) - psi2d(ir,nz-2))/dz
      d2sdrdz(ir,0) = (psi2d(ir+1,1)-psi2d(ir+1,0)-psi2d(ir-1,1)+psi2d(ir-1,0))/(2*dr*dz)
      d2sdrdz(ir,nz-1) = (psi2d(ir+1,nz-1)-psi2d(ir+1,nz-2)-psi2d(ir-1,nz-1)+psi2d(ir-1,nz-2))/(2*dr*dz)
   endfor

   for iz = 1,nz-2 do begin
      dsdr(0,iz) = (psi2d(1,iz)-psi2d(0,iz))/dr
      dsdr(nr-1,iz) = (psi2d(nr-1,iz)-psi2d(nr-2,iz))/dr
      d2sdrdz(0,iz) = (psi2d(1,iz+1)-psi2d(0,iz+1)-psi2d(1,iz-1)+psi2d(0,iz-1))/(2*dr*dz)
      d2sdrdz(nr-1,iz) = (psi2d(nr-1,iz+1)-psi2d(nr-2,iz+1)-psi2d(nr-1,iz-1)+psi2d(nr-2,iz-1))/(2*dr*dz)
   endfor

   for ir=0,nr-2 do begin
      for iz=0,nz-2 do begin
         index = iz + nz*ir + nr*nz*it
         b = [ psi2d(ir,iz),        psi2d(ir+1,iz),        psi2d(ir,iz+1),        psi2d(ir+1,iz+1), $
               dsdr(ir,iz)*dr,      dsdr(ir+1,iz)*dr,      dsdr(ir,iz+1)*dr,      dsdr(ir+1,iz+1)*dr, $
               dsdz(ir,iz)*dz,      dsdz(ir+1,iz)*dz,      dsdz(ir,iz+1)*dz,      dsdz(ir+1,iz+1)*dz, $
               d2sdrdz(ir,iz)*dr*dz,d2sdrdz(ir+1,iz)*dr*dz,d2sdrdz(ir,iz+1)*dr*dz,d2sdrdz(ir+1,iz+1)*dr*dz]
         coeff = transpose(invert(a))##b
         cbz(index,*,0) = coeff(0:3)
         cbz(index,*,1) = coeff(4:7)
         cbz(index,*,2) = coeff(8:11)
         cbz(index,*,3) = coeff(12:15)
      endfor
   endfor

   psi2d = reform(bphi2d(*,it,*))

   dsdr = (shift(psi2d,-1,0)-shift(psi2d,1,0))/(2*dr)
   dsdz = (shift(psi2d,0,-1)-shift(psi2d,0,1))/(2*dz)
   d2sdrdz = (shift(psi2d,-1,-1) - shift(psi2d,1,-1) - shift(psi2d,-1,1) + shift(psi2d,1,1))/(4*dr*dz)

   for ir = 1,nr-2 do begin
      dsdz(ir,0) = (psi2d(ir,1) - psi2d(ir,0))/dz
      dsdz(ir,nz-1) = (psi2d(ir,nz-1) - psi2d(ir,nz-2))/dz
      d2sdrdz(ir,0) = (psi2d(ir+1,1)-psi2d(ir+1,0)-psi2d(ir-1,1)+psi2d(ir-1,0))/(2*dr*dz)
      d2sdrdz(ir,nz-1) = (psi2d(ir+1,nz-1)-psi2d(ir+1,nz-2)-psi2d(ir-1,nz-1)+psi2d(ir-1,nz-2))/(2*dr*dz)
   endfor

   for iz = 1,nz-2 do begin
      dsdr(0,iz) = (psi2d(1,iz)-psi2d(0,iz))/dr
      dsdr(nr-1,iz) = (psi2d(nr-1,iz)-psi2d(nr-2,iz))/dr
      d2sdrdz(0,iz) = (psi2d(1,iz+1)-psi2d(0,iz+1)-psi2d(1,iz-1)+psi2d(0,iz-1))/(2*dr*dz)
      d2sdrdz(nr-1,iz) = (psi2d(nr-1,iz+1)-psi2d(nr-2,iz+1)-psi2d(nr-1,iz-1)+psi2d(nr-2,iz-1))/(2*dr*dz)
   endfor

   for ir=0,nr-2 do begin
      for iz=0,nz-2 do begin
         index = iz + nz*ir + nr*nz*it
         b = [ psi2d(ir,iz),        psi2d(ir+1,iz),        psi2d(ir,iz+1),        psi2d(ir+1,iz+1), $
               dsdr(ir,iz)*dr,      dsdr(ir+1,iz)*dr,      dsdr(ir,iz+1)*dr,      dsdr(ir+1,iz+1)*dr, $
               dsdz(ir,iz)*dz,      dsdz(ir+1,iz)*dz,      dsdz(ir,iz+1)*dz,      dsdz(ir+1,iz+1)*dz, $
               d2sdrdz(ir,iz)*dr*dz,d2sdrdz(ir+1,iz)*dr*dz,d2sdrdz(ir,iz+1)*dr*dz,d2sdrdz(ir+1,iz+1)*dr*dz]
         coeff = transpose(invert(a))##b
         cbt(index,*,0) = coeff(0:3)
         cbt(index,*,1) = coeff(4:7)
         cbt(index,*,2) = coeff(8:11)
         cbt(index,*,3) = coeff(12:15)
      endfor
   endfor
endfor

print,'elapsed time: ',systime(1)-tstart,' sec'

return,{cbr:cbr,cbz:cbz,cbt:cbt}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_psi_bicub_coeffs_old,g
;;; This version just drops the edge grid cells

nr = g.mw
nz = g.mh

dr = g.r(1)-g.r(0)
dz = g.z(1)-g.z(0)

ip_sign=-g.cpasma/abs(g.cpasma)
psi2d=ip_sign*g.psirz

dsdr = (shift(psi2d,-1,0)-shift(psi2d,1,0))/(2*dr)
dsdz = (shift(psi2d,0,-1)-shift(psi2d,0,1))/(2*dz)
d2sdrdz = (shift(psi2d,-1,-1) - shift(psi2d,1,-1) - shift(psi2d,-1,1) + shift(psi2d,1,1))/(4*dr*dz)

a = get_bicub_mat()

c=dblarr(nr*nz,4,4)
for ir=0,nr-2 do begin
    for iz=0,nz-2 do begin
        index = iz + nz*ir
        b = [ psi2d(ir,iz),        psi2d(ir+1,iz),        psi2d(ir,iz+1),        psi2d(ir+1,iz+1), $
              dsdr(ir,iz)*dr,      dsdr(ir+1,iz)*dr,      dsdr(ir,iz+1)*dr,      dsdr(ir+1,iz+1)*dr, $
              dsdz(ir,iz)*dz,      dsdz(ir+1,iz)*dz,      dsdz(ir,iz+1)*dz,      dsdz(ir+1,iz+1)*dz, $
              d2sdrdz(ir,iz)*dr*dz,d2sdrdz(ir+1,iz)*dr*dz,d2sdrdz(ir,iz+1)*dr*dz,d2sdrdz(ir+1,iz+1)*dr*dz]
        coeff = transpose(invert(a))##b
        c(index,*,0) = coeff(0:3)
        c(index,*,1) = coeff(4:7)
        c(index,*,2) = coeff(8:11)
        c(index,*,3) = coeff(12:15)
    endfor
endfor

return,c
end

;Version 25-oct-94
;file readg.pro
;
;Modifications:
;4/19/92: comment out lines using the content of "ecase" because an
;         idl bug introduced in version 2.2.2 causes this string to
;         be blank when read from an unformatted file.  See locations
;         near ^^^^^^^^^
;12/4/92: fix the previous bug by reading ecase as a byte array and then
;         converting it to a string.
;         Modify the open statement for unformatted reads so that this
;         routine will work unchanged on the vax.
;1/13/93: modify to handle filenames on a vax.
;10/25/94: adopt version from Cary Forest.  increases grid max size to 129
;          and adds R,Z tags to geqdsk structure.
;
;         change use of arguments: if shot and time are provided as integers
;         a standard g file name is generated.  If a string filename is
;         provided, the file name is used unchanged.
;         This allows special
;         g file name formats to be used for various applications
;         and full pathnames to be used.
;
;10/6/03: S.A. Sabbagh: Add variables to support reading rotation results
;1/11/05: R. Maingi modified line ~ 400 to trap error in dbase reads (temporary)
;----------------------------------------------------------------------
;Read the content of a G eqdsk file and return it in an idl structure.
;
;Calling format:
; g = readg(shot,time)
; or
; g = readg(filename)
;
;Arguments:
;
; shot: integer value giving the shot number.
; time: integer value giving the time in milliseconds.
; filename: string variable equal to the g eqdsk file name 
;
;Returned values:
; g : A structure called geqdsk containing the entire content of 
;     the G eqdsk file.
;     The structure is defined below.
;
; g.error : indicates whether there was an error in reading the file.
;           true (or odd) indicates that there was an error.
;           false or 0 indicates no error.
;
; g.shot:
; g.time: these contain the shot number and time provided as function
;         arguments. If the function finshes properly, these contain the
;         shot and time from the g file.
;
;----------------------------------------------------------------------
;This routine is based on the VMS file sys$d3efit:rdgeqk.for, version
;of 22-MAR-1991.  
;Differences are:
;1. The read of "header" from the G file is not implemented.
;2. To determine if the limiter and boundary data are present in the
;   file the efit version number is checked for both formatted and
;   unformatted files rather than using the end of file for the formatted
;   case.
;3. Only the default directory is search for the G file (there is no
;   central storage location for eqdsk files assumed).
;
;Here is some documentation from rdgeqk.for:
;
;c----------------------------------------------------------------------
;c--  CASE   : descriptive character strings                          --
;c--  XDIM   : radial dimension of rectangular grid in m              --
;c--  ZDIM   : vertical dimension in m                                --
;c--  RGRID(1):minimum major radius in m of grid                      --
;c--  ZMID   : vertical position of box center in m                   --
;c--  BCENTR : toroidal field in T at RZERO                           --
;c--  RMAXIS : major radius of magnetic axis in m                     --
;c--  ZMAXIS : vertical position of magnetic axis in m                --
;c--  SSIMAG : flux at magnetic axis in v-sec/rad                     --
;c--  SSIBRY : boundary flux                                          --
;c--  CPASMA : plasma current in Amp                                  --
;c--  FPOL   : array of poloidal current function RBt at uniform flux --
;c--           grids                                                  --
;c--  PRES   : pressure array in nt/m2 at flux grids                  --
;c--  FFPRIM : radial derivative of FPOL                              --
;c--  PPRIME : radial derivative of PRES                              --
;c--  PSIRZ  : poloidal fluxes at the (R,Z) grids                     --
;c--  QPSI   : q array at flux grids                                  --
;c--  NBDRY  : number of boundary points                              --
;c--  LIMITR : number of limiter points                               --
;c--  RBDRY, ZBDRY : (R,Z) of boundary in m                           --
;c--  XLIM, YLIM : (R,Z) of limiter                                   --
;c----------------------------------------------------------------------
;
;----------------------------------------------------------------------
;
function readg_g3d,arg1,arg2,noname_struct=noname_struct,noload=noload,nosave=nosave
;
;Define the structure to be returned.  Initially the structure is filled
;with zeros.
;
limmax = 900
bdrymax = 900
;
;Initialize the error flag to say that there was an error.  This will
;be changed later if this routine executes successfully.
;

;
;----------------------------------------------------------------------
;Determine the calling format used.
;
;First, check the second argument. If it is an integer/long then get
;the time value from it and the shot number from the first argument.
;If the second argument is undefined, then
;check the first argument for a string and get the shot and time from the
;string.
;
t2=size(arg2)
t1=size(arg1)
;
if( (t2(1+t2(0)) eq 2) or (t2(1+t2(0)) eq 3) ) then begin 
;
;  arg2 =integer or longword
;
   time = arg2(0)
;
   if( (t1(1+t1(0)) eq 2) or (t1(1+t1(0)) eq 3) ) then begin  
;
;     arg1 = integer or longword
;
      shot = arg1(0)
      filename='g' + strtrim(string(shot,'(i6.6)'),2)+$
               '.'+strtrim(string(time,'(i5.5)'))
;
   endif else begin
;
      print,'If there are two arguments for readg, both must be integer/long'
      return,123
;
   endelse
;
endif else if(t2(1+t2(0)) eq 0) then begin	      
;
;  arg2 is undefined.
;
   if(t1(1+t1(0)) eq 7) then begin		  
;
;     arg1=string
;
; The file name is the same as the argument.  The shot number and time
; are extracted from the argument.
;
      on_ioerror,argproblem

      iposv = strpos(arg1,']')       ;         for  Vaxes
      iposu = strpos(arg1,'/')       ;         for  Unix

      if(iposv(0) ne -1 ) then begin
         istart = ipos(0)+1
      endif else if(iposu(0) ne -1) then begin
        ipos_ = iposu				;modification to
	while  ipos_ ne -1   do begin           ;find the last character
           ipos = ipos_				;CBF
           ipos_ = strpos(arg1,'/',ipos + 1)  
	endwhile       
        istart = ipos(n_elements(ipos)-1)+1
      endif else begin
         istart = 0
      endelse
;
      shot = long(strmid(arg1(0),istart+1,6))
      time = long(strmid(arg1(0),istart+8,5))
      filename = arg1
;
   endif else begin
;
;     arg2 undefined and arg1 is not a string
;
      print,'readg requires 2 integer/long or 1 string argument.'
      return,123
;
   endelse
;
endif else begin
;
;  arg2 is not integer/long and is not undefined.
;
   print,'readg requires 2 integer/long or 1 string argument.'
   return,123
;
endelse
;
goto,keepgo
;
argproblem:
   print,'"'+arg1(0)+'" could not be converted to a shot and time.'
   return,123
;
;----------------------------------------------------------------------
keepgo:
;
;Save the shot number and time.
;


savname = filename + '.sav'
if file_test(savname) and not keyword_set(noload) then begin
    restore,savname
    goto,finishup
endif

;
;----------------------------------------------------------------------
;Open the file for formatted reads.
;
on_ioerror,nofile
openr,lun,filename,/get_lun
fileformatted = 1
;
;Create the variables for the first read.
;
ecase = strarr(6)
idum = long(0)
mw = long(0)
mh = long(0)
;
;Read the first group of variables assuming that the file is formatted.
;If there is an error jump ahead to try to read the file assuming that
;it is unformatted.
;
on_ioerror,itsunformatted
;
readf,lun,ecase,idum,mw,mh,format='(6a8,3i4)'
goto,keepgo2
;
itsunformatted:
;
;Close the file then reopen it for unformatted reads.
;
on_ioerror,miscioerror
close,lun
on_ioerror,nofile
openr,lun,filename,/f77_unformatted,/segmented
;
;Create the variables for the first read.
;
ecase = strarr(6)
ecase(*) = '          '
bdum = bytarr(10,6)
idum = long(0)
mw = long(0)
mh = long(0)
;
on_ioerror,bigerror
;
;readu,lun,ecase,idum,mw,mh
readu,lun,bdum,idum,mw,mh
fileformatted = 0
ecase = string(bdum)
;
keepgo2:
;
;At this point we should have the values of mw and mh.  A typical error
;seems to be that these are 0 for some as yet undetermined reason.
;

if( (mw eq 0) or (mh eq 0) ) then return,g  ;g.error will be true (=1).

if keyword_set(noname_struct) then $
g={shot:long(0),time:long(0),error:long(0),$
          ecase:strarr(6),mw:long(0),mh:long(0),xdim:double(0),$
          zdim:double(0),rzero:double(0),rgrid1:double(0),zmid:double(0),$
          rmaxis:double(0),zmaxis:double(0),ssimag:double(0),ssibry:double(0),$
          bcentr:double(0),cpasma:double(0),$
          fpol:dblarr(mw),pres:dblarr(mw),ffprim:dblarr(mw),$
          pprime:dblarr(mw),psirz:dblarr(mw,mh),qpsi:dblarr(mw),$
          nbdry:long(0),limitr:long(0),bdry:dblarr(2,bdrymax),$
          lim:dblarr(2,limmax),R:dblarr(mw),Z:dblarr(mh),filename:'', $
          bfield_type:'',bicub_coeffs:dblarr(mw*mh,4,4),pn:dblarr(mw),fpol_coeffs:dblarr(8),psitor:dblarr(mw)} $
else $
g={shot:long(0),time:long(0),error:long(0),$
          ecase:strarr(6),mw:long(0),mh:long(0),xdim:double(0),$
          zdim:double(0),rzero:double(0),rgrid1:double(0),zmid:double(0),$
          rmaxis:double(0),zmaxis:double(0),ssimag:double(0),ssibry:double(0),$
          bcentr:double(0),cpasma:double(0),$
          fpol:dblarr(mw),pres:dblarr(mw),ffprim:dblarr(mw),$
          pprime:dblarr(mw),psirz:dblarr(mw,mh),qpsi:dblarr(mw),$
          nbdry:long(0),limitr:long(0),bdry:dblarr(2,bdrymax),$
          lim:dblarr(2,limmax),R:dblarr(mw),Z:dblarr(mh),filename:'',$
          bfield_type:'',bicub_coeffs:dblarr(mw*mh,4,4),pn:dblarr(mw),fpol_coeffs:dblarr(8),psitor:dblarr(mw)}
g.error = 1
g.shot = shot
g.time = time
g.bfield_type = 'geq'
g.filename = filename
;
;----------------------------------------------------------------------
;create the second group of variables.
;
xdim=double(0)
zdim=double(0)
rzero=double(0)
rgrid1=double(0)
zmid = double(0)
rmaxis = double(0)
zmaxis = double(0)
ssimag = double(0)
ssibry = double(0)
bcentr = double(0)   
cpasma = double(0)
xdum = double(0)
xdum1 = double(0)
fpol = dblarr(mw)
pres = dblarr(mw)
ffprim = dblarr(mw)
pprime = dblarr(mw)
psirz = dblarr(mw,mh)
qpsi = dblarr(mw)
;
;Read the second group of variables.
;
if(fileformatted) then begin
   readf,lun,xdim,zdim,rzero,rgrid1,zmid,format='(5e16.9)'
   readf,lun,rmaxis,zmaxis,ssimag,ssibry,bcentr,format='(5e16.9)'
   readf,lun,cpasma,ssimag,xdum,rmaxis,xdum,format='(5e16.9)'
   readf,lun,zmaxis,xdum,ssibry,xdum,xdum,format='(5e16.9)'
   readf,lun,fpol,format='(5e16.9)'
   readf,lun,pres,format='(5e16.9)'
   readf,lun,ffprim,format='(5e16.9)'
   readf,lun,pprime,format='(5e16.9)'
   readf,lun,psirz,format='(5e16.9)'
   readf,lun,qpsi,format='(5e16.9)'

endif else begin
   readu,lun,xdim,zdim,rzero,rgrid1,zmid
   readu,lun,rmaxis,zmaxis,ssimag,ssibry,bcentr
   readu,lun,cpasma,ssimag,xdum,rmaxis,xdum
   readu,lun,zmaxis,xdum,ssibry,xdum,xdum
   readu,lun,fpol
   readu,lun,pres
   readu,lun,ffprim
   readu,lun,pprime
   readu,lun,psirz
   readu,lun,qpsi
endelse
;
;----------------------------------------------------------------------
;Find the efit version number that was used to create the file.
;If it is large enough, then boundary and limiter data is present
;in the file so it can be read.
;
;Create the third group of variables.  These are needed even if the
;data cannot be read from the file.
;
nbdry = long(0)
limitr = long(0)
;
;nvernum = long(strmid(ecase(3-1),2-1,2)+strmid(ecase(2-1),4-1,2)+$
;               strmid(ecase(2-1),7-1,2))
;
;if(nvernum ge 870520) then begin
;
;Read the third group of variables.
;
   if(fileformatted) then begin
      readf,lun,nbdry,limitr,format='(2i5)'
   endif else begin
      readu,lun,nbdry,limitr
   endelse
;
;Make the fourth group of variables.
;
   if nbdry ne 0 then   bdry = dblarr(2,nbdry)
   lim = dblarr(2,limitr)
;
;Read the fourth group of variables.
;
 chardm = ''
   if(fileformatted) then begin
      if nbdry ne 0 then $
         readf,lun,bdry,format='(5e16.9)' else  $
         readf,lun, chardm
      readf,lun,lim,format='(5e16.9)'
   endif else begin
      readu,lun,bdry
      readu,lun,lim
   endelse
;
;endif
;
;----------------------------------------------------------------------
;NOTE: at this point the header information is not implemented on 
;non-VAX computers.
;
;header = string(replicate(32b,42))
;print,'header ="'+header+'"'
;;
;if(fileformatted) then begin
;   readf,lun,header,format='(a43)'
;endif else begin
;   readu,lun,header
;endelse
;;
;print,'header ="'+header+'"'
;
;----------------------------------------------------------------------
;Copy the variables into the structure.
;
g.ecase = ecase
g.mw = mw
g.mh = mh
g.xdim=xdim
g.zdim=zdim
g.rzero=rzero
g.rgrid1=rgrid1
g.zmid=zmid
g.rmaxis=rmaxis
g.zmaxis=zmaxis
g.ssimag=ssimag 
g.ssibry=ssibry 
g.bcentr=bcentr 
g.cpasma=cpasma 
g.fpol=fpol
g.pres=pres 
g.ffprim=ffprim 
g.pprime=pprime 
for i=0,mh-1 do g.psirz(0:mw-1,i) = psirz(0:mw-1,i)
g.qpsi=qpsi 
g.nbdry=nbdry 
g.limitr=limitr 
if(nbdry gt 0) then begin
   if(nbdry le bdrymax) then for i=0,nbdry-1 do $
                                    g.bdry(0:1,i)=bdry(0:1,i)
endif
if(limitr gt 0) then begin
   if(limitr le limmax) then for i=0,limitr-1 do $
                                    g.lim(0:1,i)=lim(0:1,i)
endif

;************** added g.r,g.z calc and changed signs on flux--B. Rice
dR = xdim/(mw-1)
dz = zdim/(mh-1)

for i=0,mw-1 do begin
   g.r(i) = rgrid1+i*dR
endfor

for i=0,mh-1 do begin
    g.z(i) = zmid-0.5*zdim+i*dz
endfor

g.pn = dindgen(g.mw)/(g.mw-1)
g.fpol_coeffs = poly_fit(g.pn,g.fpol,7)

g.bicub_coeffs = get_psi_bicub_coeffs(g)
dpsi = (g.pn(1)-g.pn(0))*(g.ssibry-g.ssimag)
g.psitor = [0,total(0.5*(g.qpsi(0:g.mw-2)+g.qpsi(1:g.mw-1))*dpsi,/cumulative)]
g.psitor = (g.psitor - g.psitor(0))/(g.psitor(g.mw-1) - g.psitor(0))


;  ip_sign= g.cpasma/abs(g.cpasma)
; g.ssimag= -g.ssimag*ip_sign ;change signs here so Bz comes out with right sign
;  g.ssibry= -g.ssibry*ip_sign
;  g.psirz= -g.psirz*ip_sign 

;************************
;
;----------------------------------------------------------------------
;Get the actual shot number and time from the G file data.
;
;

;
;----------------------------------------------------------------------
;The routine executed successfully if the boundary or limiter arrays
;did not overflow.
;
if( (nbdry le bdrymax) and (limitr le limmax) ) then g.error = 0
;
;----------------------------------------------------------------------
;close the file
;
alldone:
   on_ioerror,miscioerror
   close,lun
   on_ioerror,miscioerror
   free_lun,lun
   if not keyword_set(nosave) then save,filename = savname
;
finishup:

return,g
;
;----------------------------------------------------------------------
;errors
;
nofile:
   print,'In readg, the file "'+filename+'" could not be located.'
   goto,finishup
bigerror:
   print,'Error reading the file "'+filename+'" in readg.'
   goto,alldone
miscioerror:
   print,'I/O error in readg'
   goto,finishup
;
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_psi_bicub,g,r1,z1,derivs=derivs,secderiv=secderiv
dr = g.r(1)-g.r(0)
dz = g.z(1)-g.z(0)

ir = floor((r1-g.r(0))/dr)
iz = floor((z1-g.z(0))/dz)
index = iz + g.mh*ir

dir = (r1 - g.r(ir))/dr
diz = (z1 - g.z(iz))/dz

c = g.bicub_coeffs(index,*,*)

dir2 = dir^2
dir3 = dir^3
diz2 = diz^2
diz3 = diz^3

psi1 = c(*,0,0)       + c(*,1,0)*dir       + c(*,2,0)*dir2       + c(*,3,0)*dir3 + $
  c(*,0,1)*diz   + c(*,1,1)*dir*diz   + c(*,2,1)*dir2*diz   + c(*,3,1)*dir3*diz + $
  c(*,0,2)*diz2 + c(*,1,2)*dir*diz2 + c(*,2,2)*dir2*diz2 + c(*,3,2)*dir3*diz2 + $
  c(*,0,3)*diz3 + c(*,1,3)*dir*diz3 + c(*,2,3)*dir2*diz3 + c(*,3,3)*dir3*diz3

if keyword_set(derivs) or keyword_set(secderiv) then begin
    dsdr1 = c(*,1,0)       + 2*c(*,2,0)*dir       + 3*c(*,3,0)*dir2 + $
      c(*,1,1)*diz   + 2*c(*,2,1)*dir*diz   + 3*c(*,3,1)*dir2*diz + $
      c(*,1,2)*diz2 + 2*c(*,2,2)*dir*diz2 + 3*c(*,3,2)*dir2*diz2 + $
      c(*,1,3)*diz3 + 2*c(*,2,3)*dir*diz3 + 3*c(*,3,3)*dir2*diz3

    dsdz1 = c(*,0,1)         + c(*,1,1)*dir         + c(*,2,1)*dir2         + c(*,3,1)*dir3 + $
      2*c(*,0,2)*diz   + 2*c(*,1,2)*dir*diz   + 2*c(*,2,2)*dir2*diz   + 2*c(*,3,2)*dir3*diz + $
      3*c(*,0,3)*diz2 + 3*c(*,1,3)*dir*diz2 + 3*c(*,2,3)*dir2*diz2 + 3*c(*,3,3)*dir3*diz2
    dsdr1 = dsdr1/dr
    dsdz1 = dsdz1/dz
    if keyword_set(secderiv) then begin
           d2sdr21 =  2*c(*,2,0)       + 6*c(*,3,0)*dir + $
                      2*c(*,2,1)*diz   + 6*c(*,3,1)*dir*diz + $
                      2*c(*,2,2)*diz2 + 6*c(*,3,2)*dir*diz2 + $
                      2*c(*,2,3)*diz3 + 6*c(*,3,3)*dir*diz3

           d2sdz21 =  2*c(*,0,2)   + 2*c(*,1,2)*dir   + 2*c(*,2,2)*dir2   + 2*c(*,3,2)*dir3 + $
                      6*c(*,0,3)*diz + 6*c(*,1,3)*dir*diz + 6*c(*,2,3)*dir2*diz + 6*c(*,3,3)*dir3*diz
           
           d2sdrdz1 =  c(*,1,1)       + 2*c(*,2,1)*dir       + 3*c(*,3,1)*dir2 + $
                       2*c(*,1,2)*diz   + 4*c(*,2,2)*dir*diz   + 6*c(*,3,2)*dir2*diz + $
                       3*c(*,1,3)*diz2 + 6*c(*,2,3)*dir*diz2 + 9*c(*,3,3)*dir2*diz^2
    
           d2sdr21 = d2sdr21/(dr*dr)
           d2sdz21 = d2sdz21/(dz*dz)
           d2sdrdz1 = d2sdrdz1/(dr*dz)
           return,{psi:psi1,dsdr:dsdr1,dsdz:dsdz1,d2sdr2:d2sdr21,d2sdz2:d2sdz21,d2sdrdz:d2sdrdz1}
           endif else return,{psi:psi1,dsdr:dsdr1,dsdz:dsdz1}
endif else return,psi1

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_geq_bicub,g,r1,z1,derivs=derivs,zero=zero

if keyword_set(zero) then begin
    return,{br:0*r1,bz:0*r1,bphi:0*r1,dbrdr:0*r1,dbrdz:0*r1,dbzdr:0*r1,dbzdz:0*r1, $
            dbphidr:0*r1,dbphidz:0*r1}
endif else begin

dr = g.r(1)-g.r(0)
dz = g.z(1)-g.z(0)

ir = floor((r1-g.r(0))/dr)
iz = floor((z1-g.z(0))/dz)

;print,r1,z1,ir,iz

index = iz + g.mh*ir

dir = (r1 - g.r(ir))/dr
diz = (z1 - g.z(iz))/dz

c = g.bicub_coeffs(index,*,*)

dir2 = dir^2
dir3 = dir^3
diz2 = diz^2
diz3 = diz^3

psi1 = c(*,0,0)       + c(*,1,0)*dir       + c(*,2,0)*dir2       + c(*,3,0)*dir3 + $
  c(*,0,1)*diz   + c(*,1,1)*dir*diz   + c(*,2,1)*dir2*diz   + c(*,3,1)*dir3*diz + $
  c(*,0,2)*diz2 + c(*,1,2)*dir*diz2 + c(*,2,2)*dir2*diz2 + c(*,3,2)*dir3*diz2 + $
  c(*,0,3)*diz3 + c(*,1,3)*dir*diz3 + c(*,2,3)*dir2*diz3 + c(*,3,3)*dir3*diz3

dsdr1 = c(*,1,0)       + 2*c(*,2,0)*dir       + 3*c(*,3,0)*dir2 + $
  c(*,1,1)*diz   + 2*c(*,2,1)*dir*diz   + 3*c(*,3,1)*dir2*diz + $
  c(*,1,2)*diz2 + 2*c(*,2,2)*dir*diz2 + 3*c(*,3,2)*dir2*diz2 + $
  c(*,1,3)*diz3 + 2*c(*,2,3)*dir*diz3 + 3*c(*,3,3)*dir2*diz3

dsdz1 = c(*,0,1)         + c(*,1,1)*dir         + c(*,2,1)*dir2         + c(*,3,1)*dir3 + $
  2*c(*,0,2)*diz   + 2*c(*,1,2)*dir*diz   + 2*c(*,2,2)*dir2*diz   + 2*c(*,3,2)*dir3*diz + $
  3*c(*,0,3)*diz2 + 3*c(*,1,3)*dir*diz2 + 3*c(*,2,3)*dir2*diz2 + 3*c(*,3,3)*dir3*diz2

dsdr1 = dsdr1/dr
dsdz1 = dsdz1/dz

br1 = -dsdz1/r1
bz1 = dsdr1/r1

psiN = (-psi1*g.cpasma/abs(g.cpasma) - g.ssimag)/(g.ssibry-g.ssimag)
fpol = poly(psiN,g.fpol_coeffs)
bt1 = g.bcentr*g.rzero/r1
ic = where(psiN le 1.0)
if ic(0) ne -1 then bt1(ic) = fpol(ic)/r1(ic)

iroff =  where(ir lt 0 or ir gt g.mw-2,count)
if count gt 0 then begin
    indroff = iz(iroff) + g.mh*ir(iroff)
    br1(indroff) = 0.0d
    bz1(indroff) = 0.0d         ;
    bt1(indroff) = 1d10
endif

izoff =  where(iz le 1 or iz ge g.mh-2,count)
if count gt 0 then begin
    indzoff = iz(izoff) + g.mh*ir(izoff)
    br1(indzoff) = 0.0d
    bz1(indzoff) = 0.0d         ;
    bt1(indzoff) = 1d10
endif

if keyword_set(derivs) then begin
    d2sdr21 =  2*c(*,2,0)       + 6*c(*,3,0)*dir + $
      2*c(*,2,1)*diz   + 6*c(*,3,1)*dir*diz + $
      2*c(*,2,2)*diz2 + 6*c(*,3,2)*dir*diz2 + $
      2*c(*,2,3)*diz3 + 6*c(*,3,3)*dir*diz3

    d2sdz21 =  2*c(*,0,2)   + 2*c(*,1,2)*dir   + 2*c(*,2,2)*dir2   + 2*c(*,3,2)*dir3 + $
      6*c(*,0,3)*diz + 6*c(*,1,3)*dir*diz + 6*c(*,2,3)*dir2*diz + 6*c(*,3,3)*dir3*diz

    d2sdrdz1 =  c(*,1,1)       + 2*c(*,2,1)*dir       + 3*c(*,3,1)*dir2 + $
      2*c(*,1,2)*diz   + 4*c(*,2,2)*dir*diz   + 6*c(*,3,2)*dir2*diz + $
      3*c(*,1,3)*diz2 + 6*c(*,2,3)*dir*diz2 + 9*c(*,3,3)*dir2*diz^2
    
    d2sdr21 = d2sdr21/(dr*dr)
    d2sdz21 = d2sdz21/(dz*dz)
    d2sdrdz1 = d2sdrdz1/(dr*dz)

    dbrdr1 = -br1/r1 - d2sdrdz1/r1
    dbrdz1 = -d2sdz21/r1
    dbzdr1 = -bz1/r1 + d2sdr21/r1
    dbzdz1 = d2sdrdz1/r1

    dfpol = polyder(psiN,g.fpol_coeffs)
    dfpdr = dfpol*(-dsdr1*g.cpasma/abs(g.cpasma))/(g.ssibry-g.ssimag)
    dfpdz = dfpol*(-dsdz1*g.cpasma/abs(g.cpasma))/(g.ssibry-g.ssimag)
    dbtdr1 = -bt1/r1
    dbtdz1 = 0.0*dbtdr1
    if ic(0) ne -1 then dbtdr1(ic) = dbtdr1(ic) + dfpdr(ic)/r1(ic)
    if ic(0) ne -1 then dbtdz1(ic) = dfpdz(ic)/r1(ic)

    return,{br:br1,bz:bz1,bphi:bt1,dbrdr:dbrdr1,dbrdz:dbrdz1,dbzdr:dbzdr1,dbzdz:dbzdz1, $
            dbphidr:dbtdr1,dbphidz:dbtdz1}
endif else return,{br:br1,bz:bz1,bphi:bt1}

endelse
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_geq_cartes,g,x,y,z,derivs=derivs,zero=zero
x=double(x) & y=double(y) & z=double(z)
r=sqrt(x^2+y^2)
phi=atan(y,x)
if keyword_set(zero) then begin
    return,{bx:0.*x,by:0.*x,bz:0.*x,dbxdx:0.*x,dbxdy:0.*x,dbxdz:0.*x, $
            dbydx:0.*x,dbydy:0.*x,dbydz:0.*x, $
            dbzdx:0.*x,dbzdy:0.*x,dbzdz:0.*x}
endif else begin

    b=bfield_geq_bicub(g,r,z,derivs=derivs)
    bx = b.br*cos(phi) - b.bphi*sin(phi)
    by = b.br*sin(phi) + b.bphi*cos(phi)
    bz = b.bz

    if keyword_set(derivs) then begin
        drdx = x/r
        drdy = y/r
        dphidx = -y/(x^2 + y^2)
        dphidy = x/(x^2 + y^2)
        dbxdr = b.dbrdr*cos(phi) - b.dbphidr*sin(phi)
        dbydr = b.dbrdr*sin(phi) + b.dbphidr*cos(phi)
        dbxdphi = -b.br*sin(phi) - b.bphi*cos(phi)
        dbydphi = b.br*cos(phi) - b.bphi*sin(phi)
        dbxdx = dbxdr*drdx + dbxdphi*dphidx
        dbxdy = dbxdr*drdy + dbxdphi*dphidy
        dbydx = dbydr*drdx + dbydphi*dphidx
        dbydy = dbydr*drdy + dbydphi*dphidy
        dbxdz = b.dbrdz*cos(phi) - b.dbphidz*sin(phi)
        dbydz = b.dbrdz*sin(phi) + b.dbphidz*cos(phi)
        dbzdz = b.dbzdz
        dbzdx = b.dbzdr*drdx
        dbzdy = b.dbzdr*drdy
        return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz, $
                dbydx:dbydx,dbydy:dbydy,dbydz:dbydz, $
                dbzdx:dbzdx,dbzdy:dbzdy,dbzdz:dbzdz}
    endif else return,{bx:bx,by:by,bz:bz}
endelse
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_bs,x,y,z,coil,current,derivs=derivs
npts = size(coil,/dimension)
npts = npts(0)

x=double(x)
y=double(y)
z=double(z)
nobs = n_elements(x)
bx=dblarr(nobs) & by=bx & bz=bx
dbxdx=bx & dbydy=bx & dbxdy=bx & dbxdz=bx & dbydz=bx

c = current(0:npts-2)*1d-7
for i = 0L,nobs-1 do begin
bl = -coil(0:npts-2,0)+x(i)
bm = -coil(0:npts-2,1)+y(i)
bn = -coil(0:npts-2,2)+z(i)
bb = 1/sqrt(bl^2 + bm^2 + bn^2)

cl = -coil(1:npts-1,0)+x(i)
cm = -coil(1:npts-1,1)+y(i)
cn = -coil(1:npts-1,2)+z(i)
cc = 1/sqrt(cl^2 + cm^2 + cn^2)

al = bl-cl
am = bm-cm
an = bn-cn

adotb = (al*bl+am*bm+an*bn)     ;
adotc = (al*cl+am*cm+an*cn)     ;

ul = cm*bn-cn*bm                ;
um = cn*bl-cl*bn                ;
un = cl*bm-cm*bl                ;

recu = 1/(ul^2 + um^2 + un^2)   ;

w = (adotc*cc-adotb*bb)*recu    ;

bx(i) = total(c*w*ul)              ;
by(i) = total(c*w*um)              ;
bz(i) = total(c*w*un)              ;

if keyword_set(derivs) then begin
    adotbb3 = adotb*bb^3        ;
    adotcc3 = adotc*cc^3        ;
    asq = (al*al+am*am+an*an)   ;

    dbm = adotbb3*bm - am*bb + am*cc - adotcc3*cm ;
    dbn = adotbb3*bn - an*bb + an*cc - adotcc3*cn ;

    dcl2 = 2*(al*adotc-cl*asq)*w ;
    dcm2 = 2*(am*adotc-cm*asq)*w ;

    dbxdx(i) = total(c*(ul*(adotbb3*bl - al*bb + al*cc - adotcc3*cl + dcl2)*recu)) ;
    dbxdy(i) = total(c*((ul*dbm + dcl2*um)*recu - an*w)) ;
    dbydy(i) = total(c*(um*(dbm + dcm2)*recu)) ;
    dbxdz(i) = total(c*((ul*dbn + dcl2*un)*recu + am*w)) ;
    dbydz(i) = total(c*((um*dbn + dcm2*un)*recu - al*w)) ;
endif
endfor
if keyword_set(derivs) then begin
    return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz, $
            dbydx:dbxdy,dbydy:dbydy,dbydz:dbydz, $
            dbzdx:dbxdz,dbzdy:dbydz,dbzdz:-(dbxdx+dbydy)}
endif else return,{bx:bx,by:by,bz:bz}

return,0
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_bs_vec,x,y,z,coil,current,derivs=derivs
npts = size(coil,/dimension)
npts = npts(0)

x=double(x)
y=double(y)
z=double(z)
nobs = n_elements(x)
bx=dblarr(nobs) & by=bx & bz=bx
dbxdx=bx & dbydy=bx & dbxdy=bx & dbxdz=bx & dbydz=bx

;c = current(0:npts-2)*1d-7
bl = dblarr((npts-1),nobs)
bm=bl & bn=bl & cl=bl & cm=bl & cn=bl & c=bl

;blnew = cmreplicate(-coil(0:npts-2),nobs)

for i = 0L,nobs-1 do begin
    c(*,i) = current(0:npts-2)*1d-7

    bl(*,i) = -coil(0:npts-2,0)+x(i)
    bm(*,i) = -coil(0:npts-2,1)+y(i)
    bn(*,i) = -coil(0:npts-2,2)+z(i)

    cl(*,i) = -coil(1:npts-1,0)+x(i)
    cm(*,i) = -coil(1:npts-1,1)+y(i)
    cn(*,i) = -coil(1:npts-1,2)+z(i)

endfor
bb = 1/sqrt(bl^2 + bm^2 + bn^2)
cc = 1/sqrt(cl^2 + cm^2 + cn^2)
al = bl-cl
am = bm-cm
an = bn-cn

adotb = (al*bl+am*bm+an*bn)     ;
adotc = (al*cl+am*cm+an*cn)     ;

ul = cm*bn-cn*bm                ;
um = cn*bl-cl*bn                ;
un = cl*bm-cm*bl                ;

recu = 1/(ul^2 + um^2 + un^2)   ;

w = (adotc*cc-adotb*bb)*recu    ;

bx = total(c*w*ul,1)            ;
by = total(c*w*um,1)            ;
bz = total(c*w*un,1)            ;

if keyword_set(derivs) then begin
    adotbb3 = adotb*bb^3        ;
    adotcc3 = adotc*cc^3        ;
    asq = (al*al+am*am+an*an)   ;

    dbm = adotbb3*bm - am*bb + am*cc - adotcc3*cm ;
    dbn = adotbb3*bn - an*bb + an*cc - adotcc3*cn ;

    dcl2 = 2*(al*adotc-cl*asq)*w ;
    dcm2 = 2*(am*adotc-cm*asq)*w ;

    dbxdx(i) = total(c*(ul*(adotbb3*bl - al*bb + al*cc - adotcc3*cl + dcl2)*recu)) ;
    dbxdy(i) = total(c*((ul*dbm + dcl2*um)*recu - an*w)) ;
    dbydy(i) = total(c*(um*(dbm + dcm2)*recu)) ;
    dbxdz(i) = total(c*((ul*dbn + dcl2*un)*recu + am*w)) ;
    dbydz(i) = total(c*((um*dbn + dcm2*un)*recu - al*w)) ;
endif

if keyword_set(derivs) then begin
    return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz, $
            dbydx:dbxdy,dbydy:dbydy,dbydz:dbydz, $
            dbzdx:dbxdz,dbzdy:dbydz,dbzdz:-(dbxdx+dbydy)}
endif else return,{bx:bx,by:by,bz:bz}

return,0
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_bs_cyl,r,phi,z,coil,current,derivs=derivs

x=r*cos(phi)
y=r*sin(phi)
b=bfield_bs(x,y,z,coil,current,derivs=derivs)

br = b.bx*cos(phi) + b.by*sin(phi)
bphi = -b.bx*sin(phi) + b.by*cos(phi)
bz = b.bz

if keyword_set(derivs) then $
  return,{br:br,bphi:bphi,bz:bz,$
          dbxdx:b.dbxdx,dbxdy:b.dbxdy,dbxdz:b.dbxdz, $
          dbydx:b.dbxdy,dbydy:b.dbydy,dbydz:b.dbydz, $
          dbzdx:b.dbxdz,dbzdy:b.dbydz,dbzdz:-(b.dbxdx+b.dbydy)} $
  else return,{br:br,bphi:bphi,bz:bz}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function convert_m3dc1,a

nr = a.nr
nz = a.nz

dr = double(a.dr)
dz = double(a.dz)

r=double(a.r)
z=double(a.z)

br_amp = double(reform(a.br_amp,nr,nz))
br_phase = double(reform(a.br_phase,nr,nz))
bz_amp = double(reform(a.bz_amp,nr,nz))
bz_phase = double(reform(a.bz_phase,nr,nz))
bt_amp = double(reform(a.bt_amp,nr,nz))
bt_phase = double(reform(a.bt_phase,nr,nz))

return,{br_amp:br_amp,bt_amp:bt_amp,bz_amp:bz_amp,br_phase:br_phase,bt_phase:bt_phase,bz_phase:bz_phase, $
        nn:a.nn,dr:dr,dz:dz,nr:nr,nz:nz,r:r,z:z}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_m3dc1,r,phi,z,b,derivs=derivs,scale=scale,efit=efit

if not keyword_set(scale) then scale=13.37d;4.0d

;ir = double((r-b.r(0)))/b.dr
;iz = double((z-b.z(0)))/b.dz

;print,double(r-b.r(0)),b.dr

ir = interpol(dindgen(b.nr),b.r,r)
iz = interpol(dindgen(b.nz),b.z,z)

;br_amp = interpolate(reform(b.br_amp,b.nr,b.nz),ir,iz)
;br_phase = interpolate(reform(b.br_phase,b.nr,b.nz),ir,iz)

;bz_amp = interpolate(reform(b.bz_amp,b.nr,b.nz),ir,iz)
;bz_phase = interpolate(reform(b.bz_phase,b.nr,b.nz),ir,iz)

;bt_amp = interpolate(reform(b.bt_amp,b.nr,b.nz),ir,iz)
;bt_phase = interpolate(reform(b.bt_phase,b.nr,b.nz),ir,iz)


tr = ir-floor(ir)
tz = iz-floor(iz)
ir = floor(ir)
iz = floor(iz)

;bra1 = b.br_amp(ir+iz*b.nr)

;bra = reform(b.br_amp,b.nr,b.nz)
;bza = reform(b.bz_amp,b.nr,b.nz)
;bta = reform(b.bt_amp,b.nr,b.nz)
;brp = reform(b.br_phase,b.nr,b.nz)
;bzp = reform(b.bz_phase,b.nr,b.nz)
;btp = reform(b.bt_phase,b.nr,b.nz)
 
;bra12 = bra(ir,iz)
;print,bra1,bra12
 
;br_amp2 = bra(ir,iz)*(1-tr)*(1-tz) + bra(ir+1,iz)*tr*(1-tz) + bra(ir,iz+1)*(1-tr)*tz + bra(ir+1,iz+1)*tr*tz
;bz_amp2 = bza(ir,iz)*(1-tr)*(1-tz) + bza(ir+1,iz)*tr*(1-tz) + bza(ir,iz+1)*(1-tr)*tz + bza(ir+1,iz+1)*tr*tz
;bt_amp2 = bta(ir,iz)*(1-tr)*(1-tz) + bta(ir+1,iz)*tr*(1-tz) + bta(ir,iz+1)*(1-tr)*tz + bta(ir+1,iz+1)*tr*tz
;br_phase2 = brp(ir,iz)*(1-tr)*(1-tz) + brp(ir+1,iz)*tr*(1-tz) + brp(ir,iz+1)*(1-tr)*tz + brp(ir+1,iz+1)*tr*tz
;bz_phase2 = bzp(ir,iz)*(1-tr)*(1-tz) + bzp(ir+1,iz)*tr*(1-tz) + bzp(ir,iz+1)*(1-tr)*tz + bzp(ir+1,iz+1)*tr*tz
;bt_phase2 = btp(ir,iz)*(1-tr)*(1-tz) + btp(ir+1,iz)*tr*(1-tz) + btp(ir,iz+1)*(1-tr)*tz + btp(ir+1,iz+1)*tr*tz

br_amp = b.br_amp(ir+iz*b.nr)*(1-tr)*(1-tz) + b.br_amp(ir+1+iz*b.nr)*tr*(1-tz) + b.br_amp(ir+(iz+1)*b.nr)*(1-tr)*tz + b.br_amp(ir+1+(iz+1)*b.nr)*tr*tz
bz_amp = b.bz_amp(ir+iz*b.nr)*(1-tr)*(1-tz) + b.bz_amp(ir+1+iz*b.nr)*tr*(1-tz) + b.bz_amp(ir+(iz+1)*b.nr)*(1-tr)*tz + b.bz_amp(ir+1+(iz+1)*b.nr)*tr*tz
bt_amp = b.bt_amp(ir+iz*b.nr)*(1-tr)*(1-tz) + b.bt_amp(ir+1+iz*b.nr)*tr*(1-tz) + b.bt_amp(ir+(iz+1)*b.nr)*(1-tr)*tz + b.bt_amp(ir+1+(iz+1)*b.nr)*tr*tz
br_phase = b.br_phase(ir+iz*b.nr)*(1-tr)*(1-tz) + b.br_phase(ir+1+iz*b.nr)*tr*(1-tz) + b.br_phase(ir+(iz+1)*b.nr)*(1-tr)*tz + b.br_phase(ir+1+(iz+1)*b.nr)*tr*tz
bz_phase = b.bz_phase(ir+iz*b.nr)*(1-tr)*(1-tz) + b.bz_phase(ir+1+iz*b.nr)*tr*(1-tz) + b.bz_phase(ir+(iz+1)*b.nr)*(1-tr)*tz + b.bz_phase(ir+1+(iz+1)*b.nr)*tr*tz
bt_phase = b.bt_phase(ir+iz*b.nr)*(1-tr)*(1-tz) + b.bt_phase(ir+1+iz*b.nr)*tr*(1-tz) + b.bt_phase(ir+(iz+1)*b.nr)*(1-tr)*tz + b.bt_phase(ir+1+(iz+1)*b.nr)*tr*tz

;print,br_amp,br_amp2,bz_amp,bz_amp2,bt_amp,bt_amp2
;print,br_phase,br_phase2,bz_phase,bz_phase2,bt_phase,bt_phase2


;print,ir,ir2,iz,iz2


;br_amp = interpolate(b.br_amp,ir,iz)
;br_phase = interpolate(b.br_phase,ir,iz)

;bz_amp = interpolate(b.bz_amp,ir,iz)
;bz_phase = interpolate(b.bz_phase,ir,iz)

;bt_amp = interpolate(b.bt_amp,ir,iz)
;bt_phase = interpolate(b.bt_phase,ir,iz)

;br_amp = 0
;br_phase = 0

;bz_amp = 0
;bz_phase = 0

;bt_amp = 0
;bt_phase = 0

br = scale*br_amp*sin(br_phase+phi*b.nn)
bz = scale*bz_amp*sin(bz_phase+phi*b.nn)
bt = scale*bt_amp*sin(bt_phase+phi*b.nn)

if keyword_set(efit) then begin
   br_eq = b.br_eq(ir+iz*b.nr)*(1-tr)*(1-tz) + b.br_eq(ir+1+iz*b.nr)*tr*(1-tz) + b.br_eq(ir+(iz+1)*b.nr)*(1-tr)*tz + b.br_eq(ir+1+(iz+1)*b.nr)*tr*tz
   bz_eq = b.bz_eq(ir+iz*b.nr)*(1-tr)*(1-tz) + b.bz_eq(ir+1+iz*b.nr)*tr*(1-tz) + b.bz_eq(ir+(iz+1)*b.nr)*(1-tr)*tz + b.bz_eq(ir+1+(iz+1)*b.nr)*tr*tz
   bt_eq = b.bt_eq(ir+iz*b.nr)*(1-tr)*(1-tz) + b.bt_eq(ir+1+iz*b.nr)*tr*(1-tz) + b.bt_eq(ir+(iz+1)*b.nr)*(1-tr)*tz + b.bt_eq(ir+1+(iz+1)*b.nr)*tr*tz
;   print,br_eq,bz_eq,bt_eq
   br += br_eq
   bz += bz_eq
   bt += bt_eq
endif

return,{br:br,bphi:bt,bz:bz}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function read_extender_file,file=file

if not keyword_set(file) then file = 'extender_points.out'

nl = file_lines(file)

rb=dblarr(nl-1)
phib=dblarr(nl-1)
zb=dblarr(nl-1)
brb = rb
bphib=brb
bzb=brb

openr,lun,file,/get_lun
s=''
readf,lun,s
tmp=dblarr(9)
for i=0L,nl-2 do begin
   readf,lun,tmp
   rb(i) = tmp(1)
   phib(i) = tmp(2)
   zb(i) = tmp(3)
   
   brb(i) = tmp(4)*4*!dpi*1d-7
   bphib(i) =tmp(5)*4*!dpi*1d-7
   bzb(i) = tmp(6)*4*!dpi*1d-7
endfor

r = rb(uniq(rb,sort(rb)))
z = zb(uniq(zb,sort(zb)))
phi = phib(uniq(phib,sort(phib)))

nr = n_elements(r)
nt = n_elements(phi)
nz = n_elements(z)

dr = (max(r)-min(r))/(nr-1)
dt = (max(phi)-min(phi))/(nt-1)
dz = (max(z)-min(z))/(nz-1)

br=dblarr(nr,nt,nz)
bphi=br
bz=br

itot = 0L
for ir = 0,nr-1 do begin
   for it=0,nt-1 do begin
      for iz=0,nz-1 do begin
         if abs(r(ir)-rb(itot)) gt 1d-6 then print,'Uh oh: r, rb= ',r(ir),rb(itot)
         if abs(phi(it)-phib(itot)) gt 1d-6 then print,'Uh oh: phi, phib= ',phi(it),phib(itot)
         if abs(z(iz)-zb(itot)) gt 1d-6 then print,'Uh oh: z, zb= ',z(iz),zb(itot)

         br(ir,it,iz) = brb(itot)
         bphi(ir,it,iz) = bphib(itot)
         bz(ir,it,iz) = bzb(itot)
         itot += 1
      endfor
   endfor
endfor

return,{r:r,phi:phi,z:z,br:br,bphi:bphi,bz:bz,nr:nr,nt:nt,nz:nz,dr:dr,dt:dt,dz:dz,nfp:1}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function read_fieldlines_file,file

file_id = h5f_open(file)
data_id = h5d_open(file_id,'/raxis')
r = h5d_read(data_id)
data_id = h5d_open(file_id,'/zaxis')
z = h5d_read(data_id)
data_id = h5d_open(file_id,'/phiaxis')
phi = h5d_read(data_id)

data_id = h5d_open(file_id,'/B_R')
br = h5d_read(data_id)
data_id = h5d_open(file_id,'/B_Z')
bz = h5d_read(data_id)
data_id = h5d_open(file_id,'/B_PHI')
bphi = h5d_read(data_id)


dr=r(1)-r(0)
dz=z(1)-z(0)
dt=phi(1)-phi(0)

nr=n_elements(r)
nz=n_elements(z)
nt=n_elements(phi)

nfp = round(2*!dpi/max(phi))

return,{r:r,z:z,phi:phi,br:br,bz:bz,bphi:bphi,dr:dr,dz:dz,dt:dt,nr:nr,nz:nz,nt:nt,nfp:nfp}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function cleanup_fieldlines,f,ibad,jbad,kbad

nsig = 1.0
nbad = 0

fc = f

for i=1,f.nr-2 do begin
   for j=1,f.nt-2 do begin
      for k=1,f.nz-2 do begin
;for i=88,88 do begin
;   for j=30,30 do begin
;      for k=41,41 do begin
;         aver = mean([f.br(i-1,j,k),f.br(i+1,j,k),f.br(i,j-1,k),f.br(i,j+1,k),f.br(i,j,k-1),f.br(i,j,k+1)])
;         stdr = stddev([f.br(i-1,j,k),f.br(i+1,j,k),f.br(i,j-1,k),f.br(i,j+1,k),f.br(i,j,k-1),f.br(i,j,k+1)])
;         avet = mean([f.bphi(i-1,j,k),f.bphi(i+1,j,k),f.bphi(i,j-1,k),f.bphi(i,j+1,k),f.bphi(i,j,k-1),f.bphi(i,j,k+1)])
;         stdt = stddev([f.bphi(i-1,j,k),f.bphi(i+1,j,k),f.bphi(i,j-1,k),f.bphi(i,j+1,k),f.bphi(i,j,k-1),f.bphi(i,j,k+1)])
;         avez = mean([f.bz(i-1,j,k),f.bz(i+1,j,k),f.bz(i,j-1,k),f.bz(i,j+1,k),f.bz(i,j,k-1),f.bz(i,j,k+1)])
;         stdz = stddev([f.bz(i-1,j,k),f.bz(i+1,j,k),f.bz(i,j-1,k),f.bz(i,j+1,k),f.bz(i,j,k-1),f.bz(i,j,k+1)])

         aver = mean([f.br(i-1,j,k),f.br(i+1,j,k),f.br(i,j,k-1),f.br(i,j,k+1)])
         stdr = stddev([f.br(i-1,j,k),f.br(i+1,j,k),f.br(i,j,k-1),f.br(i,j,k+1)])
         avet = mean([f.bphi(i-1,j,k),f.bphi(i+1,j,k),f.bphi(i,j,k-1),f.bphi(i,j,k+1)])
         stdt = stddev([f.bphi(i-1,j,k),f.bphi(i+1,j,k),f.bphi(i,j,k-1),f.bphi(i,j,k+1)])
         avez = mean([f.bz(i-1,j,k),f.bz(i+1,j,k),f.bz(i,j,k-1),f.bz(i,j,k+1)])
         stdz = stddev([f.bz(i-1,j,k),f.bz(i+1,j,k),f.bz(i,j,k-1),f.bz(i,j,k+1)])

         isbad = (abs(f.br(i,j,k)-aver) gt nsig*stdr) or (abs(f.bphi(i,j,k)-avet) gt nsig*stdt) or (abs(f.bz(i,j,k)-avez) gt nsig*stdz)

         if isbad then begin
            if nbad lt 0.5 then begin
               ibad = i
               jbad = j
               kbad = k             
            endif else begin
               ibad = [ibad,i]
               jbad = [jbad,j]
               kbad = [kbad,k]
            endelse
            fc.br(i,j,k) = aver
            fc.bz(i,j,k) = avez
            fc.bphi(i,j,k) = avet
            nbad += 1
         endif

      endfor
   endfor
endfor

print,'Changed '+strtrim(nbad,2)+' bad point'
return,fc
end


function check_fieldlines_divB,f,ds,w=w,dsmax=dsmax

if not keyword_set(dsmax) then dsmax = 0.1

ds = dindgen(f.nr,f.nt,f.nz)

if keyword_set(w) then begin
;   p = get_vmec_surfaces(w)
;   p=replicate(p,f.nt)
   for j=0,f.nt-1 do begin
      p = get_vmec_surfaces(w,phi=f.phi(j))
      for i=0,f.nr-1 do begin
         for k=0,f.nz-1 do begin
            stmp = sqrt((f.r(i)-p.r)^2 + (f.z(k)-p.z)^2)
            ds(i,j,k) = min(stmp)
         endfor
      endfor
   endfor
endif


br=f.br
bz=f.bz
bt=f.bphi
;rbt = 0.0*bt
r3d = 0.0*f.br
for j=0,f.nt-1 do begin
   for k=0,f.nz-1 do begin
;      rbr(*,j,k) = f.r*bt(*,j,k)
      r3d(*,j,k) = f.r
   endfor
endfor

rbr = r3d*br

dr=f.dr
dz=f.dz
dt=f.dt

drbrdr = (shift(rbr,-1,0,0)-shift(rbr,1,0,0))/(2*dr)
dbzdz = (shift(bz,0,0,-1)-shift(bz,0,0,1))/(2*dz)
dbtdt = (shift(bt,0,-1,0)-shift(bt,0,1,0))/(2*dt)

divB = drbrdr/r3d + dbtdt/r3d + dbzdz

ik = where(ds gt dsmax,count)
if count gt 0 then divB(ik)=0.0d

return,divB
end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_extender_bicub,r1,phi1,z1,f,c,derivs=derivs,zero=zero

ipi = ceil(abs(phi1)/(2*!dpi/f.nfp))
phitmp = (phi1+ipi*2*!dpi/f.nfp) mod (2*!dpi/f.nfp)

dr = f.dr
dz = f.dz
dt = f.dt

ir = floor((r1-f.r(0))/dr)
iz = floor((z1-f.z(0))/dz)
it = floor((phitmp-f.phi(0))/dt)
itp1 = it+1
tt = phitmp-f.phi(it)
;if itp1 eq f.nt then itp1=0

;print,r1,z1,ir,iz

index = iz + f.nz*ir + f.nr*f.nz*it
indexp1 = iz + f.nz*ir + f.nr*f.nz*itp1

cbr = c.cbr(index,*,*)
cbz = c.cbz(index,*,*)
cbt = c.cbt(index,*,*)

cbrp1 = c.cbr(indexp1,*,*)
cbzp1 = c.cbz(indexp1,*,*)
cbtp1 = c.cbt(indexp1,*,*)


;c = g.bicub_coeffs(index,*,*)

dir = (r1 - f.r(ir))/dr
diz = (z1 - f.z(iz))/dz

;c = g.bicub_coeffs(index,*,*)

dir2 = dir^2
dir3 = dir^3
diz2 = diz^2
diz3 = diz^3

;stop

br1 = cbr(*,0,0)       + cbr(*,1,0)*dir       + cbr(*,2,0)*dir2       + cbr(*,3,0)*dir3 + $
  cbr(*,0,1)*diz   + cbr(*,1,1)*dir*diz   + cbr(*,2,1)*dir2*diz   + cbr(*,3,1)*dir3*diz + $
  cbr(*,0,2)*diz2 + cbr(*,1,2)*dir*diz2 + cbr(*,2,2)*dir2*diz2 + cbr(*,3,2)*dir3*diz2 + $
  cbr(*,0,3)*diz3 + cbr(*,1,3)*dir*diz3 + cbr(*,2,3)*dir2*diz3 + cbr(*,3,3)*dir3*diz3

bz1 = cbz(*,0,0)       + cbz(*,1,0)*dir       + cbz(*,2,0)*dir2       + cbz(*,3,0)*dir3 + $
  cbz(*,0,1)*diz   + cbz(*,1,1)*dir*diz   + cbz(*,2,1)*dir2*diz   + cbz(*,3,1)*dir3*diz + $
  cbz(*,0,2)*diz2 + cbz(*,1,2)*dir*diz2 + cbz(*,2,2)*dir2*diz2 + cbz(*,3,2)*dir3*diz2 + $
  cbz(*,0,3)*diz3 + cbz(*,1,3)*dir*diz3 + cbz(*,2,3)*dir2*diz3 + cbz(*,3,3)*dir3*diz3

bt1 = cbt(*,0,0)       + cbt(*,1,0)*dir       + cbt(*,2,0)*dir2       + cbt(*,3,0)*dir3 + $
  cbt(*,0,1)*diz   + cbt(*,1,1)*dir*diz   + cbt(*,2,1)*dir2*diz   + cbt(*,3,1)*dir3*diz + $
  cbt(*,0,2)*diz2 + cbt(*,1,2)*dir*diz2 + cbt(*,2,2)*dir2*diz2 + cbt(*,3,2)*dir3*diz2 + $
  cbt(*,0,3)*diz3 + cbt(*,1,3)*dir*diz3 + cbt(*,2,3)*dir2*diz3 + cbt(*,3,3)*dir3*diz3

br2 = cbrp1(*,0,0)       + cbrp1(*,1,0)*dir       + cbrp1(*,2,0)*dir2       + cbrp1(*,3,0)*dir3 + $
  cbrp1(*,0,1)*diz   + cbrp1(*,1,1)*dir*diz   + cbrp1(*,2,1)*dir2*diz   + cbrp1(*,3,1)*dir3*diz + $
  cbrp1(*,0,2)*diz2 + cbrp1(*,1,2)*dir*diz2 + cbrp1(*,2,2)*dir2*diz2 + cbrp1(*,3,2)*dir3*diz2 + $
  cbrp1(*,0,3)*diz3 + cbrp1(*,1,3)*dir*diz3 + cbrp1(*,2,3)*dir2*diz3 + cbrp1(*,3,3)*dir3*diz3

bz2 = cbzp1(*,0,0)       + cbzp1(*,1,0)*dir       + cbzp1(*,2,0)*dir2       + cbzp1(*,3,0)*dir3 + $
  cbzp1(*,0,1)*diz   + cbzp1(*,1,1)*dir*diz   + cbzp1(*,2,1)*dir2*diz   + cbzp1(*,3,1)*dir3*diz + $
  cbzp1(*,0,2)*diz2 + cbzp1(*,1,2)*dir*diz2 + cbzp1(*,2,2)*dir2*diz2 + cbzp1(*,3,2)*dir3*diz2 + $
  cbzp1(*,0,3)*diz3 + cbzp1(*,1,3)*dir*diz3 + cbzp1(*,2,3)*dir2*diz3 + cbzp1(*,3,3)*dir3*diz3

bt2 = cbtp1(*,0,0)       + cbtp1(*,1,0)*dir       + cbtp1(*,2,0)*dir2       + cbtp1(*,3,0)*dir3 + $
  cbtp1(*,0,1)*diz   + cbtp1(*,1,1)*dir*diz   + cbtp1(*,2,1)*dir2*diz   + cbtp1(*,3,1)*dir3*diz + $
  cbtp1(*,0,2)*diz2 + cbtp1(*,1,2)*dir*diz2 + cbtp1(*,2,2)*dir2*diz2 + cbtp1(*,3,2)*dir3*diz2 + $
  cbtp1(*,0,3)*diz3 + cbtp1(*,1,3)*dir*diz3 + cbtp1(*,2,3)*dir2*diz3 + cbtp1(*,3,3)*dir3*diz3

;dsdr1 = c.cbr(*,1,0)       + 2*c.cbr(*,2,0)*dir       + 3*c.cbr(*,3,0)*dir2 + $
;  c.cbr(*,1,1)*diz   + 2*c.cbr(*,2,1)*dir*diz   + 3*c.cbr(*,3,1)*dir2*diz + $
;  c.cbr(*,1,2)*diz2 + 2*c.cbr(*,2,2)*dir*diz2 + 3*c.cbr(*,3,2)*dir2*diz2 + $
;  c.cbr(*,1,3)*diz3 + 2*c.cbr(*,2,3)*dir*diz3 + 3*c.cbr(*,3,3)*dir2*diz3

;dsdz1 = c.cbr(*,0,1)         + c.cbr(*,1,1)*dir         + c.cbr(*,2,1)*dir2         + c.cbr(*,3,1)*dir3 + $
;  2*c.cbr(*,0,2)*diz   + 2*c.cbr(*,1,2)*dir*diz   + 2*c.cbr(*,2,2)*dir2*diz   + 2*c.cbr(*,3,2)*dir3*diz + $
;  3*c.cbr(*,0,3)*diz2 + 3*c.cbr(*,1,3)*dir*diz2 + 3*c.cbr(*,2,3)*dir2*diz2 + 3*c.cbr(*,3,3)*dir3*diz2

br = br1*(1-tt)+br2*tt
bz = bz1*(1-tt)+bz2*tt
bt = bt1*(1-tt)+bt2*tt

iroff =  where(ir lt 0 or ir gt f.nr-2,count)
if count gt 0 then begin
    indroff = iz(iroff) + f.nr*ir(iroff)
    br(indroff) = 0.0d
    bz(indroff) = 0.0d         ;
    bt(indroff) = 1d10
endif

izoff =  where(iz le 1 or iz ge f.nz-2,count)
if count gt 0 then begin
    indzoff = iz(izoff) + f.nr*ir(izoff)
    br(indzoff) = 0.0d
    bz(indzoff) = 0.0d         ;
    bt(indzoff) = 1d10
endif

;; if keyword_set(derivs) then begin
;;     d2sdr21 =  2*c(*,2,0)       + 6*c(*,3,0)*dir + $
;;       2*c(*,2,1)*diz   + 6*c(*,3,1)*dir*diz + $
;;       2*c(*,2,2)*diz2 + 6*c(*,3,2)*dir*diz2 + $
;;       2*c(*,2,3)*diz3 + 6*c(*,3,3)*dir*diz3

;;     d2sdz21 =  2*c(*,0,2)   + 2*c(*,1,2)*dir   + 2*c(*,2,2)*dir2   + 2*c(*,3,2)*dir3 + $
;;       6*c(*,0,3)*diz + 6*c(*,1,3)*dir*diz + 6*c(*,2,3)*dir2*diz + 6*c(*,3,3)*dir3*diz

;;     d2sdrdz1 =  c(*,1,1)       + 2*c(*,2,1)*dir       + 3*c(*,3,1)*dir2 + $
;;       2*c(*,1,2)*diz   + 4*c(*,2,2)*dir*diz   + 6*c(*,3,2)*dir2*diz + $
;;       3*c(*,1,3)*diz2 + 6*c(*,2,3)*dir*diz2 + 9*c(*,3,3)*dir2*diz^2
    
;;     d2sdr21 = d2sdr21/(dr*dr)
;;     d2sdz21 = d2sdz21/(dz*dz)
;;     d2sdrdz1 = d2sdrdz1/(dr*dz)

;;     dbrdr1 = -br1/r1 - d2sdrdz1/r1
;;     dbrdz1 = -d2sdz21/r1
;;     dbzdr1 = -bz1/r1 + d2sdr21/r1
;;     dbzdz1 = d2sdrdz1/r1

;;     dfpol = polyder(psiN,g.fpol_coeffs)
;;     dfpdr = dfpol*(-dsdr1*g.cpasma/abs(g.cpasma))/(g.ssibry-g.ssimag)
;;     dfpdz = dfpol*(-dsdz1*g.cpasma/abs(g.cpasma))/(g.ssibry-g.ssimag)
;;     dbtdr1 = -bt1/r1
;;     dbtdz1 = 0.0*dbtdr1
;;     if ic(0) ne -1 then dbtdr1(ic) = dbtdr1(ic) + dfpdr(ic)/r1(ic)
;;     if ic(0) ne -1 then dbtdz1(ic) = dfpdz(ic)/r1(ic)

;;     return,{br:br1,bz:bz1,bphi:bt1,dbrdr:dbrdr1,dbrdz:dbrdz1,dbzdr:dbzdr1,dbzdz:dbzdz1, $
;;             dbphidr:dbtdr1,dbphidz:dbtdz1}
;; endif else return,{br:br1,bz:bz1,bphi:bt1}
;; endelse

return,{br:br,bz:bz,bphi:bt}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_extender,r,phi,z,a,biext=biext

if keyword_set(biext) then begin
   return,bfield_extender_bicub(r,phi,z,a,biext)
endif else begin

ipi = ceil(abs(phi)/(2*!dpi/a.nfp))

phitmp = (phi+ipi*2*!dpi/a.nfp) mod (2*!dpi/a.nfp)

ir = floor((r-a.r(0))/a.dr)
it = floor((phitmp-a.phi(0))/a.dt)
iz = floor((z-a.z(0))/a.dz)

offgrid = where(ir gt a.nr-2 or it gt a.nt-2 or iz gt a.nz-2 or ir lt 0 or it lt 0 or iz lt 0,count)
if count gt 0 then begin
;   print,[phi,phitmp,ipi]
;   print,'Offgrid, r,phi,z= ',[r,phi,z]
;   print,'Offgrid, ir, it, iz= ',[ir,it,iz]
   ir(offgrid)=0
   iz(offgrid)=0
   it(offgrid)=0
endif
;if ir gt a.nr-2 or it gt a.nt-2 or iz gt a.nz-2 then begin
;   print,'Point off of extender grid, dummy';
;stop
;   return,{br:0.0+0*r,bphi:1d10+0*r,bz:0.0+0*r}
;endif else begin
tr = r-a.r(ir)
tt = phitmp-a.phi(it)
tz = z-a.z(iz)

;   print,a.r(ir),a.phi(it),a.z(iz),tr,tt,tz

   br = a.br(ir,it,iz)*(1-tr)*(1-tt)*(1-tz) + a.br(ir+1,it,iz)*(tr)*(1-tt)*(1-tz) + a.br(ir,it+1,iz)*(1-tr)*(tt)*(1-tz) + a.br(ir,it,iz+1)*(1-tr)*(1-tt)*(tz)  $
   + a.br(ir+1,it+1,iz)*(tr)*(tt)*(1-tz) + a.br(ir+1,it,iz+1)*(tr)*(1-tt)*(tz) + a.br(ir,it+1,iz+1)*(1-tr)*(tt)*(tz) + a.br(ir+1,it+1,iz+1)*(tr)*(tt)*(tz)
   bphi = a.bphi(ir,it,iz)*(1-tr)*(1-tt)*(1-tz) + a.bphi(ir+1,it,iz)*(tr)*(1-tt)*(1-tz) + a.bphi(ir,it+1,iz)*(1-tr)*(tt)*(1-tz) + a.bphi(ir,it,iz+1)*(1-tr)*(1-tt)*(tz)  $
   + a.bphi(ir+1,it+1,iz)*(tr)*(tt)*(1-tz) + a.bphi(ir+1,it,iz+1)*(tr)*(1-tt)*(tz) + a.bphi(ir,it+1,iz+1)*(1-tr)*(tt)*(tz) + a.bphi(ir+1,it+1,iz+1)*(tr)*(tt)*(tz) 
   bz = a.bz(ir,it,iz)*(1-tr)*(1-tt)*(1-tz) + a.bz(ir+1,it,iz)*(tr)*(1-tt)*(1-tz) + a.bz(ir,it+1,iz)*(1-tr)*(tt)*(1-tz) + a.bz(ir,it,iz+1)*(1-tr)*(1-tt)*(tz)  $
   + a.bz(ir+1,it+1,iz)*(tr)*(tt)*(1-tz) + a.bz(ir+1,it,iz+1)*(tr)*(1-tt)*(tz) + a.bz(ir,it+1,iz+1)*(1-tr)*(tt)*(tz) + a.bz(ir+1,it+1,iz+1)*(tr)*(tt)*(tz) 

;print,'1: ',[r,phi,z,br,bphi,bz]
;print,'2: ',[ir,it,iz,tr,tt,tz]

if count gt 0 then begin
   br(offgrid) = 0.0
   bphi(offgrid) = 1d10
   bz(offgrid) = 0.0
endif

return,{br:br,bphi:bphi,bz:bz,ir:ir,it:it,iz:iz,tr:tr,tt:tt,tz:tz}
endelse
end



function bfield_extender_ave,r,phi,z,a

x=r*cos(phi)
y=r*sin(phi)
dx = 0.2*a.dr
x1=x+dx
x2=x-dx
y3=y+dx
y4=y-dx
z5=z+dx
z6=z-dx

r1 = sqrt(x1^2+y^2)
r2 = sqrt(x2^2+y^2)
r3 = sqrt(x^2+y3^2)
r4 = sqrt(x^2+y4^2)
phi1 = atan(y,x1)
phi2 = atan(y,x2)
phi3 = atan(y3,x)
phi4 = atan(y4,x)

b1=bfield_extender(r1,phi1,z,a)
b2=bfield_extender(r2,phi2,z,a)
b3=bfield_extender(r3,phi3,z,a)
b4=bfield_extender(r4,phi4,z,a)
b5=bfield_extender(r,phi,z5,a)
b6=bfield_extender(r,phi,z6,a)

br = (b1.br + b2.br + b3.br + b4.br + b5.br + b6.br)/6.0d
bphi = (b1.bphi + b2.bphi + b3.bphi + b4.bphi + b5.bphi + b6.bphi)/6.0d
bz = (b1.bz + b2.bz + b3.bz + b4.bz + b5.bz + b6.bz)/6.0d

return,{br:br,bphi:bphi,bz:bz}
end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function find_xpt,g,second=second,refine=refine,tol=tol

if not keyword_set(tol) then tol=1e-6

ik = where (g.bdry(0,*) gt 1e-4)
rb = g.bdry(0,ik)
zb = g.bdry(1,ik)

b=bfield_geq_bicub(g,rb,zb)
bpol = sqrt(b.br^2+b.bz^2)

junk = min(bpol,ix)
rx = rb(ix)
zx = zb(ix)
print,bpol(ix),rx,zx


if keyword_set(refine) then begin
err = 1e-1
n1 = 100
bp = dblarr(100,100)
rg = bp
zg = bp
br=bp
bz=bp

while err gt tol do begin
    de = err
     rt = rx-0.5*de+de*dindgen(100)/99
     zt = zx-0.5*de+de*dindgen(100)/99
     for i = 0,n1-1 do begin
         b=bfield_geq_bicub(g,rt,replicate(zt(i),n1))
         rg(i,*) = rt
         zg(i,*) = replicate(zt(i),n1)
         bp(i,*) = sqrt(b.br^2+b.bz^2)
     endfor
     err = sqrt((rg(ix)-rx)^2 + (zg(ix)-zx)^2)
     junk=min(bp,ix)
     rx = rg(ix)
     zx = zg(ix)
     bpx = bp(ix)
;     print,err
 endwhile
;junk=min(bp,ix)
print,bpx,rx,zx
endif
rx2=rx
zx2=zx

if keyword_set(second) then begin
rx2=rx
zx2=-zx
err = 1e-1
n1 = 100
bp = dblarr(n1,n1)
rg = bp
zg = bp
br=bp
bz=bp

while err gt tol do begin
    de = err
     rt = rx2-0.5*de+de*dindgen(n1)/(n1-1)
     zt = zx2-0.5*de+de*dindgen(n1)/(n1-1)
     for i = 0,n1-1 do begin
         b=bfield_geq_bicub(g,rt,replicate(zt(i),n1))
         rg(i,*) = rt
         zg(i,*) = replicate(zt(i),n1)
         bp(i,*) = sqrt(b.br^2+b.bz^2)
     endfor
     err = sqrt((rg(ix)-rx2)^2 + (zg(ix)-zx2)^2)
     junk=min(bp,ix)
     rx2 = rg(ix)
     zx2 = zg(ix)
     bpx2 = bp(ix)
;     print,err,rx2,zx2,bpx2
;err=1e-8
 endwhile
;junk=min(bp,ix)
print,bpx2,rx2,zx2
endif

return,{rx:rx,zx:zx,rx2:rx2,zx2:zx2}
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_d3d_icoils,taper,ntorpts=ntorpts,cenzero=cenzero
if n_elements(taper) eq 1 then taper = taper*[1,-1,1,-1,1,-1,1,-1,1,-1,1,-1]
if not keyword_set(ntorpts) then ntorpts = 5
npts = 2*ntorpts+1
phicens = double([-032.7,-087.3,-152.7,-207.3,-272.7,-327.3])*!dpi/180.0d
phiext = double(51.72)*!dpi/180.0d
R = double([2.184,2.394])
Z = double([1.012,0.504])

coil = dblarr(12*npts,3)
current = dblarr(12*npts)

for i = 0,5 do begin
    phi = phiext*dindgen(ntorpts)/(ntorpts-1)+phicens(i)
    if keyword_set(cenzero) then phi += -0.5d*phiext

    coil(i*npts:i*npts+ntorpts-1,0) = R(0)*cos(phi)
    coil(i*npts:i*npts+ntorpts-1,1) = R(0)*sin(phi)
    coil(i*npts:i*npts+ntorpts-1,2) = Z(0)
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,0) = R(1)*cos(reverse(phi))
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,1) = R(1)*sin(reverse(phi))
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,2) = Z(1)
    coil((i+1)*npts-1,*) = coil(i*npts,*)
    current(i*npts:(i+1)*npts-1) = taper(i)
    current((i+1)*npts-1) = 0.0d
    
    coil((i+6)*npts:(i+7)*npts-1,0) = reverse(coil(i*npts:(i+1)*npts-1,0))
    coil((i+6)*npts:(i+7)*npts-1,1) = reverse(coil(i*npts:(i+1)*npts-1,1))
    coil((i+6)*npts:(i+7)*npts-1,2) = -reverse(coil(i*npts:(i+1)*npts-1,2))
    current((i+6)*npts:(i+7)*npts-1) = taper(i+6)
    current((i+6)*npts-1) = 0.0d


endfor

return,{coil:coil,current:current}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_d3d_ccoils,taper,ntorpts=ntorpts
if n_elements(taper) eq 1 then taper = taper*[1,-1,1,-1,1,-1]
if not keyword_set(ntorpts) then ntorpts = 5
npts = 2*ntorpts+1
phicens = double([0.0,60.,120.,180.,240.,300.])*!dpi/180.0d
phiext = double(59.0)*!dpi/180.0d
R = double([3.2,3.2])
Z = double([0.8,-0.8])

coil = dblarr(6*npts,3)
current = dblarr(6*npts)

for i = 0,5 do begin
    phi = phiext*dindgen(ntorpts)/(ntorpts-1)+phicens(i)

    coil(i*npts:i*npts+ntorpts-1,0) = R(0)*cos(phi)
    coil(i*npts:i*npts+ntorpts-1,1) = R(0)*sin(phi)
    coil(i*npts:i*npts+ntorpts-1,2) = Z(0)
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,0) = R(1)*cos(reverse(phi))
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,1) = R(1)*sin(reverse(phi))
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,2) = Z(1)
    coil((i+1)*npts-1,*) = coil(i*npts,*)
    current(i*npts:(i+1)*npts-1) = taper(i)
    current((i+1)*npts-1) = 0.0d
endfor

return,{coil:coil,current:current}
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_nstx_rwmcoils,taper,ntorpts=ntorpts
if n_elements(taper) eq 1 then taper = taper*[1,-1,1,-1,1,-1]
if not keyword_set(ntorpts) then ntorpts = 5
npts = 2*ntorpts+1
phicens = double([0.0,60.,120.,180.,240.,300.])*!dpi/180.0d
phiext = double(56.0)*!dpi/180.0d
R = double([1.76,1.76])
Z = double([0.4826,-0.4826])

coil = dblarr(6*npts,3)
current = dblarr(6*npts)

for i = 0,5 do begin
    phi = phiext*dindgen(ntorpts)/(ntorpts-1)+phicens(i)

    coil(i*npts:i*npts+ntorpts-1,0) = R(0)*cos(phi)
    coil(i*npts:i*npts+ntorpts-1,1) = R(0)*sin(phi)
    coil(i*npts:i*npts+ntorpts-1,2) = Z(0)
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,0) = R(1)*cos(reverse(phi))
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,1) = R(1)*sin(reverse(phi))
    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,2) = Z(1)
    coil((i+1)*npts-1,*) = coil(i*npts,*)
    current(i*npts:(i+1)*npts-1) = taper(i)
    current((i+1)*npts-1) = 0.0d
endfor

return,{coil:coil,current:current}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_nstx_ncc,taper,ntorpts=ntorpts,even=even,ltaper=ltaper
if n_elements(taper) eq 1 then taper = taper*[1,1,-1,-1,1,1,-1,-1,1,1,-1,-1]
if keyword_set(even) then lfac = 1.0d else lfac=-1.0d
if not keyword_set(ltaper) then ltaper=taper

xpu=[1.44641,1.47433,1.48366,1.47433,1.44641,1.40607,1.36573,1.32538,1.35096,1.35951,1.35096,1.32538,1.36572,1.40607,1.44641]
ypu=[-0.330343,-0.166154,0,0.166154,0.330343,0.321129,0.311915,0.302701,0.15225,0,-0.15225,-0.302701,-0.311914,-0.32113,-0.330343]
zpu=[0.586503,0.586503,0.586503,0.586503,0.586503,0.712036,0.83757,0.963103,0.963103,0.963103,0.963103,0.963103,0.83757,0.712036,0.586503]

xpl=[1.32538,1.35096,1.35951,1.35096,1.32538,1.36573,1.40607,1.44641,1.47433,1.48366,1.47433,1.44641,1.40607,1.36572,1.32538]
ypl=[-0.302701,-0.15225,0,0.15225,0.302701,0.311915,0.321129,0.330343,0.166154,0,-0.166154,-0.330343,-0.32113,-0.311914,-0.302701]
zpl=[-0.963103,-0.963103,-0.963103,-0.963103,-0.963103,-0.83757,-0.712036,-0.586503,-0.586503,-0.586503,-0.586503,-0.586503,-0.712036,-0.83757,-0.963103]

rpu = sqrt(xpu^2+ypu^2)
phipu = atan(ypu,xpu)

rpl = sqrt(xpl^2+ypl^2)
phipl = atan(ypl,xpl)

npts = n_elements(rpu)

phicens = (30*dindgen(12)+15.0d)*!dpi/180.0d
coil = dblarr(24*npts,3)
current = dblarr(24*npts)

for i=0,11 do begin
   coil(i*npts:(i+1)*npts-1,0) = rpu*cos(phipu+phicens(i))
   coil(i*npts:(i+1)*npts-1,1) = rpu*sin(phipu+phicens(i))
   coil(i*npts:(i+1)*npts-1,2) = zpu
   current(i*npts:(i+1)*npts-2) = taper(i)

   coil((i+12)*npts:(i+12+1)*npts-1,0) = rpl*cos(phipl+phicens(i))
   coil((i+12)*npts:(i+12+1)*npts-1,1) = rpl*sin(phipl+phicens(i))
   coil((i+12)*npts:(i+12+1)*npts-1,2) = zpl
   current((i+12)*npts:(i+12+1)*npts-2) = lfac*ltaper(i)
endfor

return,{coil:coil,current:current}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_nstx_tfcoils,curr=curr
if not keyword_set(curr) then curr = 55e3
restore,'/home/6zc/VMEC_RUNS/NSTX/MGRID_COILS/NSTX_3d_coil_sensor_system.dat'
c = coils.wd0
s=size(c.tf.x,/dimension)
ntfcoils = s(0)
print,ntfcoils
ntfpts = s(1)
Rtf = sqrt(c.tf.x(0,*)^2+c.tf.y(0,*)^2)
Ztf = coils.wd0.tf.z(0,*)
phitf = 2*!dpi*findgen(ntfcoils)/ntfcoils

Rtf(ntfpts-1) = Rtf(0)
Ztf(ntfpts-1) = Ztf(0)

coil = dblarr(ntfcoils*ntfpts,3)
current = dblarr(ntfcoils*ntfpts)

for i = 0,ntfcoils-1 do begin
    coil(i*ntfpts:(i+1)*ntfpts-1,0) = Rtf*cos(phitf(i))
    coil(i*ntfpts:(i+1)*ntfpts-1,1) = Rtf*sin(phitf(i))
    coil(i*ntfpts:(i+1)*ntfpts-1,2) = Ztf
    current(i*ntfpts:(i+1)*ntfpts-1) = curr
    current((i+1)*ntfpts-1) = 0.0d
    
;    coil(i*npts:i*npts+ntorpts-1,0) = R(0)*cos(phi)
;    coil(i*npts:i*npts+ntorpts-1,1) = R(0)*sin(phi)
;    coil(i*npts:i*npts+ntorpts-1,2) = Z(0)
;    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,0) = R(1)*cos(reverse(phi))
;    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,1) = R(1)*sin(reverse(phi))
;    coil(i*npts+ntorpts:i*npts+2*ntorpts-1,2) = Z(1)
;    coil((i+1)*npts-1,*) = coil(i*npts,*)
;    current(i*npts:(i+1)*npts-1) = taper(i)
;    current((i+1)*npts-1) = 0.0d
endfor

return,{coil:coil,current:current}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_nstx_pf5_error,current
; Current is PF5 amps per turn
nturns=24
restore,'~/idl/bfield/NSTX_FILES/NSTX_3d_coil_sensor_system.dat'
;help,coils.wd0.pf5u,/str

;z=mycolors()
;plot,coils.wd0.pf5l.x,coils.wd0.pf5l.y,/iso
;oplot,coils.wd1.pf5l.x,coils.wd1.pf5l.y,color=z.brick
nphi = n_elements(coils.wd0.pf5u.x)
newcoil = dblarr(4*nphi,3)
phinds = dindgen(nphi)
newcoil(phinds,0) = coils.wd0.pf5u.x
newcoil(phinds,1) = coils.wd0.pf5u.y
newcoil(phinds,2) = coils.wd0.pf5u.z
newcoil(nphi+phinds,0) = coils.wd0.pf5l.x
newcoil(nphi+phinds,1) = coils.wd0.pf5l.y
newcoil(nphi+phinds,2) = coils.wd0.pf5l.z
newcoil(2*nphi+phinds,0) = coils.wd1.pf5u.x
newcoil(2*nphi+phinds,1) = coils.wd1.pf5u.y
newcoil(2*nphi+phinds,2) = coils.wd1.pf5u.z
newcoil(3*nphi+phinds,0) = coils.wd1.pf5l.x
newcoil(3*nphi+phinds,1) = coils.wd1.pf5l.y
newcoil(3*nphi+phinds,2) = coils.wd1.pf5l.z
phinds = dindgen(nphi-1)
newcurrent=dblarr(4*nphi)
newcurrent(phinds) = -nturns*current
newcurrent(nphi+phinds) = -nturns*current
newcurrent(2*nphi+phinds) = nturns*current
newcurrent(3*nphi+phinds) = nturns*current

return,{coil:newcoil,current:newcurrent}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function read_mgrid_file,filename

openr,1,filename
s = ''
readf,1,s
;print,s
readf,1,s
;print,s
readf,1,s
;print,s

count = 0L
while s ne 'end' do begin
readf,1,s
count = count+1
endwhile
count = count-1
;print,count

close,1
openr,1,filename
s = ''
readf,1,s
;print,s
readf,1,s
;print,s
readf,1,s
;print,s
tmp = dblarr(4)
for i=0L,count-1 do begin
    readf,1,tmp
    if n_elements(current) eq 0 then begin
        coil = tmp(0:2)
        current = tmp(3)
    endif else begin
        coil = [[coil],[tmp(0:2)]]
        current = [current,tmp(3)]
    endelse
endfor

close,1
return,{coil:transpose(coil),current:current}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_b_elements,p,derivs=derivs,degree=degree

if not keyword_set(degree) then degree = 4

x=double(p(0))
y=double(p(1))
z=double(p(2))

nelem = [3,8,15,24,35]
nelem = nelem(floor(degree))

bx=dblarr(nelem)
by=bx & bz=bx

if keyword_set(derivs) then dbxdx=bx & dbxdy=bx & dbxdz=bx & dbydx=bx & dbydy=bx & dbydz=bx & dbzdx=bx & dbzdy=bx & dbzdz=bx 

; Constant elements:
; 0. a000=1 : bx=1, db*d* = 0
; 1. b000=1 : by=1; db*d* = 0
; 2. c000=1 : bz=1; db*d* = 0
bx(0) = 1
by(1) = 1
bz(2) = 1

if degree ge 1 then begin
; Linear elements
; 3. b001=1, c010=1: by=z, bz=y, dbydz=1, dbzdy=1
; 4. a001=1, c100=1: bx=z, bz=x, dbxdz=1, dbzdx=1
; 5. a010=1, b100=1, bx=y, by=x, dbxdy=1, dbydx=1
; 6. a100=1, c001=-1, bx=x, bz=-z, dbxdx=1, dbzdz=-1
; 7. a100=1, b010=-1, bx=x, by=-y, dbxdx=1, dbydy=-1
    by(3)=z & bz(3)=y
    bx(4)=z & bz(4)=x
    bx(5)=y & by(5)=x
    bx(6)=x & bz(6)=-z
    bx(7)=x & by(7)=-y
    if keyword_set(derivs) then begin
        dbydz(3) = 1 & dbzdy(3)=1
        dbxdz(4) = 1 & dbzdx(4)=1
        dbxdy(5) = 1 & dbydx(5)=1
        dbxdx(6) = 1 & dbzdz(6)=-1
        dbxdx(7) = 1 & dbydy(7)=-1
    endif
endif

if degree ge 2 then begin
; Quadratic elements
; 8. a011=1, b101=1,c110=1: bx=yz,by=xz,bz=xy, 
;                           dbxdy=z,dbxdz=y,dbydx=z,dbydz=x,dbzdx=y,dbzdy=x
; 9. a020=1,a200=-1,b110=2: bx=y^2-x^2, by=2xy
;                           dbxdx=-2x, dbxdy=2y, dbydx=2y,dbydy=2x
;10. a002=1,a200=-1,c101=2: bx=z^2-x^2,bz=2xz
;                           dbxdx=-2x, dbxdz=2z, dbzdx=2z,dbzdz=2x
;11. b002=1,b020=-1,c011=2: by=z^2-y^2,bz=2yz
;                           dbydy=-2y, dbydz=2z, dbzdy=2z, dbzdz=2y
;12. a110=2,b200=1,b020=-1: bx=2xy, by=x^2-y^2
;                           dbxdx=2y, dbxdy=2x, dbydx=2x, dbydy=-2y
;13. b011=2,c020=1,c002=-1: by=2yz, bz=y^2-z^2
;                           dbydy=2z, dbydz=2y, dbzdy=2y, dbzdz=-2z
;14. a101=2,c200=1,c002=-1: bx=2xz, bz=x^2-z^2
;                           dbxdx=2z, dbxdz=2x, dbzdx=2x, dbzdz=-2z
;    bx()= & by()= & bz()=
    bx(8)=y*z & by(8)=x*z & bz(8)=x*y
    bx(9)= y^2-x^2 & by(9)= 2*x*y
    bx(10)= z^2-x^2 & bz(10)=2*x*z
    by(11)= z^2-y^2 & bz(11)= 2*y*z
    bx(12)= 2*x*y & by(12)= x^2-y^2
    by(13)= 2*y*z & bz(13)= y^2-z^2
    bx(14)= 2*x*z & bz(14)= x^2-z^2
    if keyword_set(derivs) then begin
        dbxdy(8)=z & dbxdz(8)=y & dbydx(8)=z & dbydz(8)=x & dbzdx(8)=y & dbzdy(8)=x
        dbxdx(9)=-2*x & dbxdy(9)=2*y & dbydx(9)=2*y & dbydy(9)=2*x
        dbxdx(10)=-2*x & dbxdz(10)=2*z & dbzdx(10)=2*z & dbzdz(10)=2*x
        dbydy(11)=-2*y & dbydz(11)=2*z & dbzdy(11)=2*z & dbzdz(11)=2*y
        dbxdx(12)=2*y & dbxdy(12)=2*x & dbydx(12)=2*x & dbydy(12)=-2*y
        dbydy(13)=2*z & dbydz(13)=2*y & dbzdy(13)=2*y & dbzdz(13)=-2*z
        dbxdx(14)=2*z & dbxdz(14)=2*x & dbzdx(14)=2*x & dbzdz(14)=-2*z
    endif
endif

if degree ge 3 then begin
; Cubic elements
;15. a111=6,b021=-3,b201=3,c030=-1,c210=3:
;      bx=6xyz, by=-3y^2z+3x^2z, bz=-y^3+3x^2y
;      dbxdx=6yz, dbxdy=6xz, dbxdz=6xy, dbydx=6xz, dbydy=-6xz,
;      dbydz=-3y^2+3x^2, dbzdx=6xy, dbzdy=-3y2+3x^2
;16. a111=6,b003=-1,b201=3,c210=3,c012=-3
;      bx=6xyz, by=-z^3+3x^2z, bz=3x^2y-3yz^2
;      dbxdx=6yz, dbxdy=6xz, dbxdz=6xy, dbybx=6xz, dbydz=-3z^2+3x^2
;      dbzdx=6xy dbzdy=3x^2-3z^2, dbzdz=-6yz
;17. a210=3,a012=-3,b300=1,b102=-3,c111=-6
;      bx=3x^2*y-3yz^2,by=x^3-3xz^2,bz=-6xyz
;      dbxdx=6xy, dbxdy=3x^2-3z^2, dbxdz=-6yz, dbydx=3x^2-3z^2
;      dbydz=-6xz, dbzdx=-6yz, dbzdy=-6xz, dbzdz=-6xy
;18. a030=1,a012=-3,b120=3,b102=-3,c111=-6
;      bx=y^3-3yz^2, by=3xy^2-3xz^2,bz=-6xyz
;      dbxdy=3y^2-3z^2,dbxdz=-6yz,dbydx=3y^2-3z^2,
;      dbydy=6xy,dbydz=-6xz, dbzdx=-6yz, dbzdy=-6xz, dbzdz=-6xy
;19. a300=2,a120=-3,a102=-3,b210=-3,b012=3,c201=-3,c021=3
;      bx=2x^3-3xy^2-3xz^2, by=-3x^2y+3yz^2, bz=-3x^2z+3y^2z
;      dbxdx=6x^2-3y^2-3z^2, dbxdy=-6xy, dbxdz=-6xz, dbydx=-6xy
;      dbydy=-3x^2+3z^2, dbydz=6yz, dbzdx=-6xz, dbzdy=6yz,
;      dbzdz=-3x^2+3y^2
;20. a102=3,a120=-3,b030=2,b012=-3,b210=-3,c201=3,c021=-3
;      bx=3xz^2-3xy^2, by=2y^3-3yz^2-3x^2y, bz=3x^2z-3y^2z
;      dbxdx=3z^2-3y^2, dbxdy=-6xy, dbxdz=6xz, dbydx=-6xy,
;      dbydy=6y^2-3z^2-3x^2, dbydz=-6yz, dbzdx=6xz, dbzdy=-6yz,
;      dbzdz=3x^2-3y^2
;21. a120=3,a102=-3,b210=3,b012=-3,c003=2,c201=-3,c021=-3
;      bx=3xy^2-3xz^2, by=3x^2y-3yz^2, bz=2z^3-3x^2z-3y^2z
;      dbxdx=3y^2-3z^2, dbxdy=6xy, dbxdz=-6xz, dbydx=6xy, dbydy=3x^2-3z^2
;      dbydz=-6yz, dbzdx=-6xz, dbzdy=-6yz, dbzdz=6z^2-3x^2-3y^2
;22. a003=1,a021=-3,b111=-6,c102=3,c120=-3
;      bx=z^3-3y^2z, by=-6xyz, bz=3xz^2-3xy^2
;      dbxdy=-6yz, dbxdz=3z^2-3y^2, dbydx=-6yz, dbydy=-6xz,
;      dbydz=-6xy, dbzdx=3z^2-3y^2, dbzdy=-6xy, dbzdz=6xz
;23. a201=-3,a021=3,b111=6,c300=-1,c120=3
;      bx=-3x^2z+3y^2z, by=6xyz, bz=-x^3+3xy^2
;      dbxdx=-6xz, dbxdy=6yz, dbxdz=-3x^2+3y^2, dbydx=6yz
;      dbydy=6xz, dbydz=6xy, dbzdx=-3x^2+3y^2, dbzdy=6xy
    bx(15)=6*x*y*z & by(15)=-3*y^2*z+3*x^2*z & bz(15)=-y^3+3*x^2*y
    bx(16)=6*x*y*z & by(16)=-z^3+3*x^2*z & bz(16)=3*x^2*y-3*y*z^2
    bx(17)=3*x^2*y-3*y*z^2 & by(17)=x^3-3*x*z^2 & bz(17)=-6*x*y*z
    bx(18)=y^3-3*y*z^2 & by(18)=3*x*y^2-3*x*z^2 & bz(18)=-6*x*y*z
    bx(19)=2*x^3-3*x*y^2-3*x*z^2 & by(19)=-3*x^2*y+3*y*z^2 & bz(19)=-3*x^2*z+3*y^2*z
    bx(20)=3*x*z^2-3*x*y^2 & by(20)=2*y^3-3*y*z^2-3*x^2*y & bz(20)=3*x^2*z-3*y^2*z
    bx(21)=3*x*y^2-3*x*z^2 & by(21)=3*x^2*y-3*y*z^2 & bz(21)=2*z^3-3*x^2*z-3*y^2*z
    bx(22)=z^3-3*y^2*z & by(22)=-6*x*y*z & bz(22)=3*x*z^2-3*x*y^2
    bx(23)=-3*x^2*z+3*y^2*z & by(23)=6*x*y*z & bz(23)=-x^3+3*x*y^2
    if keyword_set(derivs) then begin
        dbxdx(15)=6*y*z & dbxdy(15)=6*x*z & dbxdz(15)=6*x*y & dbydx(15)=6*x*z & dbydy(15)=-6*y*z
        dbydz(15)=-3*y^2+3*x^2 & dbzdx(15)=6*x*y & dbzdy(15)=-3*y^2+3*x^2
        dbxdx(16)=6*y*z & dbxdy(16)=6*x*z & dbxdz(16)=6*x*y & dbydx(16)=6*x*z & dbydz(16)=-3*z^2+3*x^2
        dbzdx(16)=6*x*y & dbzdy(16)=3*x^2-3*z^2 & dbzdz(16)=-6*y*z
        dbxdx(17)=6*x*y & dbxdy(17)=3*x^2-3*z^2 & dbxdz(17)=-6*y*z & dbydx(17)=3*x^2-3*z^2
        dbydz(17)=-6*x*z & dbzdx(17)=-6*y*z & dbzdy(17)=-6*x*z & dbzdz(17)=-6*x*y
        dbxdy(18)=3*y^2-3*z^2 & dbxdz(18)=-6*y*z & dbydx(18)=3*y^2-3*z^2 
        dbydy(18)=6*x*y &dbydz(18)=-6*x*z & dbzdx(18)=-6*y*z & dbzdy(18)=-6*x*z & dbzdz(18)=-6*x*y
        dbxdx(19)=6*x^2-3*y^2-3*z^2 & dbxdy(19)=-6*x*y & dbxdz(19)=-6*x*z & dbydx(19)=-6*x*y
        dbydy(19)=-3*x^2+3*z^2 & dbydz(19)=6*y*z & dbzdx(19)=-6*x*z & dbzdy(19)=6*y*z & dbzdz(19)=-3*x^2+3*y^2
        dbxdx(20)=3*z^2-3*y^2 & dbxdy(20)=-6*x*y & dbxdz(20)=6*x*z & dbydx(20)=-6*x*y
        dbydy(20)=6*y^2-3*z^2-3*x^2 & dbydz(20)=-6*y*z & dbzdx(20)=6*x*z & dbzdy(20)=-6*y*z & dbzdz(20)=3*x^2-3*y^2
        dbxdx(21)=3*y^2-3*z^2 & dbxdy(21)=6*x*y & dbxdz(21)=-6*x*z & dbydx(21)=6*x*y & dbydy(21)=3*x^2-3*z^2
        dbydz(21)=-6*y*z & dbzdx(21)=-6*x*z & dbzdy(21)=-6*y*z & dbzdz(21)=6*z^2-3*x^2-3*y^2
        dbxdy(22)=-6*y*z & dbxdz(22)=3*z^2-3*y^2 & dbydx(22)=-6*y*z & dbydy(22)=-6*x*z 
        dbydz(22)=-6*x*y & dbzdx(22)=3*z^2-3*y^2 & dbzdy(22)=-6*x*y & dbzdz(22)=6*x*z
        dbxdx(23)=-6*x*z & dbxdy(23)=6*y*z & dbxdz(23)=-3*x^2+3*y^2 & dbydx(23)=6*y*z
        dbydy(23)=6*x*z & dbydz(23)=6*x*y & dbzdx(23)=-3*x^2+3*y^2 & dbzdy(23)=6*x*y
    endif 
endif

if degree ge 4 then begin
; Quartic elements
;24. a400=1,a202=-3,a220=-3,a022=3,b310=-2,b112=6,c301=-2,c121=6
;     bx=x^4-3x^2z^2-3x^2y^2+3y^2z^2, by=-2x^3y+6xyz^2, bz=-2x^3z+6xy^2z
;     dbxdx=4x^3-6xy^2-6xz^2, dbxdy=-6x^2y+6yz^2, dbxdz=-6x^2z+6y^2z
;     dbydx=-6x^2y+6yz^2, dbydy=-2x^3+6xz^2, dbydz=12xyz
;     dbzdx=-6x^2z+6y^2z, dbzdy=12xyz, dbzdz=-2x^2+6xy^2
;25. a040=1,a202=3,a220=-3,a022=-3,b310=-2,b130=4,b112=-6,c301=2,c121=-6
;     bx=y^4+3x^2z^2-3x^2y^2-3y^2z^2, by=-2x^3y+4xy^3-6xyz^2, bz=2x^3z-6xy^2z
;     dbxdx=6xz^2-6xy^2, dbxdy=4y^3-6x^2y-6yz^2, dbxdz=6x^2z-6y^2z
;     dbydx=-6x^2y+4y^3-6yz^2, dbydy=-2x^3+12xy^2-6xz^2, dbydz=-12xyz
;     dbzdx=6x^2z-6y^2z, dbzdy=-12xyz, dbzdz=2x^3-6xy^2
;26. a004=1,a202=-3,a220=3,a022=-3,b310=2,b112=-6,c301=-2,c121=-6,c103=4
;     bx=z^4-3x^2z^2+3x^2y^2-3y^2z^2, by=2x^3y-6xyz^2, bz=-2x^3z-6xy^2z+4xz^3
;     dbxdx=-6xz^2+6xy^2, dbxdy=6x^2y-6yz^2, dbxdz=4z^3-6x^2z-6y^2z
;     dbydx=6xy^2-6yz^2, dbydy=2x^3-6xz^2, dbydz=-12xyz
;     dbzdx=-6x^2z-6y^2z+4z^3, dbzdy=-12xyz, dbzdz=-2x^3-6xy^2+12xz^2
;27. a112=-6,a310=4,a130=-2,b400=1,b202=-3,b022=3,b220=-3,c211=-6,c031=2
;     bx=-6xyz^2+4x^3y-2xy^3, by=x^4-3x^2z^2+3y^2z^2-3x^2y^2, bz=-6x^2yz+2y^3z
;     dbxdx=-6yz^2+12x^2y-2y^3, dbxdy=-6xz^2+4x^3-6xy^2, dbxdz=-12xyz
;     dbydx=4x^3-6xz^2-6xy^2, dbydy=6yz^2-6x^2y, dbydz=-6x^2z+6y^2z
;     dbzdx=-12xyz, dbzdy=-6x^2z+6y^2z, dbzdz=-6x^2y+2y^3
;28. a112=6,a130=-2,b040=1,b202=3,b022=-3,b220=-3,c211=6,c031=-2
;     bx=6xyz^2-2xy^3, by=y^4+3x^2z^2-3y^2z^2-3x^2y^2, bz=6x^2yz-2y^3z
;     dbxdx=6yz^2-2y^3, dbxdy=6xz^2-6xy^2, dxdz=12xyz
;     dbydx=6xz^2-6xy^2, dbydy=4y^3-6yz^2-6x^2y, dbydz=6x^2z-6y^2z
;     dbzdx=12xyz, dbzdy=6x^2z-6y^2z, dbzdz=6x^2y-2y^3
;29. a112=-6,a130=2,b004=1,b202=-3,b022=-3,b220=3,c013=4,c211=-6,c031=-2
;     bx=-6xyz^2+2xy^3, by=z^4-3x^2z^2-3y^2z^2+3x^2y^2, bz=4yz^3-6x^2yz-2y^3z
;     dbxdx=-6yz^2+2y^3, dbxdy=-6xz^2+6xy^2, dbxdz=-12xyz
;     dbydx=-6xz^2+6xy^2, dbydy=-6yz^2+6xy^2, dbydz=4z^3-6x^2z-6y^2z
;     dbzdx=-12xyz, dbzdy=4z^3-6x^2z-6y^2z, dbzdz=12yz^2-6x^2y-2y^3
;30. a103=-2,a121=-6,a301=4,b013=2,b211=-6,c400=1,c220=-3,c202=-3,c022=3
;     bx=-2xz^3-6xy^2z+4x^3z, by=2yz^3-6x^2yz, bz=x^4-3x^2y^2-3x^2z^2+3y^2z^2
;     dbxdx=-2z^3-6y^2z+12x^2z, dbxdy=-12xyz, dbxdz=-6xz^2-6xy^2+4x^3
;     dbydx=-12xyz, dbydy=2z^3-6x^2y, dbydz=6yz^2-6x^2y
;     dbzdx=4x^3-6xy^2-6xz^2, dbzdy=-6x^2y+6yz^2, dbzdz=-6x^2z+6y^2z
;31. a103=2,a121=-6,b013=-2,b031=4,b211=-6,c040=1,c220=-3,c202=3,c022=-3
;     bx=2xz^3-6xy^2z, by=-2yz^3+4y^3z-6x^2yz, bz=y^4-3x^2y^2+3x^2z^2-3y^2z^2
;     dbxdx=2z^3-6y^2z, dbxdy=-12xyz, dbxdz=6xz^2-6xy^2
;     dbydx=-12xyz, dbydy=-2z^3+12y^2z-6x^2z, dbydz=-6yz^2+4y^3-6x^2y
;     dbzdx=-6xy^2+6xz^2, dbzdy=4y^3-6x^2y-6yz^2, dbzdz=6x^2z-6y^2z
;32. a103=-2,a121=6,b013=-2,b211=6,c004=1,c220=3,c202=-3,c022=-3
;     bx=-2xz^3+6xy^2z, by=-2yz^3+6x^2yz, bz=z^4+3x^2y^2-3x^2z^2-3y^2z^2
;     dbxdx=-2z^3+6y^2z, dbxdy=12xyz, dbxdz=-6xz^2+6xy^2
;     dbydx=12xyz, dbydy=-2z^3+6x^2y, dbydz=-6yz^2+6x^2y
;     dbzdx=6xy^2-6xz^2, dbzdy=6x^2y-6yz^2, dbzdz=4z^3-6x^2z-6y^2z
;33. a211=3,a013=-1,b103=-1,b301=1,c112=-3,c310=1
;     bx=3x^2yz-yz^3, by=-xz^3+x^3z, bz=-3xyz^2+x^3y
;     dbxdx=6xyz, dbxdy=3x^2z-z^3, dbxdz=3x^2y-3yz^2
;     dbydx=-z^3+3x^2z, dbydz=-3xz^2+x^3
;     dbzdx=-3yz^2+3x^2y, dbzdy=-3xz^2+x^3, dbzdz=-6xyz
;34. a031=1,a013=-1,b121=3,b103=-1,c112=-3 c130=1
;     bx=y^3z-yz^3, by=3xy^2z-xz^3, bz=-3xyz^2+xy^3
;     dbxdy=3y^2z-z^3, dbxdz=y^3-3yz^2
;     dbydx=3y^2z-z^3, dbydy=6xyz, dbydz=3xy^2-3xz^2
;     dbzdx=-3yz^2+y^3, dbzdy=-3xz^2+3xy^2, dbzdz=-6xyz
    bx(24)=x^4-3*x^2*z^2-3*x^2*y^2+3*y^2*z^2 & by(24)=-2*x^3*y+6*x*y*z^2 & bz(24)=-2*x^3*z+6*x*y^2*z
    bx(25)=y^4+3*x^2*z^2-3*x^2*y^2-3*y^2*z^2 & by(25)=-2*x^3*y+4*x*y^3-6*x*y*z^2 & bz(25)=2*x^3*z-6*x*y^2*z
    bx(26)=z^4-3*x^2*z^2+3*x^2*y^2-3*y^2*z^2 & by(26)=2*x^3*y-6*x*y*z^2 & bz(26)=-2*x^3*z-6*x*y^2*z+4*x*z^3
    bx(27)=-6*x*y*z^2+4*x^3*y-2*x*y^3 & by(27)=x^4-3*x^2*z^2+3*y^2*z^2-3*x^2*y^2 & bz(27)=-6*x^2*y*z+2*y^3*z
    bx(28)=6*x*y*z^2-2*x*y^3 & by(28)=y^4+3*x^2*z^2-3*y^2*z^2-3*x^2*y^2 & bz(28)=6*x^2*y*z-2*y^3*z
    bx(29)=-6*x*y*z^2+2*x*y^3 & by(29)=z^4-3*x^2*z^2-3*y^2*z^2+3*x^2*y^2 & bz(29)=4*y*z^3-6*x^2*y*z-2*y^3*z
    bx(30)=-2*x*z^3-6*x*y^2*z+4*x^3*z & by(30)=2*y*z^3-6*x^2*y*z & bz(30)=x^4-3*x^2*y^2-3*x^2*z^2+3*y^2*z^2
    bx(31)=2*x*z^3-6*x*y^2*z & by(31)=-2*y*z^3+4*y^3*z-6*x^2*y*z & bz(31)=y^4-3*x^2*y^2+3*x^2*z^2-3*y^2*z^2
    bx(32)=-2*x*z^3+6*x*y^2*z & by(32)=-2*y*z^3+6*x^2*y*z & bz(32)=z^4+3*x^2*y^2-3*x^2*z^2-3*y^2*z^2
    bx(33)=3*x^2*y*z-y*z^3 & by(33)=-x*z^3+x^3*z & bz(33)=-3*x*y*z^2+x^3*y
    bx(34)=y^3*z-y*z^3 & by(34)=3*x*y^2*z-x*z^3 & bz(34)=-3*x*y*z^2+x*y^3
    if keyword_set(derivs) then begin
        dbxdx(24)=4*x^3-6*x*y^2-6*x*z^2 & dbxdy(24)=-6*x^2*y+6*y*z^2 & dbxdz(24)=-6*x^2*z+6*y^2*z
        dbydx(24)=-6*x^2*y+6*y*z^2 & dbydy(24)=-2*x^3+6*x*z^2 & dbydz(24)=12*x*y*z
        dbzdx(24)=-6*x^2*z+6*y^2*z & dbzdy(24)=12*x*y*z & dbzdz(24)=-2*x^3+6*x*y^2
        dbxdx(25)=6*x*z^2-6*x*y^2 & dbxdy(25)=4*y^3-6*x^2*y-6*y*z^2 & dbxdz(25)=6*x^2*z-6*y^2*z
        dbydx(25)=-6*x^2*y+4*y^3-6*y*z^2 & dbydy(25)=-2*x^3+12*x*y^2-6*x*z^2 & dbydz(25)=-12*x*y*z
        dbzdx(25)=6*x^2*z-6*y^2*z & dbzdy(25)=-12*x*y*z & dbzdz(25)=2*x^3-6*x*y^2
        dbxdx(26)=-6*x*z^2+6*x*y^2 & dbxdy(26)=6*x^2*y-6*y*z^2 & dbxdz(26)=4*z^3-6*x^2*z-6*y^2*z
        dbydx(26)=6*x^2*y-6*y*z^2 & dbydy(26)=2*x^3-6*x*z^2 & dbydz(26)=-12*x*y*z
        dbzdx(26)=-6*x^2*z-6*y^2*z+4*z^3 & dbzdy(26)=-12*x*y*z & dbzdz(26)=-2*x^3-6*x*y^2+12*x*z^2
        dbxdx(27)=-6*y*z^2+12*x^2*y-2*y^3 & dbxdy(27)=-6*x*z^2+4*x^3-6*x*y^2 & dbxdz(27)=-12*x*y*z
        dbydx(27)=4*x^3-6*x*z^2-6*x*y^2 & dbydy(27)=6*y*z^2-6*x^2*y & dbydz(27)=-6*x^2*z+6*y^2*z
        dbzdx(27)=-12*x*y*z & dbzdy(27)=-6*x^2*z+6*y^2*z & dbzdz(27)=-6*x^2*y+2*y^3
        dbxdx(28)=6*y*z^2-2*y^3 & dbxdy(28)=6*x*z^2-6*x*y^2 & dbxdz(28)=12*x*y*z
        dbydx(28)=6*x*z^2-6*x*y^2 & dbydy(28)=4*y^3-6*y*z^2-6*x^2*y & dbydz(28)=6*x^2*z-6*y^2*z
        dbzdx(28)=12*x*y*z & dbzdy(28)=6*x^2*z-6*y^2*z & dbzdz(28)=6*x^2*y-2*y^3
        dbxdx(29)=-6*y*z^2+2*y^3 & dbxdy(29)=-6*x*z^2+6*x*y^2 & dbxdz(29)=-12*x*y*z
        dbydx(29)=-6*x*z^2+6*x*y^2 & dbydy(29)=-6*y*z^2+6*x^2*y & dbydz(29)=4*z^3-6*x^2*z-6*y^2*z
        dbzdx(29)=-12*x*y*z & dbzdy(29)=4*z^3-6*x^2*z-6*y^2*z & dbzdz(29)=12*y*z^2-6*x^2*y-2*y^3
        dbxdx(30)=-2*z^3-6*y^2*z+12*x^2*z & dbxdy(30)=-12*x*y*z & dbxdz(30)=-6*x*z^2-6*x*y^2+4*x^3
        dbydx(30)=-12*x*y*z & dbydy(30)=2*z^3-6*x^2*z & dbydz(30)=6*y*z^2-6*x^2*y
        dbzdx(30)=4*x^3-6*x*y^2-6*x*z^2 & dbzdy(30)=-6*x^2*y+6*y*z^2 & dbzdz(30)=-6*x^2*z+6*y^2*z
        dbxdx(31)=2*z^3-6*y^2*z & dbxdy(31)=-12*x*y*z & dbxdz(31)=6*x*z^2-6*x*y^2
        dbydx(31)=-12*x*y*z & dbydy(31)=-2*z^3+12*y^2*z-6*x^2*z & dbydz(31)=-6*y*z^2+4*y^3-6*x^2*y
        dbzdx(31)=-6*x*y^2+6*x*z^2 & dbzdy(31)=4*y^3-6*x^2*y-6*y*z^2 & dbzdz(31)=6*x^2*z-6*y^2*z
        dbxdx(32)=-2*z^3+6*y^2*z & dbxdy(32)=12*x*y*z & dbxdz(32)=-6*x*z^2+6*x*y^2
        dbydx(32)=12*x*y*z & dbydy(32)=-2*z^3+6*x^2*z & dbydz(32)=-6*y*z^2+6*x^2*y
        dbzdx(32)=6*x*y^2-6*x*z^2 & dbzdy(32)=6*x^2*y-6*y*z^2 & dbzdz(32)=4*z^3-6*x^2*z-6*y^2*z
        dbxdx(33)=6*x*y*z & dbxdy(33)=3*x^2*z-z^3 & dbxdz(33)=3*x^2*y-3*y*z^2
        dbydx(33)=-z^3+3*x^2*z & dbydz(33)=-3*x*z^2+x^3
        dbzdx(33)=-3*y*z^2+3*x^2*y & dbzdy(33)=-3*x*z^2+x^3 & dbzdz(33)=-6*x*y*z
        dbxdy(34)=3*y^2*z-z^3 & dbxdz(34)=y^3-3*y*z^2
        dbydx(34)=3*y^2*z-z^3 & dbydy(34)=6*x*y*z & dbydz(34)=3*x*y^2-3*x*z^2
        dbzdx(34)=-3*y*z^2+y^3 & dbzdy(34)=-3*x*z^2+3*x*y^2 & dbzdz(34)=-6*x*y*z
    endif
endif

if keyword_set(derivs) then return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz,dbydx:dbydx,dbydy:dbydy,dbydz:dbydz,$
                                    dbzdx:dbzdx,dbzdy:dbzdy,dbzdz:dbzdz} else return,{bx:bx,by:by,bz:bz}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_b_elements_vec,p,derivs=derivs,degree=degree

if not keyword_set(degree) then degree = 4

x=double(p(*,0))
y=double(p(*,1))
z=double(p(*,2))

np=n_elements(x)

nelem = [3,8,15,24,35]
nelem = nelem(floor(degree))

bx=dblarr(np,nelem)
by=bx & bz=bx

if keyword_set(derivs) then dbxdx=bx & dbxdy=bx & dbxdz=bx & dbydx=bx & dbydy=bx & dbydz=bx & dbzdx=bx & dbzdy=bx & dbzdz=bx 

; Constant elements:
; 0. a000=1 : bx=1, db*d* = 0
; 1. b000=1 : by=1; db*d* = 0
; 2. c000=1 : bz=1; db*d* = 0
bx(*,0) = 1
by(*,1) = 1
bz(*,2) = 1

if degree ge 1 then begin
; Linear elements
; 3. b001=1, c010=1: by=z, bz=y, dbydz=1, dbzdy=1
; 4. a001=1, c100=1: bx=z, bz=x, dbxdz=1, dbzdx=1
; 5. a010=1, b100=1, bx=y, by=x, dbxdy=1, dbydx=1
; 6. a100=1, c001=-1, bx=x, bz=-z, dbxdx=1, dbzdz=-1
; 7. a100=1, b010=-1, bx=x, by=-y, dbxdx=1, dbydy=-1
    by(*,3)=z & bz(*,3)=y
    bx(*,4)=z & bz(*,4)=x
    bx(*,5)=y & by(*,5)=x
    bx(*,6)=x & bz(*,6)=-z
    bx(*,7)=x & by(*,7)=-y
    if keyword_set(derivs) then begin
        dbydz(*,3) = 1 & dbzdy(*,3)=1
        dbxdz(*,4) = 1 & dbzdx(*,4)=1
        dbxdy(*,5) = 1 & dbydx(*,5)=1
        dbxdx(*,6) = 1 & dbzdz(*,6)=-1
        dbxdx(*,7) = 1 & dbydy(*,7)=-1
    endif
endif

if degree ge 2 then begin
; Quadratic elements
; 8. a011=1, b101=1,c110=1: bx=yz,by=xz,bz=xy, 
;                           dbxdy=z,dbxdz=y,dbydx=z,dbydz=x,dbzdx=y,dbzdy=x
; 9. a020=1,a200=-1,b110=2: bx=y^2-x^2, by=2xy
;                           dbxdx=-2x, dbxdy=2y, dbydx=2y,dbydy=2x
;10. a002=1,a200=-1,c101=2: bx=z^2-x^2,bz=2xz
;                           dbxdx=-2x, dbxdz=2z, dbzdx=2z,dbzdz=2x
;11. b002=1,b020=-1,c011=2: by=z^2-y^2,bz=2yz
;                           dbydy=-2y, dbydz=2z, dbzdy=2z, dbzdz=2y
;12. a110=2,b200=1,b020=-1: bx=2xy, by=x^2-y^2
;                           dbxdx=2y, dbxdy=2x, dbydx=2x, dbydy=-2y
;13. b011=2,c020=1,c002=-1: by=2yz, bz=y^2-z^2
;                           dbydy=2z, dbydz=2y, dbzdy=2y, dbzdz=-2z
;14. a101=2,c200=1,c002=-1: bx=2xz, bz=x^2-z^2
;                           dbxdx=2z, dbxdz=2x, dbzdx=2x, dbzdz=-2z
;    bx(*,)= & by(*,)= & bz(*,)=
    bx(*,8)=y*z & by(*,8)=x*z & bz(*,8)=x*y
    bx(*,9)= y^2-x^2 & by(*,9)= 2*x*y
    bx(*,10)= z^2-x^2 & bz(*,10)=2*x*z
    by(*,11)= z^2-y^2 & bz(*,11)= 2*y*z
    bx(*,12)= 2*x*y & by(*,12)= x^2-y^2
    by(*,13)= 2*y*z & bz(*,13)= y^2-z^2
    bx(*,14)= 2*x*z & bz(*,14)= x^2-z^2
    if keyword_set(derivs) then begin
        dbxdy(*,8)=z & dbxdz(*,8)=y & dbydx(*,8)=z & dbydz(*,8)=x & dbzdx(*,8)=y & dbzdy(*,8)=x
        dbxdx(*,9)=-2*x & dbxdy(*,9)=2*y & dbydx(*,9)=2*y & dbydy(*,9)=2*x
        dbxdx(*,10)=-2*x & dbxdz(*,10)=2*z & dbzdx(*,10)=2*z & dbzdz(*,10)=2*x
        dbydy(*,11)=-2*y & dbydz(*,11)=2*z & dbzdy(*,11)=2*z & dbzdz(*,11)=2*y
        dbxdx(*,12)=2*y & dbxdy(*,12)=2*x & dbydx(*,12)=2*x & dbydy(*,12)=-2*y
        dbydy(*,13)=2*z & dbydz(*,13)=2*y & dbzdy(*,13)=2*y & dbzdz(*,13)=-2*z
        dbxdx(*,14)=2*z & dbxdz(*,14)=2*x & dbzdx(*,14)=2*x & dbzdz(*,14)=-2*z
    endif
endif

if degree ge 3 then begin
; Cubic elements
;15. a111=6,b021=-3,b201=3,c030=-1,c210=3:
;      bx=6xyz, by=-3y^2z+3x^2z, bz=-y^3+3x^2y
;      dbxdx=6yz, dbxdy=6xz, dbxdz=6xy, dbydx=6xz, dbydy=-6xz,
;      dbydz=-3y^2+3x^2, dbzdx=6xy, dbzdy=-3y2+3x^2
;16. a111=6,b003=-1,b201=3,c210=3,c012=-3
;      bx=6xyz, by=-z^3+3x^2z, bz=3x^2y-3yz^2
;      dbxdx=6yz, dbxdy=6xz, dbxdz=6xy, dbybx=6xz, dbydz=-3z^2+3x^2
;      dbzdx=6xy dbzdy=3x^2-3z^2, dbzdz=-6yz
;17. a210=3,a012=-3,b300=1,b102=-3,c111=-6
;      bx=3x^2*y-3yz^2,by=x^3-3xz^2,bz=-6xyz
;      dbxdx=6xy, dbxdy=3x^2-3z^2, dbxdz=-6yz, dbydx=3x^2-3z^2
;      dbydz=-6xz, dbzdx=-6yz, dbzdy=-6xz, dbzdz=-6xy
;18. a030=1,a012=-3,b120=3,b102=-3,c111=-6
;      bx=y^3-3yz^2, by=3xy^2-3xz^2,bz=-6xyz
;      dbxdy=3y^2-3z^2,dbxdz=-6yz,dbydx=3y^2-3z^2,
;      dbydy=6xy,dbydz=-6xz, dbzdx=-6yz, dbzdy=-6xz, dbzdz=-6xy
;19. a300=2,a120=-3,a102=-3,b210=-3,b012=3,c201=-3,c021=3
;      bx=2x^3-3xy^2-3xz^2, by=-3x^2y+3yz^2, bz=-3x^2z+3y^2z
;      dbxdx=6x^2-3y^2-3z^2, dbxdy=-6xy, dbxdz=-6xz, dbydx=-6xy
;      dbydy=-3x^2+3z^2, dbydz=6yz, dbzdx=-6xz, dbzdy=6yz,
;      dbzdz=-3x^2+3y^2
;20. a102=3,a120=-3,b030=2,b012=-3,b210=-3,c201=3,c021=-3
;      bx=3xz^2-3xy^2, by=2y^3-3yz^2-3x^2y, bz=3x^2z-3y^2z
;      dbxdx=3z^2-3y^2, dbxdy=-6xy, dbxdz=6xz, dbydx=-6xy,
;      dbydy=6y^2-3z^2-3x^2, dbydz=-6yz, dbzdx=6xz, dbzdy=-6yz,
;      dbzdz=3x^2-3y^2
;21. a120=3,a102=-3,b210=3,b012=-3,c003=2,c201=-3,c021=-3
;      bx=3xy^2-3xz^2, by=3x^2y-3yz^2, bz=2z^3-3x^2z-3y^2z
;      dbxdx=3y^2-3z^2, dbxdy=6xy, dbxdz=-6xz, dbydx=6xy, dbydy=3x^2-3z^2
;      dbydz=-6yz, dbzdx=-6xz, dbzdy=-6yz, dbzdz=6z^2-3x^2-3y^2
;22. a003=1,a021=-3,b111=-6,c102=3,c120=-3
;      bx=z^3-3y^2z, by=-6xyz, bz=3xz^2-3xy^2
;      dbxdy=-6yz, dbxdz=3z^2-3y^2, dbydx=-6yz, dbydy=-6xz,
;      dbydz=-6xy, dbzdx=3z^2-3y^2, dbzdy=-6xy, dbzdz=6xz
;23. a201=-3,a021=3,b111=6,c300=-1,c120=3
;      bx=-3x^2z+3y^2z, by=6xyz, bz=-x^3+3xy^2
;      dbxdx=-6xz, dbxdy=6yz, dbxdz=-3x^2+3y^2, dbydx=6yz
;      dbydy=6xz, dbydz=6xy, dbzdx=-3x^2+3y^2, dbzdy=6xy
    bx(*,15)=6*x*y*z & by(*,15)=-3*y^2*z+3*x^2*z & bz(*,15)=-y^3+3*x^2*y
    bx(*,16)=6*x*y*z & by(*,16)=-z^3+3*x^2*z & bz(*,16)=3*x^2*y-3*y*z^2
    bx(*,17)=3*x^2*y-3*y*z^2 & by(*,17)=x^3-3*x*z^2 & bz(*,17)=-6*x*y*z
    bx(*,18)=y^3-3*y*z^2 & by(*,18)=3*x*y^2-3*x*z^2 & bz(*,18)=-6*x*y*z
    bx(*,19)=2*x^3-3*x*y^2-3*x*z^2 & by(*,19)=-3*x^2*y+3*y*z^2 & bz(*,19)=-3*x^2*z+3*y^2*z
    bx(*,20)=3*x*z^2-3*x*y^2 & by(*,20)=2*y^3-3*y*z^2-3*x^2*y & bz(*,20)=3*x^2*z-3*y^2*z
    bx(*,21)=3*x*y^2-3*x*z^2 & by(*,21)=3*x^2*y-3*y*z^2 & bz(*,21)=2*z^3-3*x^2*z-3*y^2*z
    bx(*,22)=z^3-3*y^2*z & by(*,22)=-6*x*y*z & bz(*,22)=3*x*z^2-3*x*y^2
    bx(*,23)=-3*x^2*z+3*y^2*z & by(*,23)=6*x*y*z & bz(*,23)=-x^3+3*x*y^2
    if keyword_set(derivs) then begin
        dbxdx(*,15)=6*y*z & dbxdy(*,15)=6*x*z & dbxdz(*,15)=6*x*y & dbydx(*,15)=6*x*z & dbydy(*,15)=-6*y*z
        dbydz(*,15)=-3*y^2+3*x^2 & dbzdx(*,15)=6*x*y & dbzdy(*,15)=-3*y^2+3*x^2
        dbxdx(*,16)=6*y*z & dbxdy(*,16)=6*x*z & dbxdz(*,16)=6*x*y & dbydx(*,16)=6*x*z & dbydz(*,16)=-3*z^2+3*x^2
        dbzdx(*,16)=6*x*y & dbzdy(*,16)=3*x^2-3*z^2 & dbzdz(*,16)=-6*y*z
        dbxdx(*,17)=6*x*y & dbxdy(*,17)=3*x^2-3*z^2 & dbxdz(*,17)=-6*y*z & dbydx(*,17)=3*x^2-3*z^2
        dbydz(*,17)=-6*x*z & dbzdx(*,17)=-6*y*z & dbzdy(*,17)=-6*x*z & dbzdz(*,17)=-6*x*y
        dbxdy(*,18)=3*y^2-3*z^2 & dbxdz(*,18)=-6*y*z & dbydx(*,18)=3*y^2-3*z^2 
        dbydy(*,18)=6*x*y &dbydz(*,18)=-6*x*z & dbzdx(*,18)=-6*y*z & dbzdy(*,18)=-6*x*z & dbzdz(*,18)=-6*x*y
        dbxdx(*,19)=6*x^2-3*y^2-3*z^2 & dbxdy(*,19)=-6*x*y & dbxdz(*,19)=-6*x*z & dbydx(*,19)=-6*x*y
        dbydy(*,19)=-3*x^2+3*z^2 & dbydz(*,19)=6*y*z & dbzdx(*,19)=-6*x*z & dbzdy(*,19)=6*y*z & dbzdz(*,19)=-3*x^2+3*y^2
        dbxdx(*,20)=3*z^2-3*y^2 & dbxdy(*,20)=-6*x*y & dbxdz(*,20)=6*x*z & dbydx(*,20)=-6*x*y
        dbydy(*,20)=6*y^2-3*z^2-3*x^2 & dbydz(*,20)=-6*y*z & dbzdx(*,20)=6*x*z & dbzdy(*,20)=-6*y*z & dbzdz(*,20)=3*x^2-3*y^2
        dbxdx(*,21)=3*y^2-3*z^2 & dbxdy(*,21)=6*x*y & dbxdz(*,21)=-6*x*z & dbydx(*,21)=6*x*y & dbydy(*,21)=3*x^2-3*z^2
        dbydz(*,21)=-6*y*z & dbzdx(*,21)=-6*x*z & dbzdy(*,21)=-6*y*z & dbzdz(*,21)=6*z^2-3*x^2-3*y^2
        dbxdy(*,22)=-6*y*z & dbxdz(*,22)=3*z^2-3*y^2 & dbydx(*,22)=-6*y*z & dbydy(*,22)=-6*x*z 
        dbydz(*,22)=-6*x*y & dbzdx(*,22)=3*z^2-3*y^2 & dbzdy(*,22)=-6*x*y & dbzdz(*,22)=6*x*z
        dbxdx(*,23)=-6*x*z & dbxdy(*,23)=6*y*z & dbxdz(*,23)=-3*x^2+3*y^2 & dbydx(*,23)=6*y*z
        dbydy(*,23)=6*x*z & dbydz(*,23)=6*x*y & dbzdx(*,23)=-3*x^2+3*y^2 & dbzdy(*,23)=6*x*y
    endif 
endif

if degree ge 4 then begin
; Quartic elements
;24. a400=1,a202=-3,a220=-3,a022=3,b310=-2,b112=6,c301=-2,c121=6
;     bx=x^4-3x^2z^2-3x^2y^2+3y^2z^2, by=-2x^3y+6xyz^2, bz=-2x^3z+6xy^2z
;     dbxdx=4x^3-6xy^2-6xz^2, dbxdy=-6x^2y+6yz^2, dbxdz=-6x^2z+6y^2z
;     dbydx=-6x^2y+6yz^2, dbydy=-2x^3+6xz^2, dbydz=12xyz
;     dbzdx=-6x^2z+6y^2z, dbzdy=12xyz, dbzdz=-2x^2+6xy^2
;25. a040=1,a202=3,a220=-3,a022=-3,b310=-2,b130=4,b112=-6,c301=2,c121=-6
;     bx=y^4+3x^2z^2-3x^2y^2-3y^2z^2, by=-2x^3y+4xy^3-6xyz^2
;     bz=2x^3z-6xy^2z
;     dbxdx=6xz^2-6xy^2, dbxdy=4y^3-6x^2y-6yz^2, dbxdz=6x^2z-6y^2z
;     dbydx=-6x^2y+4y^3-6yz^2, dbydy=-2x^3+12xy^2-6xz^2, dbydz=-12xyz
;     dbzdx=6x^2z-6y^2z, dbzdy=-12xyz, dbzdz=2x^3-6xy^2
;26. a004=1,a202=-3,a220=3,a022=-3,b310=2,b112=-6,c301=-2,c121=-6,c103=4
;     bx=z^4-3x^2z^2+3x^2y^2-3y^2z^2, by=2x^3y-6xyz^2, bz=-2x^3z-6xy^2z+4xz^3
;     dbxdx=-6xz^2+6xy^2, dbxdy=6x^2y-6yz^2, dbxdz=4z^3-6x^2z-6y^2z
;     dbydx=6xy^2-6yz^2, dbydy=2x^3-6xz^2, dbydz=-12xyz
;     dbzdx=-6x^2z-6y^2z+4z^3, dbzdy=-12xyz, dbzdz=-2x^3-6xy^2+12xz^2
;27. a112=-6,a310=4,a130=-2,b400=1,b202=-3,b022=3,b220=-3,c211=-6,c031=2
;     bx=-6xyz^2+4x^3y-2xy^3, by=x^4-3x^2z^2+3y^2z^2-3x^2y^2, bz=-6x^2yz+2y^3z
;     dbxdx=-6yz^2+12x^2y-2y^3, dbxdy=-6xz^2+4x^3-6xy^2, dbxdz=-12xyz
;     dbydx=4x^3-6xz^2-6xy^2, dbydy=6yz^2-6x^2y, dbydz=-6x^2z+6y^2z
;     dbzdx=-12xyz, dbzdy=-6x^2z+6y^2z, dbzdz=-6x^2y+2y^3
;28. a112=6,a130=-2,b040=1,b202=3,b022=-3,b220=-3,c211=6,c031=-2
;     bx=6xyz^2-2xy^3, by=y^4+3x^2z^2-3y^2z^2-3x^2y^2, bz=6x^2yz-2y^3z
;     dbxdx=6yz^2-2y^3, dbxdy=6xz^2-6xy^2, dxdz=12xyz
;     dbydx=6xz^2-6xy^2, dbydy=4y^3-6yz^2-6x^2y, dbydz=6x^2z-6y^2z
;     dbzdx=12xyz, dbzdy=6x^2z-6y^2z, dbzdz=6x^2y-2y^3
;29. a112=-6,a130=2,b004=1,b202=-3,b022=-3,b220=3,c013=4,c211=-6,c031=-2
;     bx=-6xyz^2+2xy^3, by=z^4-3x^2z^2-3y^2z^2+3x^2y^2, bz=4yz^3-6x^2yz-2y^3z
;     dbxdx=-6yz^2+2y^3, dbxdy=-6xz^2+6xy^2, dbxdz=-12xyz
;     dbydx=-6xz^2+6xy^2, dbydy=-6yz^2+6xy^2, dbydz=4z^3-6x^2z-6y^2z
;     dbzdx=-12xyz, dbzdy=4z^3-6x^2z-6y^2z, dbzdz=12yz^2-6x^2y-2y^3
;30. a103=-2,a121=-6,a301=4,b013=2,b211=-6,c400=1,c220=-3,c202=-3,c022=3
;     bx=-2xz^3-6xy^2z+4x^3z, by=2yz^3-6x^2yz, bz=x^4-3x^2y^2-3x^2z^2+3y^2z^2
;     dbxdx=-2z^3-6y^2z+12x^2z, dbxdy=-12xyz, dbxdz=-6xz^2-6xy^2+4x^3
;     dbydx=-12xyz, dbydy=2z^3-6x^2y, dbydz=6yz^2-6x^2y
;     dbzdx=4x^3-6xy^2-6xz^2, dbzdy=-6x^2y+6yz^2, dbzdz=-6x^2z+6y^2z
;31. a103=2,a121=-6,b013=-2,b031=4,b211=-6,c040=1,c220=-3,c202=3,c022=-3
;     bx=2xz^3-6xy^2z, by=-2yz^3+4y^3z-6x^2yz, bz=y^4-3x^2y^2+3x^2z^2-3y^2z^2
;     dbxdx=2z^3-6y^2z, dbxdy=-12xyz, dbxdz=6xz^2-6xy^2
;     dbydx=-12xyz, dbydy=-2z^3+12y^2z-6x^2z, dbydz=-6yz^2+4y^3-6x^2y
;     dbzdx=-6xy^2+6xz^2, dbzdy=4y^3-6x^2y-6yz^2, dbzdz=6x^2z-6y^2z
;32. a103=-2,a121=6,b013=-2,b211=6,c004=1,c220=3,c202=-3,c022=-3
;     bx=-2xz^3+6xy^2z, by=-2yz^3+6x^2yz, bz=z^4+3x^2y^2-3x^2z^2-3y^2z^2
;     dbxdx=-2z^3+6y^2z, dbxdy=12xyz, dbxdz=-6xz^2+6xy^2
;     dbydx=12xyz, dbydy=-2z^3+6x^2y, dbydz=-6yz^2+6x^2y
;     dbzdx=6xy^2-6xz^2, dbzdy=6x^2y-6yz^2, dbzdz=4z^3-6x^2z-6y^2z
;33. a211=3,a013=-1,b103=-1,b301=1,c112=-3,c310=1
;     bx=3x^2yz-yz^3, by=-xz^3+x^3z, bz=-3xyz^2+x^3y
;     dbxdx=6xyz, dbxdy=3x^2z-z^3, dbxdz=3x^2y-3yz^2
;     dbydx=-z^3+3x^2z, dbydz=-3xz^2+x^3
;     dbzdx=-3yz^2+3x^2y, dbzdy=-3xz^2+x^3, dbzdz=-6xyz
;34. a031=1,a013=-1,b121=3,b103=-1,c112=-3 c130=1
;     bx=y^3z-yz^3, by=3xy^2z-xz^3, bz=-3xyz^2+xy^3
;     dbxdy=3y^2z-z^3, dbxdz=y^3-3yz^2
;     dbydx=3y^2z-z^3, dbydy=6xyz, dbydz=3xy^2-3xz^2
;     dbzdx=-3yz^2+y^3, dbzdy=-3xz^2+3xy^2, dbzdz=-6xyz
    bx(*,24)=x^4-3*x^2*z^2-3*x^2*y^2+3*y^2*z^2 & by(*,24)=-2*x^3*y+6*x*y*z^2 & bz(*,24)=-2*x^3*z+6*x*y^2*z
    bx(*,25)=y^4+3*x^2*z^2-3*x^2*y^2-3*y^2*z^2 & by(*,25)=-2*x^3*y+4*x*y^3-6*x*y*z^2 & bz(*,25)=2*x^3*z-6*x*y^2*z
    bx(*,26)=z^4-3*x^2*z^2+3*x^2*y^2-3*y^2*z^2 & by(*,26)=2*x^3*y-6*x*y*z^2 & bz(*,26)=-2*x^3*z-6*x*y^2*z+4*x*z^3
    bx(*,27)=-6*x*y*z^2+4*x^3*y-2*x*y^3 & by(*,27)=x^4-3*x^2*z^2+3*y^2*z^2-3*x^2*y^2 & bz(*,27)=-6*x^2*y*z+2*y^3*z
    bx(*,28)=6*x*y*z^2-2*x*y^3 & by(*,28)=y^4+3*x^2*z^2-3*y^2*z^2-3*x^2*y^2 & bz(*,28)=6*x^2*y*z-2*y^3*z
    bx(*,29)=-6*x*y*z^2+2*x*y^3 & by(*,29)=z^4-3*x^2*z^2-3*y^2*z^2+3*x^2*y^2 & bz(*,29)=4*y*z^3-6*x^2*y*z-2*y^3*z
    bx(*,30)=-2*x*z^3-6*x*y^2*z+4*x^3*z & by(*,30)=2*y*z^3-6*x^2*y*z & bz(*,30)=x^4-3*x^2*y^2-3*x^2*z^2+3*y^2*z^2
    bx(*,31)=2*x*z^3-6*x*y^2*z & by(*,31)=-2*y*z^3+4*y^3*z-6*x^2*y*z & bz(*,31)=y^4-3*x^2*y^2+3*x^2*z^2-3*y^2*z^2
    bx(*,32)=-2*x*z^3+6*x*y^2*z & by(*,32)=-2*y*z^3+6*x^2*y*z & bz(*,32)=z^4+3*x^2*y^2-3*x^2*z^2-3*y^2*z^2
    bx(*,33)=3*x^2*y*z-y*z^3 & by(*,33)=-x*z^3+x^3*z & bz(*,33)=-3*x*y*z^2+x^3*y
    bx(*,34)=y^3*z-y*z^3 & by(*,34)=3*x*y^2*z-x*z^3 & bz(*,34)=-3*x*y*z^2+x*y^3
    if keyword_set(derivs) then begin
        dbxdx(*,24)=4*x^3-6*x*y^2-6*x*z^2 & dbxdy(*,24)=-6*x^2*y+6*y*z^2 & dbxdz(*,24)=-6*x^2*z+6*y^2*z
        dbydx(*,24)=-6*x^2*y+6*y*z^2 & dbydy(*,24)=-2*x^3+6*x*z^2 & dbydz(*,24)=12*x*y*z
        dbzdx(*,24)=-6*x^2*z+6*y^2*z & dbzdy(*,24)=12*x*y*z & dbzdz(*,24)=-2*x^3+6*x*y^2
        dbxdx(*,25)=6*x*z^2-6*x*y^2 & dbxdy(*,25)=4*y^3-6*x^2*y-6*y*z^2 & dbxdz(*,25)=6*x^2*z-6*y^2*z
        dbydx(*,25)=-6*x^2*y+4*y^3-6*y*z^2 & dbydy(*,25)=-2*x^3+12*x*y^2-6*x*z^2 & dbydz(*,25)=-12*x*y*z
        dbzdx(*,25)=6*x^2*z-6*y^2*z & dbzdy(*,25)=-12*x*y*z & dbzdz(*,25)=2*x^3-6*x*y^2
        dbxdx(*,26)=-6*x*z^2+6*x*y^2 & dbxdy(*,26)=6*x^2*y-6*y*z^2 & dbxdz(*,26)=4*z^3-6*x^2*z-6*y^2*z
        dbydx(*,26)=6*x^2*y-6*y*z^2 & dbydy(*,26)=2*x^3-6*x*z^2 & dbydz(*,26)=-12*x*y*z
        dbzdx(*,26)=-6*x^2*z-6*y^2*z+4*z^3 & dbzdy(*,26)=-12*x*y*z & dbzdz(*,26)=-2*x^3-6*x*y^2+12*x*z^2
        dbxdx(*,27)=-6*y*z^2+12*x^2*y-2*y^3 & dbxdy(*,27)=-6*x*z^2+4*x^3-6*x*y^2 & dbxdz(*,27)=-12*x*y*z
        dbydx(*,27)=4*x^3-6*x*z^2-6*x*y^2 & dbydy(*,27)=6*y*z^2-6*x^2*y & dbydz(*,27)=-6*x^2*z+6*y^2*z
        dbzdx(*,27)=-12*x*y*z & dbzdy(*,27)=-6*x^2*z+6*y^2*z & dbzdz(*,27)=-6*x^2*y+2*y^3
        dbxdx(*,28)=6*y*z^2-2*y^3 & dbxdy(*,28)=6*x*z^2-6*x*y^2 & dbxdz(*,28)=12*x*y*z
        dbydx(*,28)=6*x*z^2-6*x*y^2 & dbydy(*,28)=4*y^3-6*y*z^2-6*x^2*y & dbydz(*,28)=6*x^2*z-6*y^2*z
        dbzdx(*,28)=12*x*y*z & dbzdy(*,28)=6*x^2*z-6*y^2*z & dbzdz(*,28)=6*x^2*y-2*y^3
        dbxdx(*,29)=-6*y*z^2+2*y^3 & dbxdy(*,29)=-6*x*z^2+6*x*y^2 & dbxdz(*,29)=-12*x*y*z
        dbydx(*,29)=-6*x*z^2+6*x*y^2 & dbydy(*,29)=-6*y*z^2+6*x^2*y & dbydz(*,29)=4*z^3-6*x^2*z-6*y^2*z
        dbzdx(*,29)=-12*x*y*z & dbzdy(*,29)=4*z^3-6*x^2*z-6*y^2*z & dbzdz(*,29)=12*y*z^2-6*x^2*y-2*y^3
        dbxdx(*,30)=-2*z^3-6*y^2*z+12*x^2*z & dbxdy(*,30)=-12*x*y*z & dbxdz(*,30)=-6*x*z^2-6*x*y^2+4*x^3
        dbydx(*,30)=-12*x*y*z & dbydy(*,30)=2*z^3-6*x^2*z & dbydz(*,30)=6*y*z^2-6*x^2*y
        dbzdx(*,30)=4*x^3-6*x*y^2-6*x*z^2 & dbzdy(*,30)=-6*x^2*y+6*y*z^2 & dbzdz(*,30)=-6*x^2*z+6*y^2*z
        dbxdx(*,31)=2*z^3-6*y^2*z & dbxdy(*,31)=-12*x*y*z & dbxdz(*,31)=6*x*z^2-6*x*y^2
        dbydx(*,31)=-12*x*y*z & dbydy(*,31)=-2*z^3+12*y^2*z-6*x^2*z & dbydz(*,31)=-6*y*z^2+4*y^3-6*x^2*y
        dbzdx(*,31)=-6*x*y^2+6*x*z^2 & dbzdy(*,31)=4*y^3-6*x^2*y-6*y*z^2 & dbzdz(*,31)=6*x^2*z-6*y^2*z
        dbxdx(*,32)=-2*z^3+6*y^2*z & dbxdy(*,32)=12*x*y*z & dbxdz(*,32)=-6*x*z^2+6*x*y^2
        dbydx(*,32)=12*x*y*z & dbydy(*,32)=-2*z^3+6*x^2*z & dbydz(*,32)=-6*y*z^2+6*x^2*y
        dbzdx(*,32)=6*x*y^2-6*x*z^2 & dbzdy(*,32)=6*x^2*y-6*y*z^2 & dbzdz(*,32)=4*z^3-6*x^2*z-6*y^2*z
        dbxdx(*,33)=6*x*y*z & dbxdy(*,33)=3*x^2*z-z^3 & dbxdz(*,33)=3*x^2*y-3*y*z^2
        dbydx(*,33)=-z^3+3*x^2*z & dbydz(*,33)=-3*x*z^2+x^3
        dbzdx(*,33)=-3*y*z^2+3*x^2*y & dbzdy(*,33)=-3*x*z^2+x^3 & dbzdz(*,33)=-6*x*y*z
        dbxdy(*,34)=3*y^2*z-z^3 & dbxdz(*,34)=y^3-3*y*z^2
        dbydx(*,34)=3*y^2*z-z^3 & dbydy(*,34)=6*x*y*z & dbydz(*,34)=3*x*y^2-3*x*z^2
        dbzdx(*,34)=-3*y*z^2+y^3 & dbzdy(*,34)=-3*x*z^2+3*x*y^2 & dbzdz(*,34)=-6*x*y*z
    endif
endif

if keyword_set(derivs) then return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz,dbydx:dbydx,dbydy:dbydy,dbydz:dbydz,$
                                    dbzdx:dbzdx,dbzdy:dbzdy,dbzdz:dbzdz} else return,{bx:bx,by:by,bz:bz}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_b_elements_timesf,p,f,derivs=derivs,degree=degree

if not keyword_set(degree) then degree = 4

x=double(p(*,0))
y=double(p(*,1))
z=double(p(*,2))

np=n_elements(x)

nelem = [3,8,15,24,35]
nelem = nelem(floor(degree))

bx=dblarr(np)
by=bx & bz=bx

if keyword_set(derivs) then dbxdx=bx & dbxdy=bx & dbxdz=bx & dbydx=bx & dbydy=bx & dbydz=bx & dbzdx=bx & dbzdy=bx & dbzdz=bx 

; Constant elements:
; 0. a000=1 : bx=1, db*d* = 0
; 1. b000=1 : by=1; db*d* = 0
; 2. c000=1 : bz=1; db*d* = 0
bx += f(*,0)
by += f(*,1)
bz += f(*,2)

if degree ge 1 then begin
; Linear elements
; 3. b001=1, c010=1: by=z, bz=y, dbydz=1, dbzdy=1
; 4. a001=1, c100=1: bx=z, bz=x, dbxdz=1, dbzdx=1
; 5. a010=1, b100=1, bx=y, by=x, dbxdy=1, dbydx=1
; 6. a100=1, c001=-1, bx=x, bz=-z, dbxdx=1, dbzdz=-1
; 7. a100=1, b010=-1, bx=x, by=-y, dbxdx=1, dbydy=-1
    by+=z*f(*,3) & bz+=y*f(*,3)
    bx+=z*f(*,4) & bz+=x*f(*,4)
    bx+=y*f(*,5) & by+=x*f(*,5)
    bx+=x*f(*,6) & bz+=-z*f(*,6)
    bx+=x*f(*,7) & by+=-y*f(*,7)
    if keyword_set(derivs) then begin
        dbydz += f(*,3) & dbzdy+=f(*,3)
        dbxdz += f(*,4) & dbzdx+=f(*,4)
        dbxdy += f(*,5) & dbydx+=f(*,5)
        dbxdx += f(*,6) & dbzdz+=-f(*,6)
        dbxdx += f(*,7) & dbydy+=-f(*,7)
    endif
endif

if degree ge 2 then begin
; Quadratic elements
; 8. a011=1, b101=1,c110=1: bx=yz,by=xz,bz=xy, 
;                           dbxdy=z,dbxdz=y,dbydx=z,dbydz=x,dbzdx=y,dbzdy=x
; 9. a020=1,a200=-1,b110=2: bx=y^2-x^2, by=2xy
;                           dbxdx=-2x, dbxdy=2y, dbydx=2y,dbydy=2x
;10. a002=1,a200=-1,c101=2: bx=z^2-x^2,bz=2xz
;                           dbxdx=-2x, dbxdz=2z, dbzdx=2z,dbzdz=2x
;11. b002=1,b020=-1,c011=2: by=z^2-y^2,bz=2yz
;                           dbydy=-2y, dbydz=2z, dbzdy=2z, dbzdz=2y
;12. a110=2,b200=1,b020=-1: bx=2xy, by=x^2-y^2
;                           dbxdx=2y, dbxdy=2x, dbydx=2x, dbydy=-2y
;13. b011=2,c020=1,c002=-1: by=2yz, bz=y^2-z^2
;                           dbydy=2z, dbydz=2y, dbzdy=2y, dbzdz=-2z
;14. a101=2,c200=1,c002=-1: bx=2xz, bz=x^2-z^2
;                           dbxdx=2z, dbxdz=2x, dbzdx=2x, dbzdz=-2z
    bx+=y*z*f(*,8) & by+=x*z*f(*,8) & bz+=x*y*f(*,8)
    bx+= (y^2-x^2)*f(*,9) & by+= 2*x*y*f(*,9)
    bx+= (z^2-x^2)*f(*,10) & bz+=2*x*z*f(*,10)
    by+= (z^2-y^2)*f(*,11) & bz+= 2*y*z*f(*,11)
    bx+= (2*x*y)*f(*,12) & by+= (x^2-y^2)*f(*,12)
    by+= 2*y*z*f(*,13) & bz+= (y^2-z^2)*f(*,13)
    bx+= 2*x*z*f(*,14) & bz+= (x^2-z^2)*f(*,14)
    if keyword_set(derivs) then begin
        dbxdy+=z*f(*,8) & dbxdz+=y*f(*,8) & dbydx+=z*f(*,8) & dbydz+=x*f(*,8) & dbzdx+=y*f(*,8) & dbzdy+=x*f(*,8)
        dbxdx+=-2*x*f(*,9) & dbxdy+=2*y*f(*,9) & dbydx+=2*y*f(*,9) & dbydy+=2*x*f(*,9)
        dbxdx+=-2*x*f(*,10) & dbxdz+=2*z*f(*,10) & dbzdx+=2*z*f(*,10) & dbzdz+=2*x*f(*,10)
        dbydy+=-2*y*f(*,11) & dbydz+=2*z*f(*,11) & dbzdy+=2*z*f(*,11) & dbzdz+=2*y*f(*,11)
        dbxdx+=2*y*f(*,12) & dbxdy+=2*x*f(*,12) & dbydx+=2*x*f(*,12) & dbydy+=-2*y*f(*,12)
        dbydy+=2*z*f(*,13) & dbydz+=2*y*f(*,13) & dbzdy+=2*y*f(*,13) & dbzdz+=-2*z*f(*,13)
        dbxdx+=2*z*f(*,14) & dbxdz+=2*x*f(*,14) & dbzdx+=2*x*f(*,14) & dbzdz+=-2*z*f(*,14)
    endif
endif

if degree ge 3 then begin
; Cubic elements
;15. a111=6,b021=-3,b201=3,c030=-1,c210=3:
;      bx=6xyz, by=-3y^2z+3x^2z, bz=-y^3+3x^2y
;      dbxdx=6yz, dbxdy=6xz, dbxdz=6xy, dbydx=6xz, dbydy=-6xz,
;      dbydz=-3y^2+3x^2, dbzdx=6xy, dbzdy=-3y2+3x^2
;16. a111=6,b003=-1,b201=3,c210=3,c012=-3
;      bx=6xyz, by=-z^3+3x^2z, bz=3x^2y-3yz^2
;      dbxdx=6yz, dbxdy=6xz, dbxdz=6xy, dbybx=6xz, dbydz=-3z^2+3x^2
;      dbzdx=6xy dbzdy=3x^2-3z^2, dbzdz=-6yz
;17. a210=3,a012=-3,b300=1,b102=-3,c111=-6
;      bx=3x^2*y-3yz^2,by=x^3-3xz^2,bz=-6xyz
;      dbxdx=6xy, dbxdy=3x^2-3z^2, dbxdz=-6yz, dbydx=3x^2-3z^2
;      dbydz=-6xz, dbzdx=-6yz, dbzdy=-6xz, dbzdz=-6xy
;18. a030=1,a012=-3,b120=3,b102=-3,c111=-6
;      bx=y^3-3yz^2, by=3xy^2-3xz^2,bz=-6xyz
;      dbxdy=3y^2-3z^2,dbxdz=-6yz,dbydx=3y^2-3z^2,
;      dbydy=6xy,dbydz=-6xz, dbzdx=-6yz, dbzdy=-6xz, dbzdz=-6xy
;19. a300=2,a120=-3,a102=-3,b210=-3,b012=3,c201=-3,c021=3
;      bx=2x^3-3xy^2-3xz^2, by=-3x^2y+3yz^2, bz=-3x^2z+3y^2z
;      dbxdx=6x^2-3y^2-3z^2, dbxdy=-6xy, dbxdz=-6xz, dbydx=-6xy
;      dbydy=-3x^2+3z^2, dbydz=6yz, dbzdx=-6xz, dbzdy=6yz,
;      dbzdz=-3x^2+3y^2
;20. a102=3,a120=-3,b030=2,b012=-3,b210=-3,c201=3,c021=-3
;      bx=3xz^2-3xy^2, by=2y^3-3yz^2-3x^2y, bz=3x^2z-3y^2z
;      dbxdx=3z^2-3y^2, dbxdy=-6xy, dbxdz=6xz, dbydx=-6xy,
;      dbydy=6y^2-3z^2-3x^2, dbydz=-6yz, dbzdx=6xz, dbzdy=-6yz,
;      dbzdz=3x^2-3y^2
;21. a120=3,a102=-3,b210=3,b012=-3,c003=2,c201=-3,c021=-3
;      bx=3xy^2-3xz^2, by=3x^2y-3yz^2, bz=2z^3-3x^2z-3y^2z
;      dbxdx=3y^2-3z^2, dbxdy=6xy, dbxdz=-6xz, dbydx=6xy, dbydy=3x^2-3z^2
;      dbydz=-6yz, dbzdx=-6xz, dbzdy=-6yz, dbzdz=6z^2-3x^2-3y^2
;22. a003=1,a021=-3,b111=-6,c102=3,c120=-3
;      bx=z^3-3y^2z, by=-6xyz, bz=3xz^2-3xy^2
;      dbxdy=-6yz, dbxdz=3z^2-3y^2, dbydx=-6yz, dbydy=-6xz,
;      dbydz=-6xy, dbzdx=3z^2-3y^2, dbzdy=-6xy, dbzdz=6xz
;23. a201=-3,a021=3,b111=6,c300=-1,c120=3
;      bx=-3x^2z+3y^2z, by=6xyz, bz=-x^3+3xy^2
;      dbxdx=-6xz, dbxdy=6yz, dbxdz=-3x^2+3y^2, dbydx=6yz
;      dbydy=6xz, dbydz=6xy, dbzdx=-3x^2+3y^2, dbzdy=6xy
    bx+=6*x*y*z*f(*,15) & by+=(-3*y^2*z+3*x^2*z)*f(*,15) & bz+=(-y^3+3*x^2*y)*f(*,15)
    bx+=6*x*y*z*f(*,16) & by+=(-z^3+3*x^2*z)*f(*,16) & bz+=(3*x^2*y-3*y*z^2)*f(*,16)
    bx+=(3*x^2*y-3*y*z^2)*f(*,17) & by+=(x^3-3*x*z^2)*f(*,17) & bz+=-6*x*y*z*f(*,17)
    bx+=(y^3-3*y*z^2)*f(*,18) & by+=(3*x*y^2-3*x*z^2)*f(*,18) & bz+=-6*x*y*z*f(*,18)
    bx+=(2*x^3-3*x*y^2-3*x*z^2)*f(*,19) & by+=(-3*x^2*y+3*y*z^2)*f(*,19) & bz+=(-3*x^2*z+3*y^2*z)*f(*,19)
    bx+=(3*x*z^2-3*x*y^2)*f(*,20) & by+=(2*y^3-3*y*z^2-3*x^2*y)*f(*,20) & bz+=(3*x^2*z-3*y^2*z)*f(*,20)
    bx+=(3*x*y^2-3*x*z^2)*f(*,21) & by+=(3*x^2*y-3*y*z^2)*f(*,21) & bz+=(2*z^3-3*x^2*z-3*y^2*z)*f(*,21)
    bx+=(z^3-3*y^2*z)*f(*,22) & by+=-6*x*y*z*f(*,22) & bz+=(3*x*z^2-3*x*y^2)*f(*,22)
    bx+=(-3*x^2*z+3*y^2*z)*f(*,23) & by+=6*x*y*z*f(*,23) & bz+=(-x^3+3*x*y^2)*f(*,23)
    if keyword_set(derivs) then begin
        dbxdx+=6*y*z*f(*,15) & dbxdy+=6*x*z*f(*,15) & dbxdz+=6*x*y*f(*,15) & dbydx+=6*x*z*f(*,15) & dbydy+=-6*y*z*f(*,15)
        dbydz+=(-3*y^2+3*x^2)*f(*,15) & dbzdx+=6*x*y*f(*,15) & dbzdy+=(-3*y^2+3*x^2)*f(*,15)
        dbxdx+=6*y*z*f(*,16) & dbxdy+=6*x*z*f(*,16) & dbxdz+=6*x*y*f(*,16) & dbydx+=6*x*z*f(*,16) & dbydz+=(-3*z^2+3*x^2)*f(*,16)
        dbzdx+=6*x*y*f(*,16) & dbzdy+=(3*x^2-3*z^2)*f(*,16) & dbzdz+=-6*y*z*f(*,16)
        dbxdx+=6*x*y*f(*,17) & dbxdy+=(3*x^2-3*z^2)*f(*,17) & dbxdz+=-6*y*z*f(*,17) & dbydx+=(3*x^2-3*z^2)*f(*,17)
        dbydz+=-6*x*z*f(*,17) & dbzdx+=-6*y*z*f(*,17) & dbzdy+=-6*x*z*f(*,17) & dbzdz+=-6*x*y*f(*,17)
        dbxdy+=(3*y^2-3*z^2)*f(*,18) & dbxdz+=-6*y*z*f(*,18) & dbydx+=(3*y^2-3*z^2)*f(*,18)
        dbydy+=6*x*y*f(*,18) &dbydz+=-6*x*z*f(*,18) & dbzdx+=-6*y*z*f(*,18) & dbzdy+=-6*x*z*f(*,18) & dbzdz+=-6*x*y*f(*,18)
        dbxdx+=(6*x^2-3*y^2-3*z^2)*f(*,19) & dbxdy+=-6*x*y*f(*,19) & dbxdz+=-6*x*z*f(*,19) & dbydx+=-6*x*y*f(*,19)
        dbydy+=(-3*x^2+3*z^2)*f(*,19) & dbydz+=6*y*z*f(*,19) & dbzdx+=-6*x*z*f(*,19) & dbzdy+=6*y*z*f(*,19) & dbzdz+=(-3*x^2+3*y^2)*f(*,19)
        dbxdx+=(3*z^2-3*y^2)*f(*,20) & dbxdy+=-6*x*y*f(*,20) & dbxdz+=6*x*z*f(*,20) & dbydx+=-6*x*y*f(*,20)
        dbydy+=(6*y^2-3*z^2-3*x^2)*f(*,20) & dbydz+=-6*y*z*f(*,20) & dbzdx+=6*x*z*f(*,20) & dbzdy+=-6*y*z*f(*,20) & dbzdz+=(3*x^2-3*y^2)*f(*,20)
        dbxdx+=(3*y^2-3*z^2)*f(*,21) & dbxdy+=6*x*y*f(*,21) & dbxdz+=-6*x*z*f(*,21) & dbydx+=6*x*y*f(*,21) & dbydy+=(3*x^2-3*z^2)*f(*,21)
        dbydz+=-6*y*z*f(*,21) & dbzdx+=-6*x*z*f(*,21) & dbzdy+=-6*y*z*f(*,21) & dbzdz+=(6*z^2-3*x^2-3*y^2)*f(*,21)
        dbxdy+=-6*y*z*f(*,22) & dbxdz+=(3*z^2-3*y^2)*f(*,22) & dbydx+=-6*y*z*f(*,22) & dbydy+=-6*x*z*f(*,22) 
        dbydz+=-6*x*y*f(*,22) & dbzdx+=(3*z^2-3*y^2)*f(*,22) & dbzdy+=-6*x*y*f(*,22) & dbzdz+=6*x*z*f(*,22)
        dbxdx+=-6*x*z*f(*,23) & dbxdy+=6*y*z*f(*,23) & dbxdz+=(-3*x^2+3*y^2)*f(*,23) & dbydx+=6*y*z*f(*,23)
        dbydy+=6*x*z*f(*,23) & dbydz+=6*x*y*f(*,23) & dbzdx+=(-3*x^2+3*y^2)*f(*,23) & dbzdy+=6*x*y*f(*,23)
    endif 
endif

if degree ge 4 then begin
; Quartic elements
;24. a400=1,a202=-3,a220=-3,a022=3,b310=-2,b112=6,c301=-2,c121=6
;     bx=x^4-3x^2z^2-3x^2y^2+3y^2z^2, by=-2x^3y+6xyz^2, bz=-2x^3z+6xy^2z
;     dbxdx=4x^3-6xy^2-6xz^2, dbxdy=-6x^2y+6yz^2, dbxdz=-6x^2z+6y^2z
;     dbydx=-6x^2y+6yz^2, dbydy=-2x^3+6xz^2, dbydz=12xyz
;     dbzdx=-6x^2z+6y^2z, dbzdy=12xyz, dbzdz=-2x^2+6xy^2
;25. a040=1,a202=3,a220=-3,a022=-3,b310=-2,b130=4,b112=-6,c301=2,c121=-6
;     bx=y^4+3x^2z^2-3x^2y^2-3y^2z^2, by=-2x^3y+4xy^3-6xyz^2, bz=2x^3z-6xy^2z
;     dbxdx=6xz^2-6xy^2, dbxdy=4y^3-6x^2y-6yz^2, dbxdz=6x^2z-6y^2z
;     dbydx=-6x^2y+4y^3-6yz^2, dbydy=-2x^3+12xy^2-6xz^2, dbydz=-12xyz
;     dbzdx=6x^2z-6y^2z, dbzdy=-12xyz, dbzdz=2x^3-6xy^2
;26. a004=1,a202=-3,a220=3,a022=-3,b310=2,b112=-6,c301=-2,c121=-6,c103=4
;     bx=z^4-3x^2z^2+3x^2y^2-3y^2z^2, by=2x^3y-6xyz^2, bz=-2x^3z-6xy^2z+4xz^3
;     dbxdx=-6xz^2+6xy^2, dbxdy=6x^2y-6yz^2, dbxdz=4z^3-6x^2z-6y^2z
;     dbydx=6xy^2-6yz^2, dbydy=2x^3-6xz^2, dbydz=-12xyz
;     dbzdx=-6x^2z-6y^2z+4z^3, dbzdy=-12xyz, dbzdz=-2x^3-6xy^2+12xz^2
;27. a112=-6,a310=4,a130=-2,b400=1,b202=-3,b022=3,b220=-3,c211=-6,c031=2
;     bx=-6xyz^2+4x^3y-2xy^3, by=x^4-3x^2z^2+3y^2z^2-3x^2y^2, bz=-6x^2yz+2y^3z
;     dbxdx=-6yz^2+12x^2y-2y^3, dbxdy=-6xz^2+4x^3-6xy^2, dbxdz=-12xyz
;     dbydx=4x^3-6xz^2-6xy^2, dbydy=6yz^2-6x^2y, dbydz=-6x^2z+6y^2z
;     dbzdx=-12xyz, dbzdy=-6x^2z+6y^2z, dbzdz=-6x^2y+2y^3
;28. a112=6,a130=-2,b040=1,b202=3,b022=-3,b220=-3,c211=6,c031=-2
;     bx=6xyz^2-2xy^3, by=y^4+3x^2z^2-3y^2z^2-3x^2y^2, bz=6x^2yz-2y^3z
;     dbxdx=6yz^2-2y^3, dbxdy=6xz^2-6xy^2, dxdz=12xyz
;     dbydx=6xz^2-6xy^2, dbydy=4y^3-6yz^2-6x^2y, dbydz=6x^2z-6y^2z
;     dbzdx=12xyz, dbzdy=6x^2z-6y^2z, dbzdz=6x^2y-2y^3
;29. a112=-6,a130=2,b004=1,b202=-3,b022=-3,b220=3,c013=4,c211=-6,c031=-2
;     bx=-6xyz^2+2xy^3, by=z^4-3x^2z^2-3y^2z^2+3x^2y^2, bz=4yz^3-6x^2yz-2y^3z
;     dbxdx=-6yz^2+2y^3, dbxdy=-6xz^2+6xy^2, dbxdz=-12xyz
;     dbydx=-6xz^2+6xy^2, dbydy=-6yz^2+6xy^2, dbydz=4z^3-6x^2z-6y^2z
;     dbzdx=-12xyz, dbzdy=4z^3-6x^2z-6y^2z, dbzdz=12yz^2-6x^2y-2y^3
;30. a103=-2,a121=-6,a301=4,b013=2,b211=-6,c400=1,c220=-3,c202=-3,c022=3
;     bx=-2xz^3-6xy^2z+4x^3z, by=2yz^3-6x^2yz, bz=x^4-3x^2y^2-3x^2z^2+3y^2z^2
;     dbxdx=-2z^3-6y^2z+12x^2z, dbxdy=-12xyz, dbxdz=-6xz^2-6xy^2+4x^3
;     dbydx=-12xyz, dbydy=2z^3-6x^2y, dbydz=6yz^2-6x^2y
;     dbzdx=4x^3-6xy^2-6xz^2, dbzdy=-6x^2y+6yz^2, dbzdz=-6x^2z+6y^2z
;31. a103=2,a121=-6,b013=-2,b031=4,b211=-6,c040=1,c220=-3,c202=3,c022=-3
;     bx=2xz^3-6xy^2z, by=-2yz^3+4y^3z-6x^2yz, bz=y^4-3x^2y^2+3x^2z^2-3y^2z^2
;     dbxdx=2z^3-6y^2z, dbxdy=-12xyz, dbxdz=6xz^2-6xy^2
;     dbydx=-12xyz, dbydy=-2z^3+12y^2z-6x^2z, dbydz=-6yz^2+4y^3-6x^2y
;     dbzdx=-6xy^2+6xz^2, dbzdy=4y^3-6x^2y-6yz^2, dbzdz=6x^2z-6y^2z
;32. a103=-2,a121=6,b013=-2,b211=6,c004=1,c220=3,c202=-3,c022=-3
;     bx=-2xz^3+6xy^2z, by=-2yz^3+6x^2yz, bz=z^4+3x^2y^2-3x^2z^2-3y^2z^2
;     dbxdx=-2z^3+6y^2z, dbxdy=12xyz, dbxdz=-6xz^2+6xy^2
;     dbydx=12xyz, dbydy=-2z^3+6x^2y, dbydz=-6yz^2+6x^2y
;     dbzdx=6xy^2-6xz^2, dbzdy=6x^2y-6yz^2, dbzdz=4z^3-6x^2z-6y^2z
;33. a211=3,a013=-1,b103=-1,b301=1,c112=-3,c310=1
;     bx=3x^2yz-yz^3, by=-xz^3+x^3z, bz=-3xyz^2+x^3y
;     dbxdx=6xyz, dbxdy=3x^2z-z^3, dbxdz=3x^2y-3yz^2
;     dbydx=-z^3+3x^2z, dbydz=-3xz^2+x^3
;     dbzdx=-3yz^2+3x^2y, dbzdy=-3xz^2+x^3, dbzdz=-6xyz
;34. a031=1,a013=-1,b121=3,b103=-1,c112=-3 c130=1
;     bx=y^3z-yz^3, by=3xy^2z-xz^3, bz=-3xyz^2+xy^3
;     dbxdy=3y^2z-z^3, dbxdz=y^3-3yz^2
;     dbydx=3y^2z-z^3, dbydy=6xyz, dbydz=3xy^2-3xz^2
;     dbzdx=-3yz^2+y^3, dbzdy=-3xz^2+3xy^2, dbzdz=-6xyz
    bx+=(x^4-3*x^2*z^2-3*x^2*y^2+3*y^2*z^2)*f(*,24) & by+=(-2*x^3*y+6*x*y*z^2)*f(*,24) & bz+=(-2*x^3*z+6*x*y^2*z)*f(*,24)
    bx+=(y^4+3*x^2*z^2-3*x^2*y^2-3*y^2*z^2)*f(*,25) & by+=(-2*x^3*y+4*x*y^3-6*x*y*z^2)*f(*,25) & bz+=(2*x^3*z-6*x*y^2*z)*f(*,25)
    bx+=(z^4-3*x^2*z^2+3*x^2*y^2-3*y^2*z^2)*f(*,26) & by+=(2*x^3*y-6*x*y*z^2)*f(*,26) & bz+=(-2*x^3*z-6*x*y^2*z+4*x*z^3)*f(*,26)
    bx+=(-6*x*y*z^2+4*x^3*y-2*x*y^3)*f(*,27) & by+=(x^4-3*x^2*z^2+3*y^2*z^2-3*x^2*y^2)*f(*,27) & bz+=(-6*x^2*y*z+2*y^3*z)*f(*,27)
    bx+=(6*x*y*z^2-2*x*y^3)*f(*,28) & by+=(y^4+3*x^2*z^2-3*y^2*z^2-3*x^2*y^2)*f(*,28) & bz+=(6*x^2*y*z-2*y^3*z)*f(*,28)
    bx+=(-6*x*y*z^2+2*x*y^3)*f(*,29) & by+=(z^4-3*x^2*z^2-3*y^2*z^2+3*x^2*y^2)*f(*,29) & bz+=(4*y*z^3-6*x^2*y*z-2*y^3*z)*f(*,29)
    bx+=(-2*x*z^3-6*x*y^2*z+4*x^3*z)*f(*,30) & by+=(2*y*z^3-6*x^2*y*z)*f(*,30) & bz+=(x^4-3*x^2*y^2-3*x^2*z^2+3*y^2*z^2)*f(*,30)
    bx+=(2*x*z^3-6*x*y^2*z)*f(*,31) & by+=(-2*y*z^3+4*y^3*z-6*x^2*y*z)*f(*,31) & bz+=(y^4-3*x^2*y^2+3*x^2*z^2-3*y^2*z^2)*f(*,31)
    bx+=(-2*x*z^3+6*x*y^2*z)*f(*,32) & by+=(-2*y*z^3+6*x^2*y*z)*f(*,32) & bz+=(z^4+3*x^2*y^2-3*x^2*z^2-3*y^2*z^2)*f(*,32)
    bx+=(3*x^2*y*z-y*z^3)*f(*,33) & by+=(-x*z^3+x^3*z)*f(*,33) & bz+=(-3*x*y*z^2+x^3*y)*f(*,33)
    bx+=(y^3*z-y*z^3)*f(*,34) & by+=(3*x*y^2*z-x*z^3)*f(*,34) & bz+=(-3*x*y*z^2+x*y^3)*f(*,34)
    if keyword_set(derivs) then begin
        dbxdx+=(4*x^3-6*x*y^2-6*x*z^2)*f(*,24) & dbxdy+=(-6*x^2*y+6*y*z^2)*f(*,24) & dbxdz+=(-6*x^2*z+6*y^2*z)*f(*,24)
        dbydx+=(-6*x^2*y+6*y*z^2)*f(*,24) & dbydy+=(-2*x^3+6*x*z^2)*f(*,24) & dbydz+=(12*x*y*z)*f(*,24)
        dbzdx+=(-6*x^2*z+6*y^2*z)*f(*,24) & dbzdy+=(12*x*y*z)*f(*,24) & dbzdz+=(-2*x^3+6*x*y^2)*f(*,24)
        dbxdx+=(6*x*z^2-6*x*y^2)*f(*,25) & dbxdy+=(4*y^3-6*x^2*y-6*y*z^2)*f(*,25) & dbxdz+=(6*x^2*z-6*y^2*z)*f(*,25)
        dbydx+=(-6*x^2*y+4*y^3-6*y*z^2)*f(*,25) & dbydy+=(-2*x^3+12*x*y^2-6*x*z^2)*f(*,25) & dbydz+=(-12*x*y*z)*f(*,25)
        dbzdx+=(6*x^2*z-6*y^2*z)*f(*,25) & dbzdy+=(-12*x*y*z)*f(*,25) & dbzdz+=(2*x^3-6*x*y^2)*f(*,25)
        dbxdx+=(-6*x*z^2+6*x*y^2)*f(*,26) & dbxdy+=(6*x^2*y-6*y*z^2)*f(*,26) & dbxdz+=(4*z^3-6*x^2*z-6*y^2*z)*f(*,26)
        dbydx+=(6*x^2*y-6*y*z^2)*f(*,26) & dbydy+=(2*x^3-6*x*z^2)*f(*,26) & dbydz+=(-12*x*y*z)*f(*,26)
        dbzdx+=(-6*x^2*z-6*y^2*z+4*z^3)*f(*,26) & dbzdy+=(-12*x*y*z)*f(*,26) & dbzdz+=(-2*x^3-6*x*y^2+12*x*z^2)*f(*,26)
        dbxdx+=(-6*y*z^2+12*x^2*y-2*y^3)*f(*,27) & dbxdy+=(-6*x*z^2+4*x^3-6*x*y^2)*f(*,27) & dbxdz+=(-12*x*y*z)*f(*,27)
        dbydx+=(4*x^3-6*x*z^2-6*x*y^2)*f(*,27) & dbydy+=(6*y*z^2-6*x^2*y)*f(*,27) & dbydz+=(-6*x^2*z+6*y^2*z)*f(*,27)
        dbzdx+=(-12*x*y*z)*f(*,27) & dbzdy+=(-6*x^2*z+6*y^2*z)*f(*,27) & dbzdz+=(-6*x^2*y+2*y^3)*f(*,27)
        dbxdx+=(6*y*z^2-2*y^3)*f(*,28) & dbxdy+=(6*x*z^2-6*x*y^2)*f(*,28) & dbxdz+=(12*x*y*z)*f(*,28)
        dbydx+=(6*x*z^2-6*x*y^2)*f(*,28) & dbydy+=(4*y^3-6*y*z^2-6*x^2*y)*f(*,28) & dbydz+=(6*x^2*z-6*y^2*z)*f(*,28)
        dbzdx+=(12*x*y*z)*f(*,28) & dbzdy+=(6*x^2*z-6*y^2*z)*f(*,28) & dbzdz+=(6*x^2*y-2*y^3)*f(*,28)
        dbxdx+=(-6*y*z^2+2*y^3)*f(*,29) & dbxdy+=(-6*x*z^2+6*x*y^2)*f(*,29) & dbxdz+=(-12*x*y*z)*f(*,29)
        dbydx+=(-6*x*z^2+6*x*y^2)*f(*,29) & dbydy+=(-6*y*z^2+6*x^2*y)*f(*,29) & dbydz+=(4*z^3-6*x^2*z-6*y^2*z)*f(*,29)
        dbzdx+=(-12*x*y*z)*f(*,29) & dbzdy+=(4*z^3-6*x^2*z-6*y^2*z)*f(*,29) & dbzdz+=(12*y*z^2-6*x^2*y-2*y^3)*f(*,29)
        dbxdx+=(-2*z^3-6*y^2*z+12*x^2*z)*f(*,30) & dbxdy+=(-12*x*y*z)*f(*,30) & dbxdz+=(-6*x*z^2-6*x*y^2+4*x^3)*f(*,30)
        dbydx+=(-12*x*y*z)*f(*,30) & dbydy+=(2*z^3-6*x^2*z)*f(*,30) & dbydz+=(6*y*z^2-6*x^2*y)*f(*,30)
        dbzdx+=(4*x^3-6*x*y^2-6*x*z^2)*f(*,30) & dbzdy+=(-6*x^2*y+6*y*z^2)*f(*,30) & dbzdz+=(-6*x^2*z+6*y^2*z)*f(*,30)
        dbxdx+=(2*z^3-6*y^2*z)*f(*,31) & dbxdy+=(-12*x*y*z)*f(*,31) & dbxdz+=(6*x*z^2-6*x*y^2)*f(*,31)
        dbydx+=(-12*x*y*z)*f(*,31) & dbydy+=(-2*z^3+12*y^2*z-6*x^2*z)*f(*,31) & dbydz+=(-6*y*z^2+4*y^3-6*x^2*y)*f(*,31)
        dbzdx+=(-6*x*y^2+6*x*z^2)*f(*,31) & dbzdy+=(4*y^3-6*x^2*y-6*y*z^2)*f(*,31) & dbzdz+=(6*x^2*z-6*y^2*z)*f(*,31)
        dbxdx+=(-2*z^3+6*y^2*z)*f(*,32) & dbxdy+=(12*x*y*z)*f(*,32) & dbxdz+=(-6*x*z^2+6*x*y^2)*f(*,32)
        dbydx+=(12*x*y*z)*f(*,32) & dbydy+=(-2*z^3+6*x^2*z)*f(*,32) & dbydz+=(-6*y*z^2+6*x^2*y)*f(*,32)
        dbzdx+=(6*x*y^2-6*x*z^2)*f(*,32) & dbzdy+=(6*x^2*y-6*y*z^2)*f(*,32) & dbzdz+=(4*z^3-6*x^2*z-6*y^2*z)*f(*,32)
        dbxdx+=(6*x*y*z)*f(*,33) & dbxdy+=(3*x^2*z-z^3)*f(*,33) & dbxdz+=(3*x^2*y-3*y*z^2)*f(*,33)
        dbydx+=(-z^3+3*x^2*z)*f(*,33) & dbydz+=(-3*x*z^2+x^3)*f(*,33)
        dbzdx+=(-3*y*z^2+3*x^2*y)*f(*,33) & dbzdy+=(-3*x*z^2+x^3)*f(*,33) & dbzdz+=(-6*x*y*z)*f(*,33)
        dbxdy+=(3*y^2*z-z^3)*f(*,34) & dbxdz+=(y^3-3*y*z^2)*f(*,34)
        dbydx+=(3*y^2*z-z^3)*f(*,34) & dbydy+=(6*x*y*z)*f(*,34) & dbydz+=(3*x*y^2-3*x*z^2)*f(*,34)
        dbzdx+=(-3*y*z^2+y^3)*f(*,34) & dbzdy+=(-3*x*z^2+3*x*y^2)*f(*,34) & dbzdz+=(-6*x*y*z)*f(*,34)
    endif
endif

if keyword_set(derivs) then return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz,dbydx:dbydx,dbydy:dbydy,dbydz:dbydz,$
                                    dbzdx:dbzdx,dbzdy:dbzdy,dbzdz:dbzdz} else return,{bx:bx,by:by,bz:bz}
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_ab_fem_scalar,bels,bin,grid_spacing=grid_spacing

if not keyword_set(grid_spacing) then grid_spacing=1.0d

;bels(i).j(k)=jth component of kth element at ith point = (Y_j^k (p_i)
;from Paul's paper)
;bin(i).j=jth component of data at ith point (=Y_l(p_i))

npts = n_elements(bels)
nelem = n_elements(bels.bx)
a=dblarr(nelem,nelem)
b=dblarr(nelem)
for i=0,npts-1 do begin
    a += bels(i).bx#bels(i).bx + bels(i).by#bels(i).by + bels(i).bz#bels(i).bz +$
              grid_spacing*(bels(i).dbxdx#bels(i).dbxdx + bels(i).dbxdy#bels(i).dbxdy +bels(i).dbxdz#bels(i).dbxdz $
                            +bels(i).dbydy#bels(i).dbydy + bels(i).dbydz#bels(i).dbydz)
    b += bin(i).bx(0)*bels(i).bx + bin(i).by(0)*bels(i).by + bin(i).bz(0)*bels(i).bz +$
      grid_spacing*(bin(i).dbxdx(0)*bels(i).dbxdx + bin(i).dbxdy(0)*bels(i).dbxdy +bin(i).dbxdz(0)*bels(i).dbxdz $
                    +bin(i).dbydy(0)*bels(i).dbydy + bin(i).dbydz(0)*bels(i).dbydz)
endfor

return,{a:a,b:b}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_feminterp_sum,bels,f,derivs=derivs

bx = total(f*bels.bx)
by = total(f*bels.by)
bz = total(f*bels.bz)

if keyword_set(derivs) then begin
    dbxdx = total(f*bels.dbxdx)
    dbxdy = total(f*bels.dbxdy)
    dbxdz = total(f*bels.dbxdz)

    dbydx = total(f*bels.dbydx)
    dbydy = total(f*bels.dbydy)
    dbydz = total(f*bels.dbydz)

    dbzdx = total(f*bels.dbzdx)
    dbzdy = total(f*bels.dbzdy)
    dbzdz = total(f*bels.dbzdz)
    return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz,dbydx:dbydx,dbydy:dbydy,dbydz:dbydz,$
                                    dbzdx:dbzdx,dbzdy:dbzdy,dbzdz:dbzdz}
endif else return,{bx:bx,by:by,bz:bz}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_feminterp_sum_vec,bels,f,derivs=derivs

bx = total(f*bels.bx,2)
by = total(f*bels.by,2)
bz = total(f*bels.bz,2)

if keyword_set(derivs) then begin
    dbxdx = total(f*bels.dbxdx,2)
    dbxdy = total(f*bels.dbxdy,2)
    dbxdz = total(f*bels.dbxdz,2)

    dbydx = total(f*bels.dbydx,2)
    dbydy = total(f*bels.dbydy,2)
    dbydz = total(f*bels.dbydz,2)

    dbzdx = total(f*bels.dbzdx,2)
    dbzdy = total(f*bels.dbzdy,2)
    dbzdz = total(f*bels.dbzdz,2)
    return,{bx:bx,by:by,bz:bz,dbxdx:dbxdx,dbxdy:dbxdy,dbxdz:dbxdz,dbydx:dbydx,dbydy:dbydy,dbydz:dbydz,$
                                    dbzdx:dbzdx,dbzdy:dbzdy,dbzdz:dbzdz}
endif else return,{bx:bx,by:by,bz:bz}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_feminterp_tables_fields,g,rmp,nr=nr,nz=nz,nt=nt

if not keyword_set(nr) then nr=10
if not keyword_set(nz) then nz=10
if not keyword_set(nt) then nt=10

nr=long(nr)
nz=long(nz)
nt=long(nt)

r=(max(g.r)-min(g.r))*dindgen(nr)/(nr-1) + min(g.r)
z=(max(g.z)-min(g.z))*dindgen(nz)/(nz-1) + min(g.z)
phi = 2*!dpi*dindgen(nt)/(nt-1)

bels = get_b_elements([0,0,0])
nelem = n_elements(bels.bx)
f = dblarr(nr*nz*nt,nelem)
bels = replicate(bels,nr*nt*nz)

b = bfield_bs(1.0,0.0,0.0,rmp.coil,rmp.current,/deriv)
b = replicate(b,nr*nt*nz)

for ir=0L,nr-1 do begin
    for iz=0L,nz-1 do begin
        for it=0L,nt-1 do begin
            index = ir*nz*nt + iz*nt + it
            x = r(ir)*cos(phi(it))
            y = r(ir)*sin(phi(it))
            b(index) = bfield_bs(x,y,z(iz),rmp.coil,rmp.current,/derivs)
        endfor
    endfor
endfor

return,{r:r,z:z,phi:phi,nr:nr,nz:nz,nt:nt,b:b}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_feminterp_tables,g,rmp,nr=nr,nz=nz,nt=nt,degree=degree,quiet=quiet,dx=dx
tstart = systime(1)

if not keyword_set(nr) then nr=10
if not keyword_set(nz) then nz=10
if not keyword_set(nt) then nt=10

if keyword_set(dx) then begin
    nr = ceil((max(g.r)-min(g.r))/dx)
    nz = ceil((max(g.z)-min(g.z))/dx)
    nt = ceil(2*!dpi*g.rmaxis/dx)
    print,'NR, NZ, NT= ',nr,nz,nt
endif


nr=long(nr)
nz=long(nz)
nt=long(nt)

r=(max(g.r)-min(g.r))*dindgen(nr)/(nr-1) + min(g.r)
z=(max(g.z)-min(g.z))*dindgen(nz)/(nz-1) + min(g.z)
phi = 2*!dpi*dindgen(nt)/(nt-1)
;print,'Phi end points: ',phi(0),phi(nt-1)

bels = get_b_elements([0,0,0],/derivs,degree=degree)
nelem = n_elements(bels.bx)
f = dblarr(nr*nz*nt,nelem)
bels = replicate(bels,8)

x0 = dblarr(nr*nz*nt)
y0=x0 & z0=x0

bin = bfield_bs(1.0,0.0,0.0,rmp.coil,rmp.current,/derivs)
bin = replicate(bin,8)

for ir=0L,nr-2 do begin
    for iz=0L,nz-2 do begin
        for it=0L,nt-2 do begin
            index = ir*nz*nt + iz*nt + it
            r0 = 0.5*(r(ir)+r(ir+1))
            z0(index) = 0.5*(z(iz)+z(iz+1))
            phi0 = 0.5*(phi(it)+phi(it+1))
            x0(index) = r0*cos(phi0)
            y0(index) = r0*sin(phi0)
            bmod=0.0d
            dermod=0.0d
            for i=0,7 do begin
                ird = i mod 2
                izd = floor(i/4)
                itd = (i - 4*izd - ird)/2
                irc=ir+ird
                izc=iz+izd
                itc=it+itd
                x = r(irc)*cos(phi(itc))
                y = r(irc)*sin(phi(itc))
                bin(i) = bfield_bs(x,y,z(izc),rmp.coil,rmp.current,/derivs)
                bmod += (bin(i).bx^2 + bin(i).by^2 + bin(i).bz^2)
                dermod += (bin(i).dbxdx^2 + bin(i).dbxdy^2 + bin(i).dbxdz^2 + $
                               bin(i).dbydy^2 + bin(i). dbydz^2)
                bels(i) = get_b_elements([x-x0(index),y-y0(index),z(izc)-z0(index)],/derivs,degree=degree)
            endfor
;            print,ir,iz,it,bmod/dermod
            a=get_ab_fem_scalar(bels,bin)
            f(index,*)=invert(a.a)##a.b
        endfor
    endfor
endfor
if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'

return,{r:r,z:z,phi:phi,nr:nr,nz:nz,nt:nt,f:f,x0:x0,y0:y0,z0:z0}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_feminterp_fields,x,y,z,t,degree=degree

r = sqrt(x^2+y^2)
phi = (atan(y,x) + 2*!dpi) mod (2*!dpi)

ir = floor((r-t.r(0))/(t.r(1)-t.r(0)))
iz = floor((z-t.z(0))/(t.z(1)-t.z(0)))
it = floor((phi-t.phi(0))/(t.phi(1)-t.phi(0)))

index = ir*t.nz*t.nt + iz*t.nt + it
print,ir,iz,it,index

r0 = 0.5*(t.r(ir)+t.r(ir+1))
z0 = 0.5*(t.z(iz)+t.z(iz+1))
phi0 = 0.5*(t.phi(it)+t.phi(it+1))

x0 = r0*cos(phi0)
y0 = r0*sin(phi0)
print,x0,y0,z0
bels = get_b_elements([t.r(ir)*cos(t.phi(it)),t.r(ir)*sin(t.phi(it)),t.z(iz)],/derivs,degree=degree)
bels = replicate(bels,8)
bin =  t.b(index)
bin = replicate(bin,8)

for i=0,7 do begin
    ird = i mod 2
    izd = floor(i/4)
    itd = (i - 4*izd - ird)/2
    irc=ir+ird
    izc=iz+izd
    itc=it+itd
    indexc = irc*t.nz*t.nt + izc*t.nt + itc
    bin(i) = t.b(indexc)
    bels(i) = get_b_elements([t.r(irc)*cos(t.phi(itc))-x0,t.r(irc)*sin(t.phi(itc))-y0,t.z(izc)-z0],/derivs,degree=degree)
endfor
print,bin.bx

a=get_ab_fem_scalar(bels,bin)
f=invert(a.a)##a.b
;print,f
bels_t = get_b_elements([x-x0,y-y0,z-z0],/derivs)
bout = bfield_feminterp_sum(bels_t,f,/derivs)

return,bout
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test_bs_ring_derivs,rmp,t

nphi=100
ds = 0.00001

phi = 2*!dpi*dindgen(nphi)/(nphi-1)
r = 1.35 + 0*phi
z = -0.2 + 0*phi

bf = bfield_bs(r(0),0,z(0),rmp.coil,rmp.current,/derivs)
bf = replicate(bf,nphi)
;bfxm = replicate(bf,nphi)
;bfyp = replicate(bf,nphi)
;bfym = replicate(bf,nphi)
;bfzp = replicate(bf,nphi)
;bfzm = replicate(bf,nphi)

dbxdx=dblarr(nphi)
dbxdy=dbxdx & dbxdz=dbxdx & dbydx=dbxdx & dbydy=dbxdx & dbydz=dbxdx & dbzdx=dbxdx & dbzdy=dbxdx & dbzdz=dbxdx

for i=0,nphi-1 do begin
    x = r(i)*cos(phi(i))
    y = r(i)*sin(phi(i))
    bf(i) = bfield_bs(x,y,z(i),rmp.coil,rmp.current,/derivs)
    bfxp = bfield_bs(x+ds,y,z(i),rmp.coil,rmp.current,/derivs)
    bfxm = bfield_bs(x-ds,y,z(i),rmp.coil,rmp.current,/derivs)
    bfyp = bfield_bs(x,y+ds,z(i),rmp.coil,rmp.current,/derivs)
    bfym = bfield_bs(x,y-ds,z(i),rmp.coil,rmp.current,/derivs)
    bfzp = bfield_bs(x,y,z(i)+ds,rmp.coil,rmp.current,/derivs)
    bfzm = bfield_bs(x,y,z(i)-ds,rmp.coil,rmp.current,/derivs)
    dbxdx(i) = (bfxp.bx-bfxm.bx)/(2*ds)
    dbydx(i) = (bfxp.by-bfxm.by)/(2*ds)
    dbzdx(i) = (bfxp.bz-bfxm.bz)/(2*ds)
    dbxdy(i) = (bfyp.bx-bfym.bx)/(2*ds)
    dbydy(i) = (bfyp.by-bfym.by)/(2*ds)
    dbzdy(i) = (bfyp.bz-bfym.bz)/(2*ds)
    dbxdz(i) = (bfzp.bx-bfzm.bx)/(2*ds)
    dbydz(i) = (bfzp.by-bfzm.by)/(2*ds)
    dbzdz(i) = (bfzp.bz-bfzm.bz)/(2*ds)
end

col=mycolors()
pos = plot_positions(nrow=3,ncol=3)
!p.position=pos(*,0)
plot,phi,bf.dbxdx,xtickname=replicate(' ',6),ytitle='dbxdx'
oplot,phi,dbxdx,color=col.brick

!p.position=pos(*,1)
plot,phi,bf.dbxdy,xtickname=replicate(' ',6),ytitle='dbxdy',/noerase
oplot,phi,dbxdy,color=col.brick

!p.position=pos(*,2)
plot,phi,bf.dbxdz,xtitle='phi',ytitle='dbxdz',/noerase
oplot,phi,dbxdz,color=col.brick

!p.position=pos(*,3)
plot,phi,bf.dbydx,xtickname=replicate(' ',6),ytitle='dbydx',/noerase
oplot,phi,dbydx,color=col.brick

!p.position=pos(*,4)
plot,phi,bf.dbydy,xtickname=replicate(' ',6),ytitle='dbydy',/noerase
oplot,phi,dbydy,color=col.brick

!p.position=pos(*,5)
plot,phi,bf.dbydz,xtitle='phi',ytitle='dbydz',/noerase
oplot,phi,dbydz,color=col.brick

!p.position=pos(*,6)
plot,phi,bf.dbzdx,xtickname=replicate(' ',6),ytitle='dbzdx',/noerase
oplot,phi,dbzdx,color=col.brick

!p.position=pos(*,7)
plot,phi,bf.dbzdy,xtickname=replicate(' ',6),ytitle='dbzdy',/noerase
oplot,phi,dbzdy,color=col.brick

!p.position=pos(*,8)
plot,phi,bf.dbzdz,xtickname=replicate(' ',6),ytitle='dbzdz',/noerase
oplot,phi,dbzdz,color=col.brick
;plot,phi,bf.dbxdx+bf.dbydy+bf.dbzdz ,xtickname=replicate(' ',6),ytitle='Div, Curl',/noerase,yrange=[-0.001,0.001]
;oplot,phi,bf.dbydz-bf.dbzdy,color=col.brick
;oplot,phi,bf.dbxdz-bf.dbzdx,color=col.blue
;oplot,phi,bf.dbxdy-bf.dbydx,color=col.green

!p.position=0

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function read_ipec_rzphile,filename,keep_sign=keep_sign

restore,filename

if not keyword_set(keep_sign) then begin
    ebrs = -ebrs
    ebzs = -ebzs
    vbrs = -vbrs
    vbzs = -vbzs
    brs  = -brs
    bzs  = -bzs
endif
nt=nt+1
phi = [phi,phi(0)+2*!dpi]
ebrx=dblarr(nr,nz,nt)
ebzx = ebrx & ebpx=ebrx
vbrx = ebrx & vbzx = ebrx & vbpx=ebrx
brx = ebrx & bzx = ebrx & bpx=ebrx

ebrx(*,*,0:nt-2) = ebrs & ebrx(*,*,nt-1) = ebrs(*,*,0)
ebzx(*,*,0:nt-2) = ebzs & ebzx(*,*,nt-1) = ebzs(*,*,0)
ebpx(*,*,0:nt-2) = ebps & ebpx(*,*,nt-1) = ebps(*,*,0)
vbrx(*,*,0:nt-2) = vbrs & vbrx(*,*,nt-1) = vbrs(*,*,0)
vbzx(*,*,0:nt-2) = vbzs & vbzx(*,*,nt-1) = vbzs(*,*,0)
vbpx(*,*,0:nt-2) = vbps & vbpx(*,*,nt-1) = vbps(*,*,0)
brx(*,*,0:nt-2) = brs & brx(*,*,nt-1) = brs(*,*,0)
bzx(*,*,0:nt-2) = bzs & bzx(*,*,nt-1) = bzs(*,*,0)
bpx(*,*,0:nt-2) = bps & bpx(*,*,nt-1) = bps(*,*,0)


return,{nr:nr,nz:nz,nphi:nt,r:r,z:z,phi:!dpi*phi/180.,ebr:ebrx,ebz:ebzx,ebphi:ebpx, $
        vbr:vbrx,vbz:vbzx,vbphi:vbpx,br:brx,bz:bzx,bphi:bpx}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_ipec,ipec,r,phi,z,btype=btype

if not keyword_set(btype) then btype='ipec'

phi = phi mod (2*!dpi)
ineg = where(phi lt 0.0,count)
if count gt 0 then phi(ineg)=2*!dpi + phi(ineg)

ir = (r - ipec.r(0,0))/(ipec.r(1,0)-ipec.r(0,0))
iz = (z - ipec.z(0,0))/(ipec.z(0,1)-ipec.z(0,0))
it = (phi - ipec.phi(0))/(ipec.phi(1)-ipec.phi(0))
;print,it
case btype of
    'eq': begin
        br = interpolate(ipec.ebr,ir,iz,it,missing=0.0d)
        bz = interpolate(ipec.ebz,ir,iz,it,missing=0.0d)
        bphi = interpolate(ipec.ebphi,ir,iz,it,missing=0.0d)
    end
    'vacuum': begin
        br = interpolate(ipec.vbr,ir,iz,it,missing=0.0d)
        bz = interpolate(ipec.vbz,ir,iz,it,missing=0.0d)
        bphi = interpolate(ipec.vbphi,ir,iz,it,missing=0.0d)
    end
    'vacuumtot': begin
        br = interpolate(ipec.vbr+ipec.ebr,ir,iz,it,missing=0.0d)
        bz = interpolate(ipec.vbz+ipec.ebz,ir,iz,it,missing=0.0d)
        bphi = interpolate(ipec.vbphi+ipec.ebphi,ir,iz,it,missing=0.0d)
    end
    'ipec': begin
        br = interpolate(ipec.br,ir,iz,it,missing=0.0d)
        bz = interpolate(ipec.bz,ir,iz,it,missing=0.0d)
        bphi = interpolate(ipec.bphi,ir,iz,it,missing=0.0d)
    end
    'ipectot': begin
        br = interpolate(ipec.br+ipec.ebr,ir,iz,it,missing=0.0d)
        bz = interpolate(ipec.bz+ipec.ebz,ir,iz,it,missing=0.0d)
        bphi = interpolate(ipec.bphi+ipec.ebphi,ir,iz,it,missing=0.0d)
    end
    else: print,'Bad IPEC B-field type choice.  Avail type: eq,vacuum,ipec'
endcase

return,{br:br,bphi:bphi,bz:bz}
end

pro test_bfield_ipec,n,a

for i = 0,n-1 do b=bfield_ipec(a,1.0,0.0,0.0)

end

function read_siesta_bfield,file

close,/all
openr,1,file

iphi0=0L & is=0L & itht0=0L

readf,1,iphi0,is,itht0
iphi=iphi0+1 & itht=itht0+1 & is0=is-1

s=dblarr(is) & u=dblarr(itht) & v=dblarr(iphi)
r=dblarr(is,itht,iphi)
z=r & phi=r & br=r & bz=r & bphi=r

v=2*!dpi*dindgen(iphi)/iphi0
s=dindgen(is)/is0
u=2*!dpi*dindgen(itht)/itht0

tmp = dblarr(6)
;junk = dblarr(2)
;readf,1,junk
;print,junk

for k=0,iphi0-1 do begin
    for i=0,is-1 do begin
        for j=0,itht0-1 do begin
            readf,1,tmp
            r(i,j,k) = tmp(0)
            z(i,j,k) = tmp(1)
            phi(i,j,k) = tmp(2)
            br(i,j,k) = tmp(3)
            bz(i,j,k) = tmp(4)
            bphi(i,j,k) = tmp(5)
        endfor
    endfor
endfor

r(*,*,iphi0) = r(*,*,0)
z(*,*,iphi0) = z(*,*,0)
phi(*,*,iphi0) = phi(*,*,0)+ 2*!dpi
br(*,*,iphi0) = br(*,*,0)
bz(*,*,iphi0) = bz(*,*,0)
bphi(*,*,iphi0) = bphi(*,*,0)

r(*,itht0,*) = r(*,0,*)
z(*,itht0,*) = z(*,0,*)
phi(*,itht0,*) = phi(*,0,*)+ 2*!dpi
br(*,itht0,*) = br(*,0,*)
bz(*,itht0,*) = bz(*,0,*)
bphi(*,itht0,*) = bphi(*,0,*)

;r = reform(r,is*iphi*itht)
;phi = reform(phi,is*iphi*itht)
;z = reform(z,is*iphi*itht)
;br = reform(br,is*iphi*itht)

close,1

;stop

return,{r:r,z:z,phi:phi,br:br,bz:bz,bphi:bphi,s:s,u:u,v:v}
end


function in_tri,a,b,c,p

v0 = c-a
v1 = b-a
v2 = p-a

;dot00 = TRANSPOSE(v0) # v0
;dot01 = TRANSPOSE(v0) # v1
;dot02 = TRANSPOSE(v0) # v2
;dot11 = TRANSPOSE(v1) # v1
;dot12 = TRANSPOSE(v1) # v2

dot00 = v0(*,0)*v0(*,0) + v0(*,1)*v0(*,1)
dot01 = v0(*,0)*v1(*,0) + v0(*,1)*v1(*,1)
dot02 = v0(*,0)*v2(*,0) + v0(*,1)*v2(*,1)
dot11 = v1(*,0)*v1(*,0) + v1(*,1)*v1(*,1)
dot12 = v1(*,0)*v2(*,0) + v1(*,1)*v2(*,1)

invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
u = (dot11 * dot02 - dot01 * dot12) * invDenom
v = (dot00 * dot12 - dot01 * dot02) * invDenom
;print,u
ins = (u gt -1d-15 and v ge -1d-15 and u+v le 1.0+1d-12)
;stop
return,{ins:ins,u:u,v:v}
end

function get_bfield_singleslice,s,r,z,iphi


npts = n_elements(r)
ns = n_elements(s.r(*,0,0))
np = n_elements(s.r(0,*,0))

br = dblarr(npts)
bphi = dblarr(npts) + 1.0
bz = dblarr(npts)

for ir=0,npts-1 do begin
    r00 = s.r(0:ns-2,0:np-2,iphi(ir))
    r01 = s.r(0:ns-2,1:np-1,iphi(ir))
    r10 = s.r(1:ns-1,0:np-2,iphi(ir))
    r11 = s.r(1:ns-1,1:np-1,iphi(ir))

    z00 = s.z(0:ns-2,0:np-2,iphi(ir))
    z01 = s.z(0:ns-2,1:np-1,iphi(ir))
    z10 = s.z(1:ns-1,0:np-2,iphi(ir))
    z11 = s.z(1:ns-1,1:np-1,iphi(ir))

;     br00 = s.br(0:ns-2,0:ns-2,iphi(ir))
;     br01 = s.br(0:ns-2,1:ns-1,iphi(ir))
;     br10 = s.br(1:ns-1,0:ns-2,iphi(ir))
;     br11 = s.br(1:ns-1,1:ns-1,iphi(ir))

;     bz00 = s.bz(0:ns-2,0:ns-2,iphi(ir))
;     bz01 = s.bz(0:ns-2,1:ns-1,iphi(ir))
;     bz10 = s.bz(1:ns-1,0:ns-2,iphi(ir))
;     bz11 = s.bz(1:ns-1,1:ns-1,iphi(ir))

;     bphi00 = s.bphi(0:ns-2,0:ns-2,iphi(ir))
;     bphi01 = s.bphi(0:ns-2,1:ns-1,iphi(ir))
;     bphi10 = s.bphi(1:ns-1,0:ns-2,iphi(ir))
;     bphi11 = s.bphi(1:ns-1,1:ns-1,iphi(ir))

;    iencr = where(((r00 le r(ir)) or r01 le r(ir) or r10 le r(ir) or r11 le r(ir)) and (r00 gt r(ir)) or r01 gt r(ir) or r10 gt r(ir) or r11 gt r(ir))
;    iencz = where(((z00 le z(ir)) or z01 le z(ir) or z10 le z(ir) or z11 le z(ir)) and (z00 gt z(ir)) or z01 gt z(ir) or z10 gt z(ir) or z11 gt z(ir))
;    ienctot = where(iencr eq iencz)
    ienc = where(((r00 le r(ir) or r01 le r(ir) or r10 le r(ir) or r11 le r(ir)) and (r00 gt r(ir) or r01 gt r(ir) or r10 gt r(ir) or r11 gt r(ir))) and $
                 (z00 le z(ir) or z01 le z(ir) or z10 le z(ir) or z11 le z(ir)) and (z00 gt z(ir) or z01 gt z(ir) or z10 gt z(ir) or z11 gt z(ir)),count)
;stop
    if count gt 0 then begin
        inds =  array_indices(r00,ienc)
    endif else begin
;        print,'Point is off grid!'
;        return,0
    endelse

    for i=0,count-1 do begin
;        print,'Starting'
        bar = in_tri([[r00(inds(0,i),inds(1,i))],[z00(inds(0,i),inds(1,i))]],$
                     [[r11(inds(0,i),inds(1,i))],[z11(inds(0,i),inds(1,i))]],$
                     [[r10(inds(0,i),inds(1,i))],[z10(inds(0,i),inds(1,i))]],[[r(ir)],[z(ir)]])
;        print,'U,V;INS: ',bar.u,bar.v,bar.ins
        if bar.ins eq 1 then begin
            br(ir) = s.br(inds(0,i),inds(1,i),iphi(ir)) + $
              bar.u*(s.br(inds(0,i)+1,inds(1,i),iphi(ir)) - s.br(inds(0,i),inds(1,i),iphi(ir))) + $
              bar.v*(s.br(inds(0,i)+1,inds(1,i)+1,iphi(ir)) - s.br(inds(0,i),inds(1,i),iphi(ir)))
            bphi(ir) = s.bphi(inds(0,i),inds(1,i),iphi(ir)) + $
              bar.u*(s.bphi(inds(0,i)+1,inds(1,i),iphi(ir)) - s.bphi(inds(0,i),inds(1,i),iphi(ir))) + $
              bar.v*(s.bphi(inds(0,i)+1,inds(1,i)+1,iphi(ir)) - s.bphi(inds(0,i),inds(1,i),iphi(ir)))
            bz(ir) = s.bz(inds(0,i),inds(1,i),iphi(ir)) + $
              bar.u*(s.bz(inds(0,i)+1,inds(1,i),iphi(ir)) - s.bz(inds(0,i),inds(1,i),iphi(ir))) + $
              bar.v*(s.bz(inds(0,i)+1,inds(1,i)+1,iphi(ir)) - s.bz(inds(0,i),inds(1,i),iphi(ir)))
;             print,[[r00(inds(0,i),inds(1,i))],[z00(inds(0,i),inds(1,i))]],$
;                      [[r11(inds(0,i),inds(1,i))],[z11(inds(0,i),inds(1,i))]],$
;                      [[r10(inds(0,i),inds(1,i))],[z10(inds(0,i),inds(1,i))]]
;            br(ir) = br00(inds(0,i),inds(1,i)) + bar.u*(br10(inds(0,i),inds(1,i)) - br00(inds(0,i),inds(1,i))) + bar.v*(br11(inds(0,i),inds(1,i)) - br00(inds(0,i),inds(1,i)))
;            bphi(ir) = bphi00(inds(0,i),inds(1,i)) + bar.u*(bphi10(inds(0,i),inds(1,i)) - bphi00(inds(0,i),inds(1,i))) + bar.v*(bphi11(inds(0,i),inds(1,i)) - bphi00(inds(0,i),inds(1,i)))
;            bz(ir) = bz00(inds(0,i),inds(1,i)) + bar.u*(bz10(inds(0,i),inds(1,i)) - bz00(inds(0,i),inds(1,i))) + bar.v*(bz11(inds(0,i),inds(1,i)) - bz00(inds(0,i),inds(1,i)))
            break
        endif
        bar = in_tri([[r00(inds(0,i),inds(1,i))],[z00(inds(0,i),inds(1,i))]],$
                     [[r01(inds(0,i),inds(1,i))],[z01(inds(0,i),inds(1,i))]],$
                     [[r11(inds(0,i),inds(1,i))],[z11(inds(0,i),inds(1,i))]],[[r(ir)],[z(ir)]])
;        print,'U,V;INS: ',bar.u,bar.v,bar.ins
        if bar.ins eq 1 then begin
            br(ir) = s.br(inds(0,i),inds(1,i),iphi(ir)) + $
              bar.u*(s.br(inds(0,i)+1,inds(1,i)+1,iphi(ir)) - s.br(inds(0,i),inds(1,i),iphi(ir))) + $
              bar.v*(s.br(inds(0,i),inds(1,i)+1,iphi(ir)) - s.br(inds(0,i),inds(1,i),iphi(ir)))
            bphi(ir) = s.bphi(inds(0,i),inds(1,i),iphi(ir)) + $
              bar.u*(s.bphi(inds(0,i)+1,inds(1,i)+1,iphi(ir)) - s.bphi(inds(0,i),inds(1,i),iphi(ir))) + $
              bar.v*(s.bphi(inds(0,i),inds(1,i)+1,iphi(ir)) - s.bphi(inds(0,i),inds(1,i),iphi(ir)))
            bz(ir) = s.bz(inds(0,i),inds(1,i),iphi(ir)) + $
              bar.u*(s.bz(inds(0,i)+1,inds(1,i)+1,iphi(ir)) - s.bz(inds(0,i),inds(1,i),iphi(ir))) + $
              bar.v*(s.bz(inds(0,i),inds(1,i)+1,iphi(ir)) - s.bz(inds(0,i),inds(1,i),iphi(ir)))
;             print,[[r00(inds(0,i),inds(1,i))],[z00(inds(0,i),inds(1,i))]],$
;                      [[r01(inds(0,i),inds(1,i))],[z01(inds(0,i),inds(1,i))]],$
;                      [[r11(inds(0,i),inds(1,i))],[z11(inds(0,i),inds(1,i))]]
;            br(ir) = br00(inds(0,i),inds(1,i)) + bar.u*(br11(inds(0,i),inds(1,i)) - br00(inds(0,i),inds(1,i))) + bar.v*(br01(inds(0,i),inds(1,i)) - br00(inds(0,i),inds(1,i)))
;            bphi(ir) = bphi00(inds(0,i),inds(1,i)) + bar.u*(bphi11(inds(0,i),inds(1,i)) - bphi00(inds(0,i),inds(1,i))) + bar.v*(bphi01(inds(0,i),inds(1,i)) - bphi00(inds(0,i),inds(1,i)))
;            bz(ir) = bz00(inds(0,i),inds(1,i)) + bar.u*(bz11(inds(0,i),inds(1,i)) - bz00(inds(0,i),inds(1,i))) + bar.v*(bz01(inds(0,i),inds(1,i)) - bz00(inds(0,i),inds(1,i)))
            break
        endif
    endfor
    if bar.ins ne 1 then print,'I don''t think you should get here: ',r(ir),z(ir)
endfor

return,{br:br,bphi:bphi,bz:bz}
end


function bfield_siesta,s,r,phi,z

phil = phi mod (2*!dpi)

dphi = s.phi(0,0,1)-s.phi(0,0,0)
iphi = floor(phil/dphi)
tphi = phil/dphi-iphi

pt = iphi+1
iw = where(pt gt n_elements(s.r(0,0,*))-1,count)
if count gt 0 then pt(iw)=1

b1=get_bfield_singleslice(s,r,z,iphi)
b2=get_bfield_singleslice(s,r,z,pt)

br = (1-tphi)*b1.br + tphi*b2.br
bz = (1-tphi)*b1.bz + tphi*b2.bz
bphi = (1-tphi)*b1.bphi + tphi*b2.bphi

;print,phi,iphi,tphi,pt


return,{br:br,bphi:bphi,bz:bz}
end

function siestacoord_singleslice,s,r,z,iphi


npts = n_elements(r)
ns = n_elements(s.r(*,0,0))
np = n_elements(s.r(0,*,0))

s_sta = dblarr(npts)
u_sta = dblarr(npts)
v_sta = dblarr(npts)

for ir=0L,npts-1 do begin
    r00 = s.r(0:ns-2,0:np-2,iphi(ir))
    r01 = s.r(0:ns-2,1:np-1,iphi(ir))
    r10 = s.r(1:ns-1,0:np-2,iphi(ir))
    r11 = s.r(1:ns-1,1:np-1,iphi(ir))

    z00 = s.z(0:ns-2,0:np-2,iphi(ir))
    z01 = s.z(0:ns-2,1:np-1,iphi(ir))
    z10 = s.z(1:ns-1,0:np-2,iphi(ir))
    z11 = s.z(1:ns-1,1:np-1,iphi(ir))


    ienc = where(((r00 le r(ir) or r01 le r(ir) or r10 le r(ir) or r11 le r(ir)) and (r00 gt r(ir) or r01 gt r(ir) or r10 gt r(ir) or r11 gt r(ir))) and $
                 (z00 le z(ir) or z01 le z(ir) or z10 le z(ir) or z11 le z(ir)) and (z00 gt z(ir) or z01 gt z(ir) or z10 gt z(ir) or z11 gt z(ir)),count)
;stop
    if count gt 0 then begin
        inds =  array_indices(r00,ienc)
    endif else begin
;        print,'Point is off grid!'
;        return,0
    endelse
;    print,'Starting search'
    for i=0,count-1 do begin
        bar = in_tri([[r00(inds(0,i),inds(1,i))],[z00(inds(0,i),inds(1,i))]],$
                     [[r11(inds(0,i),inds(1,i))],[z11(inds(0,i),inds(1,i))]],$
                     [[r10(inds(0,i),inds(1,i))],[z10(inds(0,i),inds(1,i))]],[[r(ir)],[z(ir)]])
;        print,[bar.u,bar.v,bar.ins]
        if bar.ins eq 1 then begin
            s_sta(ir) = s.s(inds(0,i)) + $
              bar.u*(s.s(inds(0,i)+1) - s.s(inds(0,i))) + $
              bar.v*(s.s(inds(0,i)+1) - s.s(inds(0,i)))
            u_sta(ir) = s.u(inds(1,i)) + $
              bar.u*(s.u(inds(1,i)+1) - s.u(inds(1,i))) + $
              bar.v*(s.u(inds(1,i)+1) - s.u(inds(1,i)))
            v_sta(ir) = s.v(iphi(ir))
            break
        endif
        bar = in_tri([[r00(inds(0,i),inds(1,i))],[z00(inds(0,i),inds(1,i))]],$
                     [[r01(inds(0,i),inds(1,i))],[z01(inds(0,i),inds(1,i))]],$
                     [[r11(inds(0,i),inds(1,i))],[z11(inds(0,i),inds(1,i))]],[[r(ir)],[z(ir)]])
;        print,[bar.u,bar.v,bar.ins]
        if bar.ins eq 1 then begin
            s_sta(ir) = s.s(inds(0,i)) + $
              bar.u*(s.s(inds(0,i)+1) - s.s(inds(0,i))) + $
              bar.v*(s.s(inds(0,i)) - s.s(inds(0,i)))
            u_sta(ir) = s.u(inds(1,i)) + $
              bar.u*(s.u(inds(1,i)+1) - s.u(inds(1,i))) + $
              bar.v*(s.u(inds(1,i)) - s.u(inds(1,i)))
            v_sta(ir) = s.v(iphi(ir))
            break
        endif
    endfor
;    if bar.ins ne 1 then print,'I don''t think you should get here'
endfor
;stop
return,{s:s_sta,u:u_sta,v:v_sta}
end

function fluxcoord_siesta,s,r,phi,z

phil = phi mod (2*!dpi)

dphi = s.phi(0,0,1)-s.phi(0,0,0)
iphi = floor(phil/dphi)
tphi = phil/dphi-iphi

pt = iphi+1
iw = where(pt gt n_elements(s.r(0,0,*))-1,count)
if count gt 0 then pt(iw)=1

b1=siestacoord_singleslice(s,r,z,iphi)
b2=siestacoord_singleslice(s,r,z,pt)
;print,b1
;print,b2
s_sta = (1-tphi)*b1.s + tphi*b2.s
u_sta = (1-tphi)*b1.u + tphi*b2.u
v_sta = (1-tphi)*b1.v + tphi*b2.v

;print,phi,iphi,tphi,pt


return,{s:s_sta,u:u_sta,v:v_sta}
end


function convert_siesta_coords,s,a

npts = n_elements(s.r(0,*))
;npts = n_elements(s.r(*,0))
for i=0L,npts-1 do begin
    sta = fluxcoord_siesta(a,s.r(*,i),s.phi(*,i),s.z(*,i))
    s.psiN(*,i) = sta.s
;    sta = fluxcoord_siesta(a,s.r(i,*),s.phi(i,*),s.z(i,*))
;    s.psiN(i,*) = sta.s
endfor

return,s
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; function field_line_derivs_phi,phi,x
; COMMON pass_rmp, rmpon, rmpinfo

; R = x(0)
; z = x(1)

; if rmpon then B = bfield_components(R,phi,z,addrmp=rmpinfo) $
;   else B = bfield_components(R,phi,z)
; dRdphi = R*B.r/B.phi
; dzdphi = R*B.z/B.phi
; dxdphi = [dRdphi,dzdphi]
; return, dxdphi
; end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function follow_fieldlines_phi,g,rmp,pnstart,pnend,nsurfs,ntheta,phistart,dphi,nsteps,zero=zero, $
                               rstart=rstart,period=period,zstart=zstart,ipec=ipec,btype=btype,quiet=quiet, $
                               diffuse=diffuse,allsamestart=allsamestart,siesta=siesta,pnarr=pnarr,fem=fem, $
                               secfem=secfem,secrmp=secrmp,sigtheta=sigtheta,sigphi=sigphi,m3dc1=m3dc1,efit=efit, $
                               extender=extender,wout=wout,biext=biext
tstart = systime(1)
r = dblarr(nsurfs*ntheta,nsteps+1)
z=r & phi=r

if nsurfs gt 1 then $
  pnvals = (pnend-pnstart)*dindgen(nsurfs)/(nsurfs-1) + pnstart $
else pnvals = pnstart
phi(*,0) = phistart
;z(*,0) = 0.0

if keyword_set(pnarr) then begin
    pnvals=pnarr
    nsurfs = n_elements(pnarr)
    r = dblarr(nsurfs*ntheta,nsteps+1)
    z=r & phi=r
endif

;pnvals = (pnend-pnstart)*dindgen(nsurfs)/(nsurfs-1) + pnstart
;print,pnvals
if ntheta gt 1 then theta = 2*!dpi*dindgen(ntheta)/ntheta else theta=0.0
nline = 1e4
for i=0,ntheta-1 do begin
    Lmax = max([g.r(g.mw-1)-g.r(0),g.z(g.mh-1) - g.z(0)])
    rline = g.rmaxis+Lmax*cos(theta(i))*dindgen(nline)/(nline-1)
    zline = g.zmaxis+Lmax*sin(theta(i))*dindgen(nline)/(nline-1)
    psiline = (-g.cpasma/abs(g.cpasma)*get_psi_bicub(g,rline,zline)-g.ssimag)/(g.ssibry-g.ssimag)
    ipn1 = where(psiline(0:nline-2) le 1.001d and psiline(1:nline-1) gt 1.001d,count)
    if count ge 1 then ipn1 = ipn1(0)
    iline = interpol(dindgen(ipn1+1),psiline(0:ipn1),pnvals)
    r(ntheta*dindgen(nsurfs)+i,0)=interpolate(rline,iline)
;    print,'Indices = ',ntheta*dindgen(nsurfs)+i
    z(ntheta*dindgen(nsurfs)+i,0)=interpolate(zline,iline)
;    print,'HERE it is! ',interpolate(psiline,iline)
;    for j=0,nsurfs-1 do begin
endfor

if keyword_set(siesta) then begin
    for i=0,nsurfs-1 do begin
        junk = min(siesta.s-pnvals(i),/abs,ind_s)
        ts = (pnvals(i) - siesta.s(ind_s))/(siesta.s(ind_s+1)-siesta.s(ind_s))
        r1 = interpol(siesta.r(ind_s,*,0),siesta.u,theta)
        z1 = interpol(siesta.z(ind_s,*,0),siesta.u,theta)
        r2 = interpol(siesta.r(ind_s+1,*,0),siesta.u,theta)
        z2 = interpol(siesta.z(ind_s+1,*,0),siesta.u,theta)
;        print,siesta.s(ind_s),ts
        r(ntheta*i:(ntheta*(i+1)-1),0) = r1*(1.0d - ts) + r2*ts
        z(ntheta*i:(ntheta*(i+1)-1),0) = z1*(1.0d - ts) + z2*ts
;        r(ntheta*i:(ntheta*(i+1)-1),0) = interpol(siesta.r(ind_s,*,0),siesta.u,theta)
;        z(ntheta*i:(ntheta*(i+1)-1),0) = interpol(siesta.z(ind_s,*,0),siesta.u,theta)
    endfor
;    map = fluxcoord_siesta(siesta,r(*,0),phi(*,0),z(*,0))
;stop
endif

if keyword_set(wout) then begin
   phiwant = 0.0d
    for i=0,nsurfs-1 do begin
       r1 = dblarr(ntheta)
       r2=r1 & z1=r1 & z2=r1
       junk = min(wout.s-pnvals(i),/abs,ind_s)
       if ind_s eq wout.ns-1 then ind_s = ind_s-1
       ts = (pnvals(i) - wout.s(ind_s))/(wout.s(ind_s+1)-wout.s(ind_s))
       for k=0,ntheta-1 do begin
          r1(k) = total(wout.rmnc(*,ind_s)*cos(wout.xm*theta(k)-wout.xn*phiwant))
          if wout.lasym then r1(k) = r1(k) + total(wout.rmns(*,ind_s)*sin(wout.xm*theta(k)-wout.xn*phiwant))
          r2(k) = total(wout.rmnc(*,ind_s+1)*cos(wout.xm*theta(k)-wout.xn*phiwant))
          if wout.lasym then r2(k) = r2(k) + total(wout.rmns(*,ind_s+1)*sin(wout.xm*theta(k)-wout.xn*phiwant))
          z1(k) = total(wout.zmns(*,ind_s)*sin(wout.xm*theta(k)-wout.xn*phiwant))
          if wout.lasym then z1(k) = z1(k) + total(wout.zmnc(*,ind_s)*cos(wout.xm*theta(k)-wout.xn*phiwant))
          z2(k) = total(wout.zmns(*,ind_s+1)*sin(wout.xm*theta(k)-wout.xn*phiwant))
          if wout.lasym then z2(k) = z2(k) + total(wout.zmnc(*,ind_s+1)*cos(wout.xm*theta(k)-wout.xn*phiwant))
          r(ntheta*i:(ntheta*(i+1)-1),0) = r1*(1.0d - ts) + r2*ts
          z(ntheta*i:(ntheta*(i+1)-1),0) = z1*(1.0d - ts) + z2*ts
;        r(ntheta*i:(ntheta*(i+1)-1),0) = interpol(siesta.r(ind_s,*,0),siesta.u,theta)
;        z(ntheta*i:(ntheta*(i+1)-1),0) = interpol(siesta.z(ind_s,*,0),siesta.u,theta)
       endfor
    endfor
;    map = fluxcoord_siesta(siesta,r(*,0),phi(*,0),z(*,0))
;stop
endif

if keyword_set(allsamestart) then begin
    r(*,0) = r(0,0)
    z(*,0) = z(0,0)
endif

if keyword_set(rstart) and not keyword_set(zstart) then begin
    r(*,0) = pnvals
    z(*,0) = 0.
endif

if keyword_set(zstart) and not keyword_set(rstart) then begin
    r(*,0) = zstart
    z(*,0) = pnvals
endif

if keyword_set(zstart) and keyword_set(rstart) then begin
    r = dblarr(n_elements(rstart),nsteps+1)
    z=r & phi=r
    phi(*,0) = phistart
    r(*,0) = rstart
    z(*,0) = zstart
endif

if 1 then begin
for i=1L,nsteps do begin
    phi(*,i)=phi(*,i-1)+dphi

    b=bfield_geq_bicub(g,r(*,i-1),z(*,i-1),zero=zero)
;    print,'old ',b.br
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1),phi(*,i-1),z(*,i-1),rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secrmp) then begin
        b2=bfield_bs_cyl(r(*,i-1),phi(*,i-1),z(*,i-1),secrmp.coil,secrmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1),phi(*,i-1),z(*,i-1),btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
     endif
    if keyword_set(m3dc1) then begin
        b2=bfield_m3dc1(r(*,i-1),phi(*,i-1),z(*,i-1),m3dc1,efit=efit)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
     endif
    if keyword_set(extender) then begin
        b2=bfield_extender(r(*,i-1),phi(*,i-1),z(*,i-1),extender,biext=biext)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
    endif
    if keyword_set(siesta) then begin
        b2=bfield_siesta(siesta,r(*,i-1),phi(*,i-1),z(*,i-1))
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
    endif
    if keyword_set(fem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1),phi(*,i-1),z(*,i-1),fem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secfem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1),phi(*,i-1),z(*,i-1),secfem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    k1r = dphi*r(*,i-1)*b.br/b.bphi
    k1z = dphi*r(*,i-1)*b.bz/b.bphi

    b=bfield_geq_bicub(g,r(*,i-1)+0.5*k1r,z(*,i-1)+0.5*k1z,zero=zero)
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secrmp) then begin
        b2=bfield_bs_cyl(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,secrmp.coil,secrmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
     endif
    if keyword_set(m3dc1) then begin
        b2=bfield_m3dc1(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,m3dc1,efit=efit)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
     endif
    if keyword_set(extender) then begin
        b2=bfield_extender(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,extender,biext=biext)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
     endif
    if keyword_set(siesta) then begin
        b2=bfield_siesta(siesta,r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(fem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,fem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secfem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,secfem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    k2r = dphi*(r(*,i-1)+0.5*k1r)*b.br/b.bphi
    k2z = dphi*(r(*,i-1)+0.5*k1r)*b.bz/b.bphi

    b=bfield_geq_bicub(g,r(*,i-1)+0.5*k2r,z(*,i-1)+0.5*k2z,zero=zero)
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secrmp) then begin
        b2=bfield_bs_cyl(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,secrmp.coil,secrmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
     endif
    if keyword_set(m3dc1) then begin
        b2=bfield_m3dc1(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,m3dc1,efit=efit)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
     endif
    if keyword_set(extender) then begin
        b2=bfield_extender(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,extender,biext=biext)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
    endif
    if keyword_set(siesta) then begin
        b2=bfield_siesta(siesta,r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(fem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,fem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secfem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,secfem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    k3r = dphi*(r(*,i-1)+0.5*k2r)*b.br/b.bphi
    k3z = dphi*(r(*,i-1)+0.5*k2r)*b.bz/b.bphi

    b=bfield_geq_bicub(g,r(*,i-1)+k3r,z(*,i-1)+k3z,zero=zero)
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secrmp) then begin
        b2=bfield_bs_cyl(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,secrmp.coil,secrmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
     endif
    if keyword_set(m3dc1) then begin
        b2=bfield_m3dc1(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,m3dc1,efit=efit)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
     endif
    if keyword_set(extender) then begin
        b2=bfield_extender(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,extender,biext=biext)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
    endif
    if keyword_set(siesta) then begin
        b2=bfield_siesta(siesta,r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(fem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,fem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(secfem) then begin
        b2=bfield_feminterp_cyl(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,secfem)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    k4r = dphi*(r(*,i-1)+k3r)*b.br/b.bphi
    k4z = dphi*(r(*,i-1)+k3r)*b.bz/b.bphi


    r(*,i) = r(*,i-1) + k1r/6 + k2r/3 + k3r/3 + k4r/6
    z(*,i) = z(*,i-1) + k1z/6 + k2z/3 + k3z/3 + k4z/6

    if keyword_set(diffuse) then begin
;        dL = r(*,i)*(phi(*,i)-phi(*,i-1))
        dL = r(*,i)*abs(phi(*,i)-phi(*,i-1))
;        x1 = r(i)*cos(phi(i))
;        x2 = r(i-1)*cos(phi(i-1))
;        y1 = r(i)*sin(phi(i))
;        y2 = r(i-1)*sin(phi(i-1))
;        z1 = z(i)
;        z2 = z(i-1)
;        dL2 = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
;        print,dL,dL2
    b=bfield_geq_bicub(g,r(*,i),z(*,i),zero=zero)
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i),phi(*,i),z(*,i),rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i),phi(*,i),z(*,i),btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(siesta) then begin
        b2=bfield_siesta(siesta,r(*,i),phi(*,i),z(*,i))
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
        perpdir1 = crossp([b.br,b.bphi,b.bz],[0,0,1])
        perpdir1 = perpdir1/sqrt(perpdir1(0)^2 + perpdir1(1)^2 + perpdir1(2)^2)
        perpdir2 = crossp([b.br,b.bphi,b.bz],[1,0,0])
        perpdir2 = perpdir2/sqrt(perpdir2(0)^2 + perpdir2(1)^2 + perpdir2(2)^2)
        alpha = 2*!pi*(-1 + 2*randomu(seed,nsurfs*ntheta))
        diffdir = cos(alpha)*perpdir1 + sin(alpha)*perpdir2
        delta_x = sqrt(diffuse*dL)
        if keyword_set(sigtheta) then begin
           difftheta = atan(z(*,i)-g.zmaxis,r(*,i)-g.rmaxis)*180.0d/!dpi
           amptheta = (1/(sqrt(2*!dpi)*sigtheta))*exp(-difftheta^2/(2*sigtheta^2))*250.0d
           delta_x = sqrt(diffuse*amptheta*dL)
;           stop
        endif
        if keyword_set(sigphi) then begin
           sigphi_n = floor(sigphi)
           sigphi_eps = sigphi - sigphi_n
           delta_x = delta_x*(1+sigphi_eps*cos(sigphi_n*phi(*,i)))
;           stop
        endif
;        print,'delta_x',delta_x
;        print,'bphi',b.bphi
;        print,'perpdir1',perpdir1
;        print,'perpdir2',perpdir2
;        print,'Change in r:',delta_x*(cos(alpha)*perpdir1(0) + sin(alpha)*perpdir2(0))
;        print, 'Change in phi: ',delta_x*(cos(alpha)*perpdir1(1) + sin(alpha)*perpdir2(1))
;        print,'Change in z: ',delta_x*(cos(alpha)*perpdir1(2) + sin(alpha)*perpdir2(2))
        r(*,i) = r(*,i) + delta_x*(cos(alpha)*perpdir1(0) + sin(alpha)*perpdir2(0))
        phi(*,i) = phi(*,i) + delta_x*(cos(alpha)*perpdir1(1) + sin(alpha)*perpdir2(1))
        z(*,i) = z(*,i) + delta_x*(cos(alpha)*perpdir1(2) + sin(alpha)*perpdir2(2))
    endif

endfor

psiN = 0.*r
for i=0,nsurfs*ntheta-1 do begin
    psiN(i,*) = -g.cpasma/abs(g.cpasma)*get_psi_bicub(g,r(i,*),z(i,*))
endfor
;psiN = -g.cpasma/abs(g.cpasma)*get_psi_bicub(g,r,z)
psiN = (psiN-g.ssimag)/(g.ssibry-g.ssimag)

;psiN = reform(psiN,nsurfs*ntheta,nsteps+1)

theta = atan(z-g.zmaxis,r-g.rmaxis)
ineg = where(theta lt 0.0,count)
if count gt 0 then theta(ineg) = theta(ineg) + 2*!pi
if not keyword_set(period) then period=2*!dpi
nip0 = round(nsteps*abs(dphi)/period)
if nip0 gt 0 then ip0 = round(period/abs(dphi))*lindgen(nip0) else ip0=0
;print,ip0
;ip0 = floor(dindgen(floor((nsteps+1)*abs(dphi)/period))*period/dphi)
;print,ip0
endif
if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'
return,{r:r,z:z,phi:phi,psiN:psiN,theta:theta,ip0:ip0}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function follow_fieldlines_rzphi,g,rmp,rstart,zstart,phistart,dphi,nsteps,zero=zero, $
                                 period=period,ipec=ipec,btype=btype,quiet=quiet, $
                                 diffuse=diffuse
tstart = systime(1)

npts = n_elements(rstart)

r = dblarr(npts,nsteps+1)
z=r & phi=r

r(*,0) = rstart
z(*,0) = zstart
phi(*,0) = phistart

for i=1L,nsteps do begin
    phi(*,i)=phi(*,i-1)+dphi

    b=bfield_geq_bicub(g,r(*,i-1),z(*,i-1),zero=zero)
;    print,'old ',b.br
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1),phi(*,i-1),z(*,i-1),rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1),phi(*,i-1),z(*,i-1),btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
;        print,'new ',b2.br
    endif
    k1r = dphi*r(*,i-1)*b.br/b.bphi
    k1z = dphi*r(*,i-1)*b.bz/b.bphi

    b=bfield_geq_bicub(g,r(*,i-1)+0.5*k1r,z(*,i-1)+0.5*k1z,zero=zero)
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    k2r = dphi*(r(*,i-1)+0.5*k1r)*b.br/b.bphi
    k2z = dphi*(r(*,i-1)+0.5*k1r)*b.bz/b.bphi

    b=bfield_geq_bicub(g,r(*,i-1)+0.5*k2r,z(*,i-1)+0.5*k2z,zero=zero)
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    k3r = dphi*(r(*,i-1)+0.5*k2r)*b.br/b.bphi
    k3z = dphi*(r(*,i-1)+0.5*k2r)*b.bz/b.bphi

    b=bfield_geq_bicub(g,r(*,i-1)+k3r,z(*,i-1)+k3z,zero=zero)
    if abs(max(rmp.current,/abs)) gt 1d-8 then begin
        b2=bfield_bs_cyl(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,rmp.coil,rmp.current)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    if keyword_set(ipec) then begin
        b2=bfield_ipec(ipec,r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,btype=btype)
        b.br += b2.br
        b.bphi += b2.bphi
        b.bz += b2.bz
    endif
    k4r = dphi*(r(*,i-1)+k3r)*b.br/b.bphi
    k4z = dphi*(r(*,i-1)+k3r)*b.bz/b.bphi


    r(*,i) = r(*,i-1) + k1r/6 + k2r/3 + k3r/3 + k4r/6
    z(*,i) = z(*,i-1) + k1z/6 + k2z/3 + k3z/3 + k4z/6

    if keyword_set(diffuse) then begin
        dL = r(*,i)*(phi(*,i)-phi(*,i-1))
;        x1 = r(i)*cos(phi(i))
;        x2 = r(i-1)*cos(phi(i-1))
;        y1 = r(i)*sin(phi(i))
;        y2 = r(i-1)*sin(phi(i-1))
;        z1 = z(i)
;        z2 = z(i-1)
;        dL2 = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
;        print,dL,dL2
        b=bfield_geq_bicub(g,r(*,i),z(*,i),zero=zero)
        if abs(max(rmp.current,/abs)) gt 1d-8 then begin
            b2=bfield_bs_cyl(r(*,i),phi(*,i),z(*,i),rmp.coil,rmp.current)
            b.br += b2.br
            b.bphi += b2.bphi
            b.bz += b2.bz
        endif
        if keyword_set(ipec) then begin
            b2=bfield_ipec(ipec,r(*,i),phi(*,i),z(*,i),btype=btype)
            b.br += b2.br
            b.bphi += b2.bphi
            b.bz += b2.bz
        endif
        perpdir1 = crossp([b.br,b.bphi,b.bz],[0,0,1])
        perpdir1 = perpdir1/sqrt(perpdir1(0)^2 + perpdir1(1)^2 + perpdir1(2)^2)
        perpdir2 = crossp([b.br,b.bphi,b.bz],[1,0,0])
        perpdir2 = perpdir2/sqrt(perpdir2(0)^2 + perpdir2(1)^2 + perpdir2(2)^2)
        alpha = 2*!pi*(-1 + 2*randomu(seed,nsurfs*ntheta))
        diffdir = cos(alpha)*perpdir1 + sin(alpha)*perpdir2
        delta_x = sqrt(diffuse*dL)
;        print,'Change in r:',delta_x*(cos(alpha)*perpdir1(0) + sin(alpha)*perpdir2(0))
;        print, 'Change in phi: ',delta_x*(cos(alpha)*perpdir1(1) + sin(alpha)*perpdir2(1))
;        print,'Change in z: ',delta_x*(cos(alpha)*perpdir1(2) + sin(alpha)*perpdir2(2))
        r(*,i) = r(*,i) + delta_x*(cos(alpha)*perpdir1(0) + sin(alpha)*perpdir2(0))
        phi(*,i) = phi(*,i) + delta_x*(cos(alpha)*perpdir1(1) + sin(alpha)*perpdir2(1))
        z(*,i) = z(*,i) + delta_x*(cos(alpha)*perpdir1(2) + sin(alpha)*perpdir2(2))
    endif

endfor

if not keyword_set(period) then period=2*!dpi
nip0 = round(nsteps*abs(dphi)/period)
if nip0 gt 0 then ip0 = round(period/abs(dphi))*lindgen(nip0) else ip0=0
if keyword_set(psitheta) then begin 
    psiN = 0.*r
    for i=0,npts-1 do begin
        psiN(i,*) = -g.cpasma/abs(g.cpasma)*get_psi_bicub(g,r(i,*),z(i,*))
    endfor
    psiN = (psiN-g.ssimag)/(g.ssibry-g.ssimag)

    theta = atan(z-g.zmaxis,r-g.rmaxis)
    ineg = where(theta lt 0.0,count)
    if count gt 0 then theta(ineg) = theta(ineg) + 2*!pi
    if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'
    return,{r:r,z:z,phi:phi,psiN:psiN,theta:theta,ip0:ip0}
endif else begin
    for i=0,npts-1 do begin
        ipol = where(z(i,0:nsteps-1) lt zstart(i) and z(i,1:nsteps) ge zstart(i))
        tpol = phi(i,ipol)
;        stop
    endfor
    if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'
    return,{r:r,z:z,phi:phi,ip0:ip0}
endelse
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_qpunct,g,pnstart,pnend,nsurfs,phistart,dphi,nsteps


rmp=build_nstx_rwmcoils(0d3)
s=follow_fieldlines_phi(g,rmp,pnstart,pnend,nsurfs,1,phistart,dphi,nsteps)

pn = s.psin(*,0)

dth = s.theta(0:nsurfs-1,1:nsteps)-s.theta(0:nsurfs-1,0:nsteps-1)

ineg = where(dth lt 0.0,count)
if count gt 0 then dth(ineg) += 2*!dpi

theta=0*s.phi
theta(*,1:nsteps)=total(dth,2,/cumulative)
qtmp=s.phi/theta

;ipnz = s.ip0(where(s.ip0 ne 0))
ipnz = s.ip0(where(s.ip0 ge 0.5*max(s.ip0)))

;plot,s.phi(*,ipnz),qtmp(*,ipnz)

qpunct = dblarr(nsurfs)
qerr = qpunct

for i=0,nsurfs-1 do begin
    qpunct(i) = mean(qtmp(i,ipnz))
    qerr(i) = stddev(qtmp(i,ipnz))/sqrt(n_elements(ipnz))
endfor

return,{pn:pn,q:qpunct,err:qerr}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function find_vf_fixaxis,coil,R0,nfp=nfp,nsteps=nsteps,nvf=nvf,vfmin=vfmin,vfmax=vfmax,quiet=quiet

if not keyword_set(nfp) then nfp=1
if not keyword_set(nvf) then nvf=11
if not keyword_set(vfmax) then vfmax=0.5
if not keyword_set(vfmin) then vfmin=-0.5
if not keyword_set(nsteps) then nsteps = 360/nfp

tstart = systime(1)
r = dblarr(nvf,nsteps+1)
z=r & phi=r

r(*,0) = R0
dphi = 2*!dpi/(nfp*nsteps)
period = 2*!dpi/nfp

vf = (vfmax-vfmin)*dindgen(nvf)/(nvf-1) + vfmin

if 1 then begin
for i=1L,nsteps do begin
    phi(*,i)=phi(*,i-1)+dphi
    b=bfield_bs_cyl(r(*,i-1),phi(*,i-1),z(*,i-1),coil.coil,coil.current)
    b.bz += vf
    k1r = dphi*r(*,i-1)*b.br/b.bphi
    k1z = dphi*r(*,i-1)*b.bz/b.bphi

    b=bfield_bs_cyl(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,coil.coil,coil.current)
    b.bz += vf
    k2r = dphi*(r(*,i-1)+0.5*k1r)*b.br/b.bphi
    k2z = dphi*(r(*,i-1)+0.5*k1r)*b.bz/b.bphi

    b=bfield_bs_cyl(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,coil.coil,coil.current)
    b.bz += vf
    k3r = dphi*(r(*,i-1)+0.5*k2r)*b.br/b.bphi
    k3z = dphi*(r(*,i-1)+0.5*k2r)*b.bz/b.bphi

    b=bfield_bs_cyl(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,coil.coil,coil.current)
    b.bz += vf
    k4r = dphi*(r(*,i-1)+k3r)*b.br/b.bphi
    k4z = dphi*(r(*,i-1)+k3r)*b.bz/b.bphi


    r(*,i) = r(*,i-1) + k1r/6 + k2r/3 + k3r/3 + k4r/6
    z(*,i) = z(*,i-1) + k1z/6 + k2z/3 + k3z/3 + k4z/6
endfor

endif
if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'
return,{r0:r(*,0),r1:r(*,nsteps),z0:z(*,0),z1:z(*,nsteps),vf:vf}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function punctplot_vf,coil,rstart,rend,nr,phistart,dphi,nsteps,vf=vf,nfp=nfp,quiet=quiet

if not keyword_set(nfp) then nfp=1
if not keyword_set(vf) then vf=0.0

tstart = systime(1)
r = dblarr(nr,nsteps+1)
z=r & phi=r
if nr gt 1 then r(*,0) = (rend-rstart)*dindgen(nr)/(nr-1) + rstart else r(*,0) = rstart

period = 2*!dpi/nfp
for i=1L,nsteps do begin
    phi(*,i)=phi(*,i-1)+dphi
    b=bfield_bs_cyl(r(*,i-1),phi(*,i-1),z(*,i-1),coil.coil,coil.current)
    b.bz += vf
    k1r = dphi*r(*,i-1)*b.br/b.bphi
    k1z = dphi*r(*,i-1)*b.bz/b.bphi

    b=bfield_bs_cyl(r(*,i-1)+0.5*k1r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k1z,coil.coil,coil.current)
    b.bz += vf
    k2r = dphi*(r(*,i-1)+0.5*k1r)*b.br/b.bphi
    k2z = dphi*(r(*,i-1)+0.5*k1r)*b.bz/b.bphi

    b=bfield_bs_cyl(r(*,i-1)+0.5*k2r,phi(*,i-1)+0.5*dphi,z(*,i-1)+0.5*k2z,coil.coil,coil.current)
    b.bz += vf
    k3r = dphi*(r(*,i-1)+0.5*k2r)*b.br/b.bphi
    k3z = dphi*(r(*,i-1)+0.5*k2r)*b.bz/b.bphi

    b=bfield_bs_cyl(r(*,i-1)+k3r,phi(*,i-1)+dphi,z(*,i-1)+k3z,coil.coil,coil.current)
    b.bz += vf
    k4r = dphi*(r(*,i-1)+k3r)*b.br/b.bphi
    k4z = dphi*(r(*,i-1)+k3r)*b.bz/b.bphi

    r(*,i) = r(*,i-1) + k1r/6 + k2r/3 + k3r/3 + k4r/6
    z(*,i) = z(*,i-1) + k1z/6 + k2z/3 + k3z/3 + k4z/6
endfor

nip0 = round(nsteps*abs(dphi)/period)
if nip0 gt 0 then ip0 = round(period/abs(dphi))*lindgen(nip0) else ip0=0

if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'
return,{r:r,z:z,phi:phi,ip0:ip0}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function field_line_derivs_phi,phi,x,private

r = x(0)
z = x(1)
b=bfield_geq_bicub(private.g,r,z,zero=zero)
if abs(max(private.rmp.current,/abs)) gt 1d-8 then begin
    b2=bfield_bs_cyl(r,phi,z,private.rmp.coil,private.rmp.current)
    b.br += b2.br
    b.bphi += b2.bphi
    b.bz += b2.bz
endif
drdphi = r*b.br/b.bphi
dzdphi = r*b.bz/b.bphi
dsdphi = r*sqrt(b.br^2+b.bphi^2+b.bz^2)/b.bphi
dxdphi = [drdphi,dzdphi,dsdphi]

return, dxdphi
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function follow_fieldlines_totarg,g,rmp,pnstart,pnend,nsurfs,ntheta,phistart,nsteps,dphi=dphi,zero=zero, $
                                  rstart=rstart,period=period,zstart=zstart,ipec=ipec,btype=btype,quiet=quiet, $
                                  diffuse=diffuse,allsamestart=allsamestart,tol=tol
tstart = systime(1)
r = dblarr(nsurfs*ntheta,nsteps+1)
z=r & phi=r & lconn=r

if not keyword_set(tol) then tol=1d-8
if not keyword_set(dphi) then dphi = 2*!dpi

vth = sqrt(2*100*1.6d-19/9.11d-31)

if nsurfs gt 1 then $
  pnvals = (pnend-pnstart)*dindgen(nsurfs)/(nsurfs-1) + pnstart $
else pnvals = pnstart
phi(*,0) = phistart
;z(*,0) = 0.0

;pnvals = (pnend-pnstart)*dindgen(nsurfs)/(nsurfs-1) + pnstart
;print,pnvals
if ntheta gt 1 then theta = 2*!dpi*dindgen(ntheta)/ntheta else theta=0.0
nline = 1e4
for i=0,ntheta-1 do begin
    Lmax = max([g.r(g.mw-1)-g.r(0),g.z(g.mh-1) - g.z(0)])
    rline = g.rmaxis+Lmax*cos(theta(i))*dindgen(nline)/(nline-1)
    zline = g.zmaxis+Lmax*sin(theta(i))*dindgen(nline)/(nline-1)
    psiline = (-g.cpasma/abs(g.cpasma)*get_psi_bicub(g,rline,zline)-g.ssimag)/(g.ssibry-g.ssimag)
    ipn1 = where(psiline(0:nline-2) le 1.001d and psiline(1:nline-1) gt 1.001d,count)
    if count ge 1 then ipn1 = ipn1(0)
    iline = interpol(dindgen(ipn1+1),psiline(0:ipn1),pnvals)
    r(ntheta*dindgen(nsurfs)+i,0)=interpolate(rline,iline)
;    print,'Indices = ',ntheta*dindgen(nsurfs)+i
    z(ntheta*dindgen(nsurfs)+i,0)=interpolate(zline,iline)
;    print,'HERE it is! ',interpolate(psiline,iline)
;    for j=0,nsurfs-1 do begin
endfor
if keyword_set(allsamestart) then begin
    r(*,0) = r(0,0)
    z(*,0) = z(0,0)
endif

if keyword_set(rstart) and not keyword_set(zstart) then begin
    r(*,0) = pnvals
    z(*,0) = 0.
endif

if keyword_set(zstart) and not keyword_set(rstart) then begin
    r(*,0) = zstart
    z(*,0) = pnvals
endif

if keyword_set(zstart) and keyword_set(rstart) then begin
    r = dblarr(n_elements(rstart),nsteps+1)
    z=r & phi=r & lconn=r
    phi(*,0) = phistart
    r(*,0) = rstart
    z(*,0) = zstart
endif

grmp = {g:g,rmp:rmp}
if 1 then begin
    for j=0L,nsurfs*ntheta-1 do begin         ;nsurfs*ntheta-1
        xs = [r(j,0),z(j,0),0.d0]
        phis = phi(j,i-1)
        for i=1L,nsteps do begin     ;nsteps do begin
            atol = tol
            phie = phis + dphi
;            print,phis,phie
            ddeabm,'field_line_derivs_phi',phis,xs,phie,grmp,epsrel=0.d0,epsabs=atol,status=status
            while status lt 0 and atol lt 1e-3 do begin
;                print,atol
                ddeabm,'field_line_derivs_phi',phis,xs,phi(j,i-1)+dphi,grmp,epsrel=rtol,epsabs=atol,status=status,state=state               
                atol = 5*atol
            endwhile
            phi(j,i) = phis
;            print,phis,xs
           
            lconn(j,i) = xs(2)

            if keyword_set(diffuse) then begin
                dL = abs(lconn(j,i)-lconn(j,i-1))
                randmag= randomu(seed)
                randangle = randomu(seed)
                drho = sqrt(12.0*diffuse*dL/vth)*randmag
                dr = drho*cos(2*!dpi*randangle)
                dz = drho*sin(2*!dpi*randangle)
;                print,'Dr,Dz,DL,drho= ',dr,dz,dL,drho
                xs(0) += dr
                xs(1) += dz
            endif

            r(j,i) = xs(0)
            z(j,i) = xs(1)
        endfor
    endfor

    psiN = 0.*r
    for i=0,nsurfs*ntheta-1 do begin
        psiN(i,*) = -g.cpasma/abs(g.cpasma)*get_psi_bicub(g,r(i,*),z(i,*))
    endfor
;psiN = -g.cpasma/abs(g.cpasma)*get_psi_bicub(g,r,z)
    psiN = (psiN-g.ssimag)/(g.ssibry-g.ssimag)

;psiN = reform(psiN,nsurfs*ntheta,nsteps+1)

    theta = atan(z-g.zmaxis,r-g.rmaxis)
    ineg = where(theta lt 0.0,count)
    if count gt 0 then theta(ineg) = theta(ineg) + 2*!pi
    if not keyword_set(period) then period=2*!dpi
endif
if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'
return,{r:r,z:z,phi:phi,lconn:lconn,psiN:psiN,theta:theta}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function footprints,g,rmp,rstart=rstart,rend=rend,zdiv=zdiv,phirange=phirange,nr=nr,nphi=nphi,ntransits=ntransits,ipec=ipec, $
                    btype=btype,zero=zero,fem=fem,secfem=secfem,secrmp=secrmp,limlaunch=limlaunch

tstart = systime(1)

if not keyword_set(rstart) then rstart = 0.3
if not keyword_set(rend) then rend = 0.6
if not keyword_set(zdiv) then zdiv=-1.603
if not keyword_set(nr) then nr=10
if not keyword_set(nphi) then nphi=10
if not keyword_set(phirange) then phirange=[0,2*!dpi]
if not keyword_set(ntransits) then ntransits = 10

rr = (rend-rstart)*dindgen(nr)/(nr-1)+rstart
if nphi eq 1 then phi=phirange(0) else phi = (phirange(1)-phirange(0))*dindgen(nphi)/(nphi-1) + phirange(0)
psimin = dblarr(nr,nphi)
qmin = psimin
lconn = psimin
phiconn = psimin
psistart = psimin
insidesep = psimin

il = where(g.lim(0,*) gt 1e-4)
rlim = g.lim(0,il)
zlim = g.lim(1,il)
vesobj = obj_new('IDLanROI',rlim,zlim)
sepobj = obj_new('IDLanROI',g.bdry(0,*),g.bdry(1,*))

if keyword_set(limlaunch) then begin
   ilow = where(zlim lt 0.0)
   zdiv = interpol(zlim(ilow),rlim(ilow),rr)+0.001
;   stop
endif else zdiv = zdiv + 0.0*rr

for i = 0,nr-1 do begin
    for j=0,nphi-1 do begin
        s=follow_fieldlines_phi(g,rmp,0.,0.,1,1,phi(j),!dpi/180,ntransits*2*!dpi/(!dpi/180),rstart=rr(i),zstart=zdiv(i),$
                                /quiet,ipec=ipec,btype=btype,zero=zero,fem=fem,secfem=secfem,secrmp=secrmp)
;        print,i,j,s.r(0),s.phi(0),s.z(0),phi(j)
        ins = vesobj->containspoints(s.r,s.z)
        ithit = where(ins lt 0.5,count)
        if count eq 0 then ithit = n_elements(s.r)-1 else ithit = ithit(0)
        psimin(i,j) = min(s.psiN(0:ithit))
        qmin(i,j) = interpol(g.qpsi,g.pn,psimin(i,j))
        phiconn(i,j) = s.phi(ithit)
        psistart(i,j) = s.psiN(0)
        x = s.r*cos(s.phi)
        y = s.r*sin(s.phi)
        dl = sqrt((x(1:ithit)-x(0:ithit-1))^2 + $
        (y(1:ithit)-y(0:ithit-1))^2 + $
                  (s.z(1:ithit)-s.z(0:ithit-1))^2)
        lconn(i,j) = total(dl)
        core = sepobj->containspoints(s.r(0:ithit),s.z(0:ithit))
        junk = where(core gt 0.5,count)
        if count gt 0 then insidesep(i,j) = 1
    endfor
endfor
print,'elapsed time: ',systime(1)-tstart,' sec'

return,{rr:rr,phi:phi,psistart:psistart,psimin:psimin,qmin:qmin,phiconn:phiconn,lconn:lconn,insidesep:insidesep}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro plot_footprints,f,ps=ps,yrange=yrange,ystyle=ystyle

if keyword_set(ps) then begin
   set_plot,'ps'
   !p.font=0
   device,file='footprints.ps',/inches,xsize=4,ysize=3.5,xoff=.5,yoff=.5,/helvetica,/bold,/color,bits_per_pixel=8
endif

make_nice_plots,ps=ps
z=mycolors()
pos = plot_positions(nrow=1,ncol=1,vertspace=0.1)
pos(2,*) = 0.8
;print,pos(*,1)
cpos = pos
cpos(0,*) = 0.9
cpos(2,*) = 0.99

!p.position=pos(*,0)
levs = (max(f.lconn)-min(f.lconn))*dindgen(50)/49 + min(f.lconn)
contour,transpose(f.lconn),f.phi/!dpi,f.rr,levels=levs,/fill,xtitle=textoidl('\phi/\pi'),ytitle=textoidl('R (m)'),yrange=yrange,ystyle=ystyle
colorbar,range=[min(levs),max(levs)],/vertical,format='(F5.1)',position=cpos(*,0),title=textoidl('L_{conn}'),charsize=0.7

if keyword_set(ps) then begin
   device,/close
   set_plot,'x'
endif

!p.position=0
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function process_fieldlines_phi,g,s

tmp = size(s.r,/dimension)
nlines =tmp(0)
nphi = tmp(1)
il = where(g.lim(0,*) gt 1e-4)
rlim = transpose(g.lim(0,il))
zlim = transpose(g.lim(1,il))
ves_obj = Obj_New('IDLanROI', rlim,zlim)
;inside = ves_obj->ContainsPoints(s.r,s.z)
;inside = reform(inside,nlines,nphi)

rhit=0.
phihit=0.
zhit=0.
lchit=0.
rpunct=0.
zpunct=0.
lcpunct=0.

nhit=0
;x = s.r*cos(s.phi)
;y = s.r*sin(s.phi)
for i = 0,nlines-1 do begin
    x = s.r(i,*)*cos(s.phi(i,*))
    y = s.r(i,*)*sin(s.phi(i,*))
    zz = s.z(i,*)
    dl = sqrt((x(1:nphi-1)-x(0:nphi-2))^2 + (y(1:nphi-1)-y(0:nphi-2))^2 + (zz(1:nphi-1)-zz(0:nphi-2))^2)
;    dl = sqrt((x(i,1:nphi-1)-x(i,0:nphi-2))^2 + (y(i,1:nphi-1)-y(i,0:nphi-2))^2 + (s.z(i,1:nphi-1)-s.z(i,0:nphi-2))^2)
    lconntmp = [0.,total(dl,/cumulative)]
    lconnmax = total(dl)
    inside = ves_obj->containspoints(s.r(i,*),s.z(i,*))
    if inside(0) eq 0 then break
    iout = where(inside eq 0,count)
    if count gt 0 then begin
        iout = iout(0)
        nhit+=1
        if rhit(0) eq 0. then begin
;            lconn = lconntmp(iout)
            rhit = s.r(i,iout)
            zhit = s.z(i,iout)
            lchit = lconntmp(iout)
            phihit = s.phi(i,iout)
            ipk = where(s.ip0 le iout,countipk)
            if countipk gt 0 then begin
                rpunct = transpose(s.r(i,s.ip0(ipk)))
                zpunct = transpose(s.z(i,s.ip0(ipk)))
                lcpunct = lconntmp(iout) - lconntmp(s.ip0(ipk))
            endif
        endif else begin
;            lconn = lconntmp(iout)
            rhit = [rhit,s.r(i,iout)]
            zhit = [zhit,s.z(i,iout)]
            lchit = [lchit,lconntmp(iout)]
            phihit = [phihit,s.phi(i,iout)]
            ipk = where(s.ip0 le phihit,countipk)
            if countipk gt 0 then begin
                rpunct = [rpunct,transpose(s.r(i,s.ip0(ipk)))]
                zpunct = [zpunct,transpose(s.z(i,s.ip0(ipk)))]
                lcpunct = [lcpunct,lconntmp(iout)-lconntmp(s.ip0(ipk))]
            endif
        endelse
    endif else begin
        if rpunct(0) eq 0. then begin
;            lconn = lconntmp(iout)
            rpunct = transpose(s.r(i,s.ip0))
            zpunct = transpose(s.z(i,s.ip0))
            lcpunct = lconntmp(s.ip0)*0.0 + lconnmax
;            endif
        endif else begin
;            lconn = lconntmp(iout)
            rpunct = [rpunct,transpose(s.r(i,s.ip0))]
            zpunct = [zpunct,transpose(s.z(i,s.ip0))]
            lcpunct = [lcpunct,lconntmp(s.ip0)*0.0 + lconnmax]
        endelse
    endelse
endfor

obj_destroy,ves_obj
return,{nhit:nhit,rhit:rhit,phihit:phihit,zhit:zhit,lchit:lchit,rpunct:rpunct,zpunct:zpunct,lcpunct:lcpunct}
;return, {rhit:rhit,zhit:zhit,phit:phit,nhit:nhit,psihit:psihit,lconn:lconn,psistart:psistart}
end


function lconn_contours,g,rmp,extender=extender,zero=zero,r0=r0,r1=r1,nr=nr,z0=z0,z1=z1,nz=nz,nphi=nphi,m3dc1=m3dc1,biext=biext

if not keyword_set(r0) then r0=1.0
if not keyword_set(r1) then r1=1.7
if not keyword_set(z0) then z0=-1.4
if not keyword_set(z1) then z1=-0.8
if not keyword_set(nr) then nr=4
if not keyword_set(nz) then nz=4

if not keyword_set(nphi) then nphi=3600

r = dblarr(nr,nz)
z = dblarr(nr,nz)

r1d = (r1-r0)*dindgen(nr)/(nr-1) + r0
z1d = (z1-z0)*dindgen(nz)/(nz-1) + z0

for iz=0,nz-1 do begin
   r(*,iz) = r1d
   z(*,iz) = z1d(iz)
endfor

rarr = reform(r,nr*nz)
zarr = reform(z,nr*nz)

sf=follow_fieldlines_phi(g,rmp,0.5,0.5,1,1,0.0d,!dpi/180,nphi,extender=extender,zero=zero,rstart=rarr,zstart=zarr,m3dc1=m3dc1,biext=biext)
sb=follow_fieldlines_phi(g,rmp,0.5,0.5,1,1,0.0d,-!dpi/180,nphi,extender=extender,zero=zero,rstart=rarr,zstart=zarr,m3dc1=m3dc1,biext=biext)

il = where(g.lim(0,*) gt 1e-4)
rlim = transpose(g.lim(0,il))
zlim = transpose(g.lim(1,il))
ves_obj = Obj_New('IDLanROI', rlim,zlim)

insf = ves_obj->containspoints(sf.r,sf.z)
insf = reform(insf,nr*nz,nphi+1)

xf = sf.r*cos(sf.phi)
yf = sf.r*sin(sf.phi)
zf = sf.z

dlf = sqrt((shift(xf,0,-1)-xf)^2 + (shift(yf,0,-1)-yf)^2 + (shift(zf,0,-1)-zf)^2)
lcf = total(dlf,2,/cumulative)*insf

insb = ves_obj->containspoints(sb.r,sb.z)
insb = reform(insb,nr*nz,nphi+1)

xb = sb.r*cos(sb.phi)
yb = sb.r*sin(sb.phi)
zb = sb.z

dlb = sqrt((shift(xb,0,-1)-xb)^2 + (shift(yb,0,-1)-yb)^2 + (shift(zb,0,-1)-zb)^2)
lcb = total(dlb,2,/cumulative)*insb

lconn=min([[max(lcf,dimension=2)],[max(lcb,dimension=2)]],dim=2)
;stop
lconn=reform(lconn,nr,nz)

return,{r:r1d,z:z1d,lconn:lconn}
end

function decimate_fieldlines,s,ndec=ndec

if not keyword_set(ndec) then ndec=10

nf = n_elements(s.r(0,*))
ns = n_elements(s.r(*,0))

nf2 = n_elements(s.r(0,0:nf-1:ndec))
;print,nf2

s2 = {r:dblarr(ns,nf2),z:dblarr(ns,nf2),phi:dblarr(ns,nf2),psiN:dblarr(ns,nf2),theta:dblarr(ns,nf2),L:dblarr(ns,nf2)}

s2.r = s.r(*,0:nf-1:ndec)
s2.z = s.z(*,0:nf-1:ndec)
s2.phi = s.phi(*,0:nf-1:ndec)
s2.psiN = s.psiN(*,0:nf-1:ndec)
s2.theta = s.theta(*,0:nf-1:ndec)

xx = s.r*cos(s.phi)
yy = s.r*sin(s.phi)
zz = s.z
dL = sqrt((xx(*,1:nf-1)-xx(*,0:nf-2))^2 + $
          (yy(*,1:nf-1)-yy(*,0:nf-2))^2 + $
          (zz(*,1:nf-1)-zz(*,0:nf-2))^2)
L = total(dL,2,/cumulative)
L2 = dblarr(ns,nf)
L2(*,1:nf-1) = L

s2.L = L2(*,0:nf-1:ndec)

return,s2
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mag_diff,g,s,redge=redge,inside=inside

if not keyword_set(redge) then redge=0.85

tmp=size(s.r(*,*),/dim)
nlines = tmp(0)
npts = tmp(1)

if not keyword_set(inside) then begin
    il = where(g.lim(0,*) gt 1e-4)
    rlim = transpose(g.lim(0,il))
    zlim = transpose(g.lim(1,il))
    ves_obj = Obj_New('IDLanROI', rlim,zlim)
    inside = ves_obj->containspoints(s.r(*,0:npts-2),s.z(*,0:npts-2))
    inside = reform(inside,nlines,npts-1)
endif

psiNtest = 1d-3
psi0 = s.psiN(0,0)
ipsi = where(abs(s.psiN(*,0)-psi0) lt psiNtest)
psi = mean(s.psiN(ipsi,0))
npol = n_elements(ipsi)
npsi = nlines/npol
;print,npsi,psi
for i=1,npsi-1 do begin
    ipg = where(s.psiN(*,0) gt psi(i-1)+psiNtest)
    psinext = s.psiN(ipg(0))
;    junk = min(s.psiN(*,0)-psi(i-1)-psiNtest,/abs,ipnext)
;    psinext = s.psiN(ipnext,0)
    ipsinext = where(abs(s.psiN(*,0)-psinext) lt psiNtest)
    psi = [psi,mean(s.psiN(ipsinext,0))]
    ipsi = [[ipsi],[ipsinext]]
endfor
;help,ipsi
;print,npsi,psi

rhod = redge*sqrt(interpol(g.psitor,g.pn,s.psiN))
xx = s.r*cos(s.phi)
yy = s.r*sin(s.phi)
zz = s.z
dL = sqrt((xx(*,1:npts-1)-xx(*,0:npts-2))^2 + $
          (yy(*,1:npts-1)-yy(*,0:npts-2))^2 + $
          (zz(*,1:npts-1)-zz(*,0:npts-2))^2)
;dL = sqrt((s.r(*,1:npts-1)-s.r(*,0:npts-2))^2 + $
;          (s.z(*,1:npts-1)-s.z(*,0:npts-2))^2)
L = total(dL,2,/cumulative)
drsq = (rhod(*,1:npts-1)-rhod(*,intarr(npts-1)))^2
dr = (rhod(*,1:npts-1)-rhod(*,0:npts-2))
psistart = s.psiN(*,0:npts-2)

Lave = dblarr(npsi,npts-1)
drsqave = dblarr(npsi,npts-1)
vr = dblarr(npsi)
D_std = vr
D = dblarr(npsi)
D2 = D

dpsi = mean(psi(1:npsi-1)-psi(0:npsi-2))/2
for i=0,npsi-1 do begin
    Lave(i,*) = total(L(ipsi(*,i),*),1)/npol
;    help,total(L(ipsi(*,i),*),1)
    drsqave(i,*) = total(drsq(ipsi(*,i),*),1)/npol
    manydr = reform(dr(ipsi(*,i),*),npol*(npts-1))
    manydL = reform(dL(ipsi(*,i),*),npol*(npts-1))
    ins = reform(inside(ipsi(*,i),*),npol*(npts-1))
    vr(i)= mean(manydr(ins))
    D_std(i) = stddev(manydr(ins))^2/(2*mean(manydL(ins)))
    D(i) = mean(manydr(ins)^2/(2*manydL(ins)))
;    print,'N1, i= ',i,': ',n_elements(ins)

    ik = where(psistart ge psi(i)-dpsi and psistart lt psi(i)+dpsi)
    ins = inside(ik)
;    print,'N2, i= ',i,': ',n_elements(ins)
    mdr = dr(ik)
    mdl = dL(ik)
;    plot,dL(ik)
;    D2(i) = mean(dr(ik)^2/dL(ik))/2
    D2(i) = mean(mdr(ins)^2/mdl(ins))/2
;    bins = histogram(manydr,location=xbin,nbins=50)
;    window,i
;    plot,xbin,bins
endfor

; print,'Vr =',vr
; print,'std =',D_std;std^2/(2*meandL)
; print,'D =',D

; for i = nlines/2,nlines/2 do begin            ;nlines-1 do begin
;     dL = sqrt((s.r(i,1:npts-1)-s.r(i,0:npts-2))^2 + $
;               (s.z(i,1:npts-1)-s.z(i,0:npts-2))^2)
;     L = total(dL,/cumulative)
;     drsq = (rhod(i,1:npts-1)-rhod(i,0))^2
; endfor

return,{L:L,drsq:drsq,Lave:Lave,drsqave:drsqave,psi:psi,D:D,D_std:D_std,vr:vr,D2:D2,in:inside}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mag_diff_plots,g,s,redge=redge,plotit=plotit,NL=NL,ps=ps,pdf=pdf,inside=inside,leave_s=leave_s,Lfit=Lfit

if not keyword_set(redge) then redge=1.0d
if not keyword_set(NL) then NL = 1000
if not keyword_set(plotit) then plotit=-1

if keyword_set(ps) then begin
    devold = !d.name
    set_plot,'ps'
    device,file='drhosq_vs_L.ps',/color,/inches,xsize=5,ysize=9,xoffset=.75,yoffset=.75
endif


tmp=size(s.r(*,*),/dim)
nlines = tmp(0)
npts = tmp(1)

if not keyword_set(inside) then begin
    il = where(g.lim(0,*) gt 1e-4)
    rlim = transpose(g.lim(0,il))
    zlim = transpose(g.lim(1,il))
    ves_obj = Obj_New('IDLanROI', rlim,zlim)
    inside = ves_obj->containspoints(s.r(*,0:npts-2),s.z(*,0:npts-2))
    inside = reform(inside,nlines,npts-1)
endif

psiNtest = 1d-3
psi0 = s.psiN(0,0)
ipsi = where(abs(s.psiN(*,0)-psi0) lt psiNtest)
psi = mean(s.psiN(ipsi,0))
npol = n_elements(ipsi)
npsi = nlines/npol
for i=1,npsi-1 do begin
    ipg = where(s.psiN(*,0) gt psi(i-1)+psiNtest)
    psinext = s.psiN(ipg(0))
    ipsinext = where(abs(s.psiN(*,0)-psinext) lt psiNtest)
    psi = [psi,mean(s.psiN(ipsinext,0))]
    ipsi = [[ipsi],[ipsinext]]
endfor

if keyword_set(leave_s) then $
  s_siesta=psi $
else s_siesta = sqrt(interpol(g.psitor/g.psitor(n_elements(g.psitor)-2),g.pn,psi))

if keyword_set(leave_s) then $
  rhod = redge*s.psiN $
else rhod = redge*double(sqrt(interpol(g.psitor,g.pn,s.psiN)))


if where(tag_names(s) eq strupcase('L')) gt -1 then begin
    L = s.L(*,1:npts-1)
    dL = s.L(*,1:npts-1) - s.L(*,0:npts-2)
endif else begin
    xx = s.r*cos(s.phi)
    yy = s.r*sin(s.phi)
    zz = s.z
    dL = sqrt((xx(*,1:npts-1)-xx(*,0:npts-2))^2 + $
              (yy(*,1:npts-1)-yy(*,0:npts-2))^2 + $
              (zz(*,1:npts-1)-zz(*,0:npts-2))^2)
    L = total(dL,2,/cumulative)
endelse
drsq = (rhod(*,1:npts-1)-rhod(*,intarr(npts-1)))^2
;dr = (rhod(*,1:npts-1)-rhod(*,0:npts-2))
dr = (rhod(*,1:npts-1)-rhod(*,lonarr(npts-1)))

psistart = s.psiN(*,0:npts-2)

Lave = dblarr(npsi,NL)
drsqave = Lave
insave = Lave

D = dblarr(npsi)
slope = D
tdep = D
Lmin = 1e20 + dblarr(npsi)
if npsi gt 1 then dpsi = mean(psi(1:npsi-1)-psi(0:npsi-2))/2 else dpsi=1.0
for i=0,npsi-1 do begin
    Lave(i,*) = max(L(ipsi(*,i),*))*dindgen(NL)/(NL-1)
;    plot,[0,1000],[-1,-1],yrange=[0,20]
    for j=0,npol-1 do begin
        Ltmp = L(ipsi(j,i),*)
        ins = inside(ipsi(j,i),*)
        iz = where(ins eq 0,count)
        if count ne 0 then ins(iz(0):npts-2) = 0
        drsqtmp = drsq(ipsi(j,i),*)*ins
        drsqNLtmp = interpol(drsqtmp,Ltmp,Lave(i,*))
        insNLtmp = round(interpol(ins,Ltmp,Lave(i,*)))
        if keyword_set(Lfit) then Lmin(i)=Lfit else Lmin(i) = min([Lmin(i),max(Ltmp*ins)])
        iLz = where(Lave(i,*) gt max(Ltmp*ins),count)
        if count ne 0 then begin
            drsqNLtmp(iLz) = 0.0
            insNLtmp(iLz) = 0.0
        endif
        drsqave(i,*) += drsqNLtmp
        insave(i,*) += insNLtmp
;        oplot,[max(Ltmp*ins),max(Ltmp*ins)],[-1e6,1e6]
    endfor
;    oplot,Lave(i,*),insave(i,*)
    drsqave(i,*) = drsqave(i,*)/insave(i,*)
    junk = min(Lave(i,*)-0.2*Lmin(i),ifit1,/abs)
    junk = min(Lave(i,*)-0.9*Lmin(i),ifit2,/abs)
    xfit = Lave(i,ifit1:ifit2)
    coeff = poly_fit(xfit,drsqave(i,ifit1:ifit2),1,yfit=yfit)
    slope(i) = coeff(1)
    coeff2 = poly_fit(alog10(xfit),alog10(drsqave(i,ifit1:ifit2)),1)
    tdep(i) = coeff2(1)
;    stop
;    slope(i) = mean(drsqave(i,1:NL-1)/Lave(i,1:NL-1))

    if plotit eq i or keyword_set(ps) then begin
        make_nice_plots,ps=ps
        z=mycolors()
        pos = plot_positions(nrow=4,ncol=1)
;        pos(1,3) += .
;        pos(3,3) += .05
        pos(1,0) += -.12
        pos(3,0) += -.12
        pos(1,1) += -.12
        pos(3,1) += -.12
;        pos([1,3],[0,1]) += -0.05
        !p.position = pos(*,0)
        plot,Lave(i,*),drsqave(i,*),thick=0.2,ytitle=textoidl('<(\Delta\rho)^2>'), $
          title=textoidl('\psi_N^{EFIT} = '+string(psi(i),FORMAT='(F5.3)') $
                         +', s^{SIESTA} = '+string(s_siesta(i),FORMAT='(F5.3)')), $
          xtickname=replicate(' ',60),xrange=[Lave(i,0),Lave(i,NL-1)],xstyle=1
;        oplot,Ltmp,drsqtmp
;        oplot,Lave(i,*),drsqNLtmp,psym=1,color=z.brick
        oplot,Lave(i,*),drsqave(i,*),psym=4,color=z.blue
;        oplot,Lave(i,*),slope(i)*Lave(i,*),color=z.green
        oplot,[Lmin(i),Lmin(i)],[-1e6,1e6]

        !p.position = pos(*,1)
        plot,Lave(i,*),insave(i,*),xtitle='L (m)',ytitle='# FL Left',/noerase,$
          xrange=[Lave(i,0),Lave(i,NL-1)],xstyle=1,yrange=[0,npol+5]

        !p.position = pos(*,3)
        plot,Lave(i,1:NL-1),drsqave(i,1:NL-1),ytitle=textoidl('<(\Delta\rho)^2>'),/noerase, $
          /xlog,/ylog,xrange=[Lave(i,1),Lave(i,NL-1)],xstyle=1,xtitle='L (m)',$
          yrange=[1e-3*max(drsqave(i,*)),max(drsqave(i,*))],$
          xticklen=1.0,yticklen=1.0,xgridstyle=1,ygridstyle=1
        oplot,xfit,1.3*yfit,color=z.blue
        oplot,[Lmin(i),Lmin(i)],[1e-16,1e6]
        legst = [['Slope of linear fit: '+strtrim(string(slope(i),format='(E10.3)'),2)],$
                 [textoidl('<(\Delta\rho)^2> ~ L^{'+strtrim(string(tdep(i),format='(F10.2)'),2)+'}')]]
        legend,legst,pos(*,3),color=[z.blue,z.blue],charthick=2.0
        !p.position = 0

        if keyword_set(pdf) then begin
            npdf = n_elements(pdf)
            for k=0,npdf-1 do oplot,[pdf(k),pdf(k)],[1e-16,1e6],color=z.brick,linestyle=2
            pos = plot_positions(nrow=npdf,ncol=1)
            !p.multi=[0,1,3,0]
            make_nice_plots,ps=makeps
            for k=0,npdf-1 do begin
                ik = lonarr(npol)
                for mm = 0,npol-1 do begin 
                    junk = min(L(mm,*)-pdf(k),/abs,iktmp) 
                    ik(mm)=iktmp
                endfor
                nums = lindgen(npol)
                ins = where(inside(nums,ik) and ik lt npts-2)
                junk = min(Lave(i,*)-pdf(k),/abs,ikave)
;                 plot,L(lindgen(npol),ik)
;                plot,L(nums(ins),ik(ins))
                print,'   Mean,           var,            mean(drsq),     drsqave: '
                print, mean(dr(nums(ins),ik(ins))),variance(dr(nums(ins),ik(ins))),mean(dr(nums(ins),ik(ins))^2),drsqave(i,ikave)
                var=variance(dr(nums(ins),ik(ins)))
                sig = stddev(dr(nums(ins),ik(ins)))
                ydata = histogram(dr(nums(ins),ik(ins)),nbins=48,locations=xdata,min=-3*sig,max=3*sig)
                nleft = fix(n_elements(ins))
                binsize = (max(xdata)-min(xdata))/(n_elements(xdata)-1)
                ydata = double(ydata)/binsize/npol
                plot,xdata,ydata,psym=4,xtitle=textoidl('\Delta\rho'),/ylog,yrange=[0.9/binsize/npol,max(ydata)],$
                  title=textoidl('L= '+string(pdf(k),FORMAT='(E10.2)')+', Variance= '+string(var,format='(E10.3)')+$
                                     ', <(\Delta\rho)^2>= '+string(drsqave(i,ikave),format='(E10.3)')),$
                  ytitle='PDF, # Lines= '+strtrim(nleft,2),xcharsize=1.2,ycharsize=1.2
;                stop
            endfor
            !p.multi=0
        endif

    endif

    ik = where(psistart ge psi(i)-dpsi and psistart lt psi(i)+dpsi and inside)
;    ins = where(inside(ik)) gt 0.5
    mdr = dr(ik)
    mdl = dL(ik)
;    window,i
;    plot,mdl,mdr^2,/xlog,/ylog,psym=3
;    print,'Min,max dr^2: ',min(mdr^2),max(mdr^2)
;    print,'Min,max dl: ',min(mdl),max(mdl)
    D(i) = mean(mdr^2/mdl)/2

endfor

Dave = drsqave/(2*Lave)

if keyword_set(ps) then begin
    device,/close
    set_plot,'x'
    z=mycolors()
endif

return,{L:L,drsq:drsq,Lave:Lave,drsqave:drsqave,psi:psi,D:D,slope:slope,Lmin:Lmin,Dave:Dave,dr:dr,in:inside}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_pest_coords,g,ntheta=ntheta,npsi=npsi,noload=noload,nosave=nosave,quiet=quiet,pnwant=pnwant,close=close,psiedge=psiedge

;savname = g.filename + '_pest.sav'
;if file_test(savname) and not keyword_set(noload) then begin
;    restore,savname
;endif else begin

if not keyword_set(ntheta) then ntheta = 128
if not keyword_set(nrefine) then nrefine = 100
if not keyword_set(npsi) then begin
    pn = g.pn
    npsi = g.mw
endif else  pn = dindgen(npsi)/(npsi-1)
if not keyword_set(psiedge) then psiedge=0.999
pn(npsi-1) = psiedge

if keyword_set(pnwant) then begin
    npsi = n_elements(pnwant)
    pn = pnwant
endif

qpsi = interpol(g.qpsi,g.pn,pn)

rmp = build_nstx_rwmcoils(0d3)
rpest = dblarr(npsi,ntheta)
zpest = rpest & phipest=rpest
jac = rpest

f=follow_fieldlines_phi(g,rmp,0.5,0.5,1,1,0.0d,!dpi/180,500*max(g.qpsi),pnarr=pn,quiet=quiet)
nth_fl = n_elements(f.theta(0,*))

if not keyword_set(pnwant) then begin
    rpest(0,*) = g.rmaxis
    zpest(0,*) = g.zmaxis
    istart=1
endif else istart=0

;for i = 200,200 do begin
for i=istart,npsi-1 do begin
   icross=where(f.theta(i,1:nth_fl-1)-f.theta(i,0:nth_fl-2) lt 0.0)
;   print,'first; ',icross(0),i
; JMC 2/25/13: kluged in a test to see if the poloidal field is in the
; right direction
   if f.theta(i,2) lt f.theta(i,1) then begin
      icross=where(f.theta(i,1:nth_fl-1)-f.theta(i,0:nth_fl-2) gt 0.0)
      icross(0) = icross(1)
   endif
;   print,'first; ',icross(0)
;   if icross eq -1 then stop
    th_fl = f.theta(i,0:icross(0)+1)
    th_fl(icross(0)+1) += 2*!dpi
    phi_fl = f.phi(i,0:icross(0)+1)
    phi_th0 = interpol(phi_fl,th_fl,2*!dpi)
    r_fl = f.r(i,0:icross(0)+1)
    z_fl = f.z(i,0:icross(0)+1)
;        plot,phi_fl,th_fl
;        oplot,[phi_th0,phi_th0],[2*!dpi,2*!dpi],psym=5
    if keyword_set(close) then phi_pest = phi_th0*dindgen(ntheta)/(ntheta-1) $
    else phi_pest = phi_th0*dindgen(ntheta)/ntheta
    rpest(i,*) = interpol(r_fl,phi_fl,phi_pest)
    zpest(i,*) = interpol(z_fl,phi_fl,phi_pest)

    b=bfield_geq_bicub(g,rpest(i,*),zpest(i,*))
    bpol = sqrt(b.br^2+b.bz^2)
;        jac(i,*) = abs((interpol(g.qpsi,g.pn,pn(i))/(g.rzero*g.bcentr)))*bpol*(rpest(i,*)^3)
    jac(i,*) = abs(qpsi(i)/(g.rzero*g.bcentr))*bpol*(rpest(i,*)^3)
;        stop
endfor
;    if not keyword_set(nosave) then save,filename = savname
;endelse

return,{pn:pn,r:rpest,z:zpest,jac:jac,q:qpsi}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function build_feminterp_tables_pol,g,j,pnres,dpn=dpn,pnclose=pnclose,nri=nri,nrf=nrf,nro=nro,nrs=nrs,np=np,nt=nt,degree=degree,quiet=quiet,dx=dx,extend=extend
tstart = systime(1)

if not keyword_set(nri) then nri=3
if not keyword_set(nro) then nro=3
if not keyword_set(nrf) then nrf=10
if not keyword_set(nrs) then nrs=3
if not keyword_set(np) then np=10
if not keyword_set(nt) then nt=10
if not keyword_set(dpn) then dpn=0.02
if not keyword_set(pnclose) then pnclose=0.001
if not keyword_set(extend) then extend=0.5

nr = nri+nro+nrf+nrs
nrc = nri+nro+nrf

nr=long(nr)
np=long(np)
nt=long(nt)


rg = dblarr(nr,np)
zg = rg & phig = rg & png=rg

pnc = dblarr(nri+nro+nrf)
pnc(0:nri-1) = (pnres-dpn)*dindgen(nri)/nri
pnc(nri:nri+nrf-1) = 2*dpn*dindgen(nrf)/nrf + pnres-dpn
pnc(nri:nri+nrf/2-1) = (dpn-pnclose)*dindgen(nrf/2)/(nrf/2-1) + pnres-dpn
pnc(nri+nrf/2:nri+nrf-1) = (dpn-pnclose)*dindgen((nrf+1)/2)/(nrf/2) + pnres + pnclose
pnc(nri+nrf:nri+nrf+nro-1) = (0.9999-pnres-dpn)*dindgen(nro)/(nro-1) + pnres+dpn

pest = get_pest_coords(g,pnwant=pnc,ntheta=np,/close)

rg(0:nrc-1,*) = pest.r
zg(0:nrc-1,*) = pest.z
for i=0,nrc-1 do png(i,*) = pnc(i)

pf = get_pest_coords(g,pnwant=0.9999,ntheta=2000,/close)
rsep = pf.r
zsep = pf.z
thsep = (atan(zsep-g.zmaxis,rsep-g.rmaxis)+2*!dpi) mod (2*!dpi)
isort = sort(thsep)
rsep=rsep(isort)
zsep=zsep(isort)
thsep=thsep(isort)
sepobj = obj_new('IDLanROI',rsep,zsep)

rst = pest.r(nrc-1,*)
zst = pest.z(nrc-1,*)
if nrs gt 0 then begin
    for i=0,np-1 do begin
        lmag = sqrt((rst(i)-g.rmaxis)^2 + (zst(i)-g.zmaxis)^2)
        lr = (rst(i) - g.rmaxis)/lmag
        lz = (zst(i) - g.zmaxis)/lmag
        dl = 0.001 + extend*dindgen(nrs)/(nrs-1)
        rg(nrc:nr-1,i) = rst(i) + lr*dl
        zg(nrc:nr-1,i) = zst(i) + lz*dl
        png(nrc:nr-1,i) = 1.0d + dl
    endfor
endif

thg = (atan(zg-g.zmaxis,rg-g.rmaxis)+2*!dpi) mod (2*!dpi)

phi = 2*!dpi*dindgen(nt)/(nt-1)
;print,'Phi end points: ',phi(0),phi(nt-1)

bels = get_b_elements([0,0,0],/derivs,degree=degree)
nelem = n_elements(bels.bx)
f = dblarr(nr*np*nt,nelem)
bels = replicate(bels,8)

x0 = dblarr(nr*np*nt)
y0=x0 & z0=x0

bin = bfield_bs(1.0,0.0,0.0,j.coil,j.current,/derivs)
bin = replicate(bin,8)

for ir=0L,nr-2 do begin
    for ip=0L,np-2 do begin
        for it=0L,nt-2 do begin
            index = ir*np*nt + ip*nt + it
            r0 = 0.25*(rg(ir,ip)+rg(ir+1,ip)+rg(ir,ip+1)+rg(ir+1,ip+1))
            z0(index) = 0.25*(zg(ir,ip)+zg(ir+1,ip)+zg(ir,ip+1)+zg(ir+1,ip+1))
            phi0 = 0.5*(phi(it)+phi(it+1))
            x0(index) = r0*cos(phi0)
            y0(index) = r0*sin(phi0)
            bmod=0.0d
            dermod=0.0d
            for i=0,7 do begin
                ird = i mod 2
                ipd = floor(i/4)
                itd = (i - 4*ipd - ird)/2
                irc=ir+ird
                ipc=ip+ipd
                itc=it+itd
                x = rg(irc,ipc)*cos(phi(itc))
                y = rg(irc,ipc)*sin(phi(itc))
                bin(i) = bfield_bs(x,y,zg(irc,ipc),j.coil,j.current,/derivs)
                bels(i) = get_b_elements([x-x0(index),y-y0(index),zg(irc,ipc)-z0(index)],/derivs,degree=degree)
            endfor
            gs = max(abs([bin.bx,bin.by,bin.bz]))/max(abs([bin.dbxdx,bin.dbxdy,bin.dbxdz,bin.dbydy,bin.dbydz]))
            a=get_ab_fem_scalar(bels,bin,grid_spacing=gs)
            f(index,*)=invert(a.a)##a.b
        endfor
    endfor
endfor
if not keyword_set(quiet) then print,'elapsed time: ',systime(1)-tstart,' sec'

return,{rg:rg,zg:zg,thg:thg,phi:phi,png:png,nr:nr,nrc:nrc,np:np,nt:nt,f:f,x0:x0,y0:y0,z0:z0, $
        rsep:rsep,zsep:zsep,thsep:thsep,sepobj:sepobj,rax:g.rmaxis,zax:g.zmaxis, $
        r2d:g.r,z2d:g.z,pn2d:(g.psirz-g.ssimag)/(g.ssibry-g.ssimag)}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_feminterp_pol,x,y,z,t,derivs=derivs,degree=degree

if (n_elements(size(x,/dim)) gt 1) then begin
    x=transpose(x)
    y=transpose(y)
    z=transpose(z)
endif

;print,size(x,/dim)

r = sqrt(x^2+y^2)
phi = (atan(y,x) + 2*!dpi) mod (2*!dpi)
ihigh = where(phi ge (2*!dpi+t.phi(0)),count)
if count gt 0 then phi(ihigh) -= 2*!dpi
theta = (atan(z-t.zax,r-t.rax) + 2*!dpi) mod (2*!dpi)

ir = ((r-t.r2d(0))/(t.r2d(1)-t.r2d(0)))
iz = ((z-t.z2d(0))/(t.z2d(1)-t.z2d(0)))
pn = interpolate(t.pn2d,ir,iz)
;pn = interpol((g.psirz-g.ssimag)/(g.ssibry-g.ssimag),ir,iz)

it = floor((phi-t.phi(0))/(t.phi(1)-t.phi(0)))
ir = 0*it
ip = ir

rsep = interpol(t.rsep,t.thsep,theta)
zsep = interpol(t.zsep,t.thsep,theta)
dsep = sqrt((rsep-t.rax)^2 + (zsep-t.zax)^2)
dpts = sqrt((r - t.rax)^2 + (z-t.zax)^2)

ic = where(dpts le dsep, count,complement=is,ncomplement=count2)
if count gt 0 then begin
    ir(ic) = floor(interpol(dindgen(t.nrc),t.png(0:t.nrc-1,0),pn(ic)))
;    for i=0,count-1 do ip(ic(i)) = 11
    for i=0,count-1 do ip(ic(i)) = where(t.thg(ir(ic(i)),0:t.np-2) le theta(ic(i)) and t.thg(ir(ic(i)),1:t.np-1) ge theta(ic(i)))
;    for i=0,count-1 do ip(ic(i)) = floor(interpol(dindgen(t.np),t.thg(ir(ic(i)),*),theta(ic(i))))
endif 

;is = where(dpts gt dsep, count)
if count2 gt 0 then begin
;    rsep=interpol(t.rsep,t.thsep,theta(is))
;    zsep = interpol(t.zsep,t.thsep,theta(is))
    ds = sqrt((r(is)-rsep(is))^2 + (z(is)-zsep(is))^2)
    ir(is) = floor(interpol(dindgen(t.nr),t.png(*,0),ds+1.0d))
;    ir(is) = where(ds+1.0d ge t.png(0:t.nr-2) and ds+1.0d lt t.png(1:t.nr-1))
;    ir(is) = floor((ds-t.png(t.nrc-1,0))/(t.png(t.nrc+1,0)-t.png(t.nrc,0)))
    ip(is) = floor(interpol(dindgen(t.np),t.thg(t.nrc,*),theta(is)))
endif

;stop
    
index = ir*t.np*t.nt + ip*t.nt + it

if abs(t.phi(0)) ge 1e-6 then begin
    phig = phi - t.phi(0)
    xg = r*cos(phig)
    yg = r*sin(phig)
endif else begin
    xg = x
    yg = y
endelse

bout = get_b_elements_timesf([[xg-t.x0(index)],[yg-t.y0(index)],[z-t.z0(index)]],t.f(index,*),derivs=derivs,degree=degree)

if abs(t.phi(0)) ge 1e-6 then begin
    br = bout.bx*cos(phi) + bout.by*sin(phi)
    bphi = -bout.bx*sin(phi) + bout.by*cos(phi)

    bout.bx = br*cos(phi)-bphi*sin(phi)
    bout.by = br*sin(phi)+bphi*cos(phi)
endif

return,bout
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_feminterp,x,y,z,t,derivs=derivs,degree=degree

names = tag_names(t)
if strcmp(names(0),'RG',/fold_case) then bout=bfield_feminterp_pol(x,y,z,t,derivs=derivs,degree=degree) else begin
    

    if (n_elements(size(x,/dim)) gt 1) then begin
        x=transpose(x)
        y=transpose(y)
        z=transpose(z)
    endif

;print,size(x,/dim)

    r = sqrt(x^2+y^2)
    phi = (atan(y,x) + 2*!dpi) mod (2*!dpi)
    ihigh = where(phi ge (2*!dpi+t.phi(0)),count)
    if count gt 0 then phi(ihigh) -= 2*!dpi

    ir = floor((r-t.r(0))/(t.r(1)-t.r(0)))
    iz = floor((z-t.z(0))/(t.z(1)-t.z(0)))
    it = floor((phi-t.phi(0))/(t.phi(1)-t.phi(0)))

    index = ir*t.nz*t.nt + iz*t.nt + it

    if abs(t.phi(0)) ge 1e-6 then begin
        phig = phi - t.phi(0)
        xg = r*cos(phig)
        yg = r*sin(phig)
    endif else begin
        xg = x
        yg = y
    endelse

;stop

;print,x,t.x0(index),y,t.y0(index),z,t.z0(index)
;print,ir,iz,it,index

;bels = get_b_elements_vec([[x-t.x0(index)],[y-t.y0(index)],[z-t.z0(index)]],derivs=derivs,degree=degree)
;bout = bfield_feminterp_sum_vec(bels,t.f(index,*),derivs=derivs)
    bout = get_b_elements_timesf([[xg-t.x0(index)],[yg-t.y0(index)],[z-t.z0(index)]],t.f(index,*),derivs=derivs,degree=degree)
;print,t.x0(index),t.y0(index),t.z0(index)
;stop
;stop
    if abs(t.phi(0)) ge 1e-6 then begin
        br = bout.bx*cos(phi) + bout.by*sin(phi)
        bphi = -bout.bx*sin(phi) + bout.by*cos(phi)

        bout.bx = br*cos(phi)-bphi*sin(phi)
        bout.by = br*sin(phi)+bphi*cos(phi)
    endif

endelse 

return,bout
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function bfield_feminterp_cyl,r,phi,z,t,derivs=derivs

x=r*cos(phi)
y=r*sin(phi)
b=bfield_feminterp(x,y,z,t,derivs=derivs)

br = b.bx*cos(phi) + b.by*sin(phi)
bphi = -b.bx*sin(phi) + b.by*cos(phi)
bz = b.bz

if keyword_set(derivs) then $
  return,{br:br,bphi:bphi,bz:bz,$
          dbxdx:b.dbxdx,dbxdy:b.dbxdy,dbxdz:b.dbxdz, $
          dbydx:b.dbxdy,dbydy:b.dbydy,dbydz:b.dbydz, $
          dbzdx:b.dbxdz,dbzdy:b.dbydz,dbzdz:-(b.dbxdx+b.dbydy)} $
  else return,{br:br,bphi:bphi,bz:bz}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test_feminterp_ring,rmp,t,degree=degree

nphi=100

phi = 2*!dpi*dindgen(nphi)/(nphi-1)
;phi = 0.05*dindgen(nphi)/(nphi-1)
r = 1.27 + 0*phi
z = -0.2 + 0*phi

x=r*cos(phi)
y=r*sin(phi)

;bs = bfield_bs(0,0,0,rmp.coil,rmp.current,/derivs)
;bf = bfield_feminterp(r(0),0,z(0),t,degree=degree,/derivs)
;bs = replicate(bs,nphi)
;bf = replicate(bf,nphi)

;for i=0,nphi-1 do begin
;    x = r(i)*cos(phi(i))
;    y = r(i)*sin(phi(i))
;    bs(i) = bfield_bs(x,y,z(i),rmp.coil,rmp.current,/derivs)
;    bf(i) = bfield_feminterp(x,y,z(i),t,degree=degree,/derivs)
;end

bs=bfield_bs(x,y,z,rmp.coil,rmp.current,/derivs)
bf=bfield_feminterp(x,y,z,t,degree=degree,/derivs)

col=mycolors()
pos = plot_positions(nrow=3,ncol=3)
!p.position=pos(*,0)
plot,phi,bs.bx,xtickname=replicate(' ',6),ytitle='Bx'
oplot,phi,bf.bx,color=col.brick

!p.position=pos(*,1)
plot,phi,bs.by,xtickname=replicate(' ',6),ytitle='By',/noerase
oplot,phi,bf.by,color=col.brick

!p.position=pos(*,2)
plot,phi,bs.bz,xtitle='phi',ytitle='Bz',/noerase
oplot,phi,bf.bz,color=col.brick

!p.position=pos(*,3)
plot,phi,bs.dbxdx,xtickname=replicate(' ',6),ytitle='dBxdx',/noerase
oplot,phi,bf.dbxdx,color=col.brick

!p.position=pos(*,4)
plot,phi,bs.dbxdy,xtickname=replicate(' ',6),ytitle='dBxdy',/noerase
oplot,phi,bf.dbxdy,color=col.brick

!p.position=pos(*,5)
plot,phi,bs.dbxdz,xtitle='phi',ytitle='dBxdz',/noerase
oplot,phi,bf.dbxdz,color=col.brick

!p.position=pos(*,6)
plot,phi,bs.dbydy,xtickname=replicate(' ',6),ytitle='dBydy',/noerase
oplot,phi,bf.dbydy,color=col.brick

!p.position=pos(*,7)
plot,phi,bs.dbydz,xtickname=replicate(' ',6),ytitle='dBydz',/noerase
oplot,phi,bf.dbydz,color=col.brick

!p.position=pos(*,8)
plot,phi,bf.dbxdx+bf.dbydy+bf.dbzdz ,xtickname=replicate(' ',6),ytitle='Div, Curl',/noerase,yrange=[-0.00001,0.00001]
oplot,phi,bf.dbydz-bf.dbzdy,color=col.brick
oplot,phi,bf.dbxdz-bf.dbzdx,color=col.blue
oplot,phi,bf.dbxdy-bf.dbydx,color=col.green
legend,['div','dbydz','dbxdz','dbxdy'],pos(*,8),color=[col.black,col.brick,col.blue,col.green]

!p.position=0

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test_feminterp_ring_derivs,rmp,t,degree=degree

nphi=100
ds = 0.00001

phi = 2*!dpi*dindgen(nphi)/(nphi-1)
r = 1.265 + 0*phi
z = -0.2 + 0*phi

bf = bfield_feminterp(r(0),0,z(0),t,degree=degree,/derivs)
bf = replicate(bf,nphi)
;bfxm = replicate(bf,nphi)
;bfyp = replicate(bf,nphi)
;bfym = replicate(bf,nphi)
;bfzp = replicate(bf,nphi)
;bfzm = replicate(bf,nphi)

dbxdx=dblarr(nphi)
dbxdy=dbxdx & dbxdz=dbxdx & dbydx=dbxdx & dbydy=dbxdx & dbydz=dbxdx & dbzdx=dbxdx & dbzdy=dbxdx & dbzdz=dbxdx

for i=0,nphi-1 do begin
    x = r(i)*cos(phi(i))
    y = r(i)*sin(phi(i))
    bf(i) = bfield_feminterp(x,y,z(i),t,degree=degree,/derivs)
    bfxp = bfield_feminterp(x+ds,y,z(i),t,degree=degree,/derivs)
    bfxm = bfield_feminterp(x-ds,y,z(i),t,degree=degree,/derivs)
    bfyp = bfield_feminterp(x,y+ds,z(i),t,degree=degree,/derivs)
    bfym = bfield_feminterp(x,y-ds,z(i),t,degree=degree,/derivs)
    bfzp = bfield_feminterp(x,y,z(i)+ds,t,degree=degree,/derivs)
    bfzm = bfield_feminterp(x,y,z(i)-ds,t,degree=degree,/derivs)
    dbxdx(i) = (bfxp.bx-bfxm.bx)/(2*ds)
    dbydx(i) = (bfxp.by-bfxm.by)/(2*ds)
    dbzdx(i) = (bfxp.bz-bfxm.bz)/(2*ds)
    dbxdy(i) = (bfyp.bx-bfym.bx)/(2*ds)
    dbydy(i) = (bfyp.by-bfym.by)/(2*ds)
    dbzdy(i) = (bfyp.bz-bfym.bz)/(2*ds)
    dbxdz(i) = (bfzp.bx-bfzm.bx)/(2*ds)
    dbydz(i) = (bfzp.by-bfzm.by)/(2*ds)
    dbzdz(i) = (bfzp.bz-bfzm.bz)/(2*ds)
end

col=mycolors()
pos = plot_positions(nrow=3,ncol=3)
!p.position=pos(*,0)
plot,phi,bf.dbxdx,xtickname=replicate(' ',6),ytitle='dbxdx'
oplot,phi,dbxdx,color=col.brick

!p.position=pos(*,1)
plot,phi,bf.dbxdy,xtickname=replicate(' ',6),ytitle='dbxdy',/noerase
oplot,phi,dbxdy,color=col.brick

!p.position=pos(*,2)
plot,phi,bf.dbxdz,xtitle='phi',ytitle='dbxdz',/noerase
oplot,phi,dbxdz,color=col.brick

!p.position=pos(*,3)
plot,phi,bf.dbydx,xtickname=replicate(' ',6),ytitle='dbydx',/noerase
oplot,phi,dbydx,color=col.brick

!p.position=pos(*,4)
plot,phi,bf.dbydy,xtickname=replicate(' ',6),ytitle='dbydy',/noerase
oplot,phi,dbydy,color=col.brick

!p.position=pos(*,5)
plot,phi,bf.dbydz,xtitle='phi',ytitle='dbydz',/noerase
oplot,phi,dbydz,color=col.brick

!p.position=pos(*,6)
plot,phi,bf.dbzdx,xtickname=replicate(' ',6),ytitle='dbzdx',/noerase
oplot,phi,dbzdx,color=col.brick

!p.position=pos(*,7)
plot,phi,bf.dbzdy,xtickname=replicate(' ',6),ytitle='dbzdy',/noerase
oplot,phi,dbzdy,color=col.brick

!p.position=pos(*,8)
plot,phi,bf.dbzdz,xtickname=replicate(' ',6),ytitle='dbzdz',/noerase
oplot,phi,dbzdz,color=col.brick
;plot,phi,bf.dbxdx+bf.dbydy+bf.dbzdz ,xtickname=replicate(' ',6),ytitle='Div, Curl',/noerase,yrange=[-0.001,0.001]
;oplot,phi,bf.dbydz-bf.dbzdy,color=col.brick
;oplot,phi,bf.dbxdz-bf.dbzdx,color=col.blue
;oplot,phi,bf.dbxdy-bf.dbydx,color=col.green

!p.position=0

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test_feminterp_pol,rmp,t,degree=degree

nth=100
theta = 2*!dpi*dindgen(nth)/(nth-1)
r=1.0 + 0.5*cos(theta)
z=0.0 + 1.0*sin(theta)
phi = 0.0 + 0*r

bs = bfield_bs(0,0,0,rmp.coil,rmp.current,/derivs)
bf = bfield_feminterp(r(0),0,z(0),t,degree=degree,/derivs)
bs = replicate(bs,nth)
bf = replicate(bf,nth)

for i=0,nth-1 do begin
    x = r(i)*cos(phi(i))
    y = r(i)*sin(phi(i))
    bs(i) = bfield_bs(x,y,z(i),rmp.coil,rmp.current,/derivs)
    bf(i) = bfield_feminterp(x,y,z(i),t,degree=degree,/derivs)
end

col=mycolors()
pos = plot_positions(nrow=3,ncol=3)
!p.position=pos(*,0)
plot,theta,bs.bx,xtickname=replicate(' ',6),ytitle='Bx'
oplot,theta,bf.bx,color=col.brick

!p.position=pos(*,1)
plot,theta,bs.by,xtickname=replicate(' ',6),ytitle='By',/noerase
oplot,theta,bf.by,color=col.brick

!p.position=pos(*,2)
plot,theta,bs.bz,xtitle='theta',ytitle='Bz',/noerase
oplot,theta,bf.bz,color=col.brick

!p.position=pos(*,3)
plot,theta,bs.dbxdx,xtickname=replicate(' ',6),ytitle='dBxdx',/noerase
oplot,theta,bf.dbxdx,color=col.brick

!p.position=pos(*,4)
plot,theta,bs.dbxdy,xtickname=replicate(' ',6),ytitle='dBxdy',/noerase
oplot,theta,bf.dbxdy,color=col.brick

!p.position=pos(*,5)
plot,theta,bs.dbxdz,xtitle='theta',ytitle='dBxdz',/noerase
oplot,theta,bf.dbxdz,color=col.brick

!p.position=pos(*,6)
plot,theta,bs.dbydy,xtickname=replicate(' ',6),ytitle='dBydy',/noerase
oplot,theta,bf.dbydy,color=col.brick

!p.position=pos(*,7)
plot,theta,bs.dbydz,xtickname=replicate(' ',6),ytitle='dBydz',/noerase
oplot,theta,bf.dbydz,color=col.brick

!p.position=pos(*,8)
plot,theta,bf.dbxdx+bf.dbydy+bf.dbzdz ,xtickname=replicate(' ',6),ytitle='Div, Curl',/noerase,yrange=[-0.001,0.001]
oplot,theta,bf.dbydz-bf.dbzdy,color=col.brick
oplot,theta,bf.dbxdz-bf.dbzdx,color=col.blue
oplot,theta,bf.dbxdy-bf.dbydx,color=col.green

!p.position=0

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test_feminterp_pn,g,pn,rmp,t,degree=degree,sol=sol,xrange=xrange

nth=1000

if not keyword_set(sol) then begin
    p=get_pest_coords(g,pnwant=pn,ntheta=nth)
    r=p.r
    z=p.z
endif else begin
    r0=build_nstx_rwmcoils(0d3)
    f=follow_fieldlines_phi(g,r0,sol,sol,1,1,0.0d,!dpi/nth,nth)
    r = f.r
    z = f.z
endelse

theta = 2*!dpi*dindgen(nth)/(nth-1)

phi = 0.0 + 0*r

bs = bfield_bs(0,0,0,rmp.coil,rmp.current,/derivs)
bf = bfield_feminterp(r(0),0,z(0),t,degree=degree,/derivs)
bs = replicate(bs,nth)
bf = replicate(bf,nth)

for i=0,nth-1 do begin
    x = r(i)*cos(phi(i))
    y = r(i)*sin(phi(i))
    bs(i) = bfield_bs(x,y,z(i),rmp.coil,rmp.current,/derivs)
    bf(i) = bfield_feminterp(x,y,z(i),t,degree=degree,/derivs)
end

col=mycolors()
pos = plot_positions(nrow=3,ncol=3)
!p.position=pos(*,0)
plot,theta,bs.bx,xtickname=replicate(' ',6),ytitle='Bx',xrange=xrange
oplot,theta,bf.bx,color=col.brick

!p.position=pos(*,1)
plot,theta,bs.by,xtickname=replicate(' ',6),ytitle='By',/noerase,xrange=xrange
oplot,theta,bf.by,color=col.brick

!p.position=pos(*,2)
plot,theta,bs.bz,xtitle='theta',ytitle='Bz',/noerase,xrange=xrange
oplot,theta,bf.bz,color=col.brick

!p.position=pos(*,3)
plot,theta,bs.dbxdx,xtickname=replicate(' ',6),ytitle='dBxdx',/noerase,xrange=xrange
oplot,theta,bf.dbxdx,color=col.brick

!p.position=pos(*,4)
plot,theta,bs.dbxdy,xtickname=replicate(' ',6),ytitle='dBxdy',/noerase,xrange=xrange
oplot,theta,bf.dbxdy,color=col.brick

!p.position=pos(*,5)
plot,theta,bs.dbxdz,xtitle='theta',ytitle='dBxdz',/noerase,xrange=xrange
oplot,theta,bf.dbxdz,color=col.brick

!p.position=pos(*,6)
plot,theta,bs.dbydy,xtickname=replicate(' ',6),ytitle='dBydy',/noerase,xrange=xrange
oplot,theta,bf.dbydy,color=col.brick

!p.position=pos(*,7)
plot,theta,bs.dbydz,xtickname=replicate(' ',6),ytitle='dBydz',/noerase,xrange=xrange
oplot,theta,bf.dbydz,color=col.brick

!p.position=pos(*,8)
plot,theta,bf.dbxdx+bf.dbydy+bf.dbzdz ,xtickname=replicate(' ',6),ytitle='Div, Curl',/noerase,yrange=[-0.001,0.001],xrange=xrange
oplot,theta,bf.dbydz-bf.dbzdy,color=col.brick
oplot,theta,bf.dbxdz-bf.dbzdx,color=col.blue
oplot,theta,bf.dbxdy-bf.dbydx,color=col.green

!p.position=0

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_pest_coords_old,g,ntheta=ntheta,nrefine=nrefine,noload=noload,nosave=nosave

savname = g.filename + '_pest.sav'
if file_test(savname) and not keyword_set(noload) then begin
    restore,savname
endif else begin

    if not keyword_set(ntheta) then ntheta = 128
    if not keyword_set(nrefine) then nrefine = 100

    rmp = build_nstx_rwmcoils(0d3)
    rpest = dblarr(g.mw,ntheta+1)
    zpest = rpest & phipest=rpest

    for i = 0,g.mw-1 do begin
        f=follow_fieldlines_phi(g,rmp,g.pn(i),g.pn(i),1,1,0.d,2.2*!dpi*g.qpsi(i)/(nrefine*ntheta),nrefine*ntheta)
        if f.z(1) gt f.z(0) then icross = where(f.z(2:nrefine*ntheta) ge f.z(0) and f.z(1:nrefine*ntheta-1) lt f.z(0))
        ind1 = max([icross-10,0])
        ind2 = min([icross+10,nrefine*ntheta])
        phimax = interpol(f.phi(ind1:ind2),f.z(ind1:ind2),f.z(0))
        fpest = follow_fieldlines_phi(g,rmp,g.pn(i),g.pn(i),1,1,0.0,phimax/(nrefine*ntheta),nrefine*ntheta)
        rpest(i,*) = fpest.r(0:nrefine*ntheta:nrefine)
        zpest(i,*) = fpest.z(0:nrefine*ntheta:nrefine)
        phipest(i,*) = fpest.phi(0:nrefine*ntheta:nrefine)
    endfor
    if not keyword_set(nosave) then save,filename = savname
endelse

return,{r:rpest,z:zpest,phi:phipest}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_b_spectrum,g,rmp,mmax=mmax,nn=nn,nsurf=nsurf,ntheta=ntheta,nphi=nphi,nojac=nojac,plot=plot,yrange=yrange,makeps=makeps,$
                        brmax=brmax,mplot=mplot,pnwant=pnwant,rmp2=rmp2,fem=fem,secfem=secfem,psiedge=psiedge,nardon=nardon,m3dc1=m3dc1
  tstart = systime(1)

  if keyword_set(makeps) then begin
     set_plot,'ps'
     !p.font=0
     device,file='b_spectrum.ps',/color,bits_per_pixel=8,/inches,xsize=6,ysize=5,xoffset=.5,yoffset=.5,/helvetica,/bold
  endif


  if not keyword_set(nsurf) then nsurf=2
  if not keyword_set(mmax) then mmax=60
  if not keyword_set(nn) then nn=3
  if not keyword_set(ntheta) then ntheta=2*mmax
  if not keyword_set(nphi) then nphi=8*nn
  if keyword_set(pnwant) then nsurf=n_elements(pnwant)

  phi = 2*!dpi*dindgen(nphi)/nphi
  theta = 2*!dpi*dindgen(ntheta)/ntheta
  dphi = phi(1)-phi(0)
  dth = theta(1)-theta(0)

  marr = dindgen(2*mmax+1) - mmax

  p=get_pest_coords(g,npsi=nsurf,/quiet,ntheta=ntheta,pnwant=pnwant,psiedge=psiedge)

  dsdr = 0*p.r
  dsdz= dsdr
  psin = dsdr
  bnorm = dblarr(ntheta,nphi)
  jac = dsdr

  br_c = dblarr(nsurf,2*mmax+1)
  br_s = br_c
  S = dblarr(nsurf)

  for i=0,nsurf-1 do begin
     psi = get_psi_bicub(g,p.r(i,*),p.z(i,*),/deriv)
     dsdr(i,*)=psi.dsdr
     dsdz(i,*)=psi.dsdz
     psin(i,*) = (-psi.psi-g.ssimag)/(g.ssibry-g.ssimag)
  endfor
  psin12 = sqrt(psin)
  ds12dr = ((-dsdr-g.ssimag)/(g.ssibry-g.ssimag))/(2*psin12)
  ds12dz = ((-dsdz-g.ssimag)/(g.ssibry-g.ssimag))/(2*psin12)
  ds12ave = sqrt(ds12dr^2 + ds12dz^2)
  rn = dsdr/sqrt(dsdr^2+dsdz^2)
  zn = dsdz/sqrt(dsdr^2+dsdz^2)

  if keyword_set(pnwant) then istart=0 else istart=1

  for i=istart,nsurf-1 do begin
     S(i) = total(p.jac(i,*),2)*2*!dpi*dth
;    b=bfield_geq_bicub(g,p.r(i,*),p.z(i,*))

     if not keyword_set(nardon) then begin
        for k=0,nphi-1 do begin
           if keyword_set(fem) then b=bfield_feminterp_cyl(p.r(i,*),phi(k)+0.0*p.r(i,*),p.z(i,*),fem) $
           else if keyword_set(m3dc1) then b=bfield_m3dc1(p.r(i,*),phi(k)+0.0*p.r(i,*),p.z(i,*),m3dc1) $
           else b=bfield_bs_cyl(p.r(i,*),phi(k)+0.0*p.r(i,*),p.z(i,*),rmp.coil,rmp.current)
;            help,b.br
           bnorm(*,k) = b.br*rn(i,*)+b.bz*zn(i,*)
           if keyword_set(rmp2) then begin
              b2=bfield_bs_cyl(p.r(i,*),phi(k)+0.0*p.r(i,*),p.z(i,*),rmp2.coil,rmp2.current)   
              bnorm(*,k) += b2.br*rn(i,*)+b2.bz*zn(i,*)
           endif
           if keyword_set(secfem) then begin
              b2=bfield_feminterp_cyl(p.r(i,*),phi(k)+0.0*p.r(i,*),p.z(i,*),secfem)   
              bnorm(*,k) += b2.br*rn(i,*)+b2.bz*zn(i,*)
           endif
           if not keyword_set(nojac) then bnorm(*,k) *= p.jac(i,*)
;            window,k
;            col=mycolors()
;            plot,theta,bnorm(*,k)
;            oplot,theta,bnorm(*,k),color=col.brick
;            stop
        endfor
     endif else begin
        for k=0,nphi-1 do begin
           if keyword_set(fem) then b=bfield_feminterp_cyl(p.r(i,*),phi(k)+0.0*p.r(i,*),p.z(i,*),fem) $
           else b=bfield_bs_cyl(p.r(i,*),phi(k)+0.0*p.r(i,*),p.z(i,*),rmp.coil,rmp.current)
           bax = bfield_geq_bicub(g,p.r(i,*),p.z(i,*))
           bnorm(*,k) = (b.br*ds12dr(i,*) + b.bz*ds12dz(i,*))/(bax.bphi/p.r(i,*))
        endfor
     endelse
        
     for m=0,2*mmax do begin
        for k=0,nphi-1 do begin
           br_c(i,m) += 2*dth*dphi*total(bnorm(*,k)*cos(nn*phi(k)-marr(m)*theta))/S(i)
           br_s(i,m) += 2*dth*dphi*total(bnorm(*,k)*sin(nn*phi(k)-marr(m)*theta))/S(i)
        endfor
     endfor
     if keyword_set(nardon) then begin
        dsave = mean(ds12ave(i,*))
        br_c(i,*) = 2*(br_c(i,*))/(g.rzero*dsave)
        br_s(i,*) = 2*(br_s(i,*))/(g.rzero*dsave)
     endif
  endfor

  br = sqrt(br_c^2 + br_s^2)

  if not keyword_set(pnwant) then begin

     mlow = ceil(nn*min(p.q))
     mhigh = floor(nn*maX(p.q))

     if mlow lt mhigh then mres = dindgen(mhigh-mlow)+ mlow else mres=mlow
     nres = n_elements(mres)
     wid = dblarr(nres)
     pnres=wid
     bres = wid

     if nres gt 1 then shear = deriv(p.pn,p.q) else shear=interpol(deriv(g.pn,g.qpsi),g.pn,p.pn)
     dpsi = abs(g.ssibry-g.ssimag)

     for i=0,nres-1 do begin
        qres = mres(i)/nn
        pnres(i) = interpol(p.pn,p.q,qres)
        sres = interpol(shear,p.pn,pnres(i))
        SAres = interpol(S,p.pn,pnres(i))
        mind = where(marr eq mres(i),count)
        if count gt 0 then begin
           bres(i) = interpol(br(*,mind),p.pn,pnres(i))
           wid(i) = sqrt((16.0/(mres(i)*dpsi))*(qres/sres)*(SAres/(4*!dpi*!dpi))*bres(i))
        endif 
;    print,'Mres,qres,pnres(i),wid= ',mres(i),interpol(p.q,p.pn,pnres(i)),pnres(i),wid(i)
     endfor


     pnchir = 0.5*(pnres(1:nres-1) + pnres(0:nres-2))
     dpchir = pnres(1:nres-1) - pnres(0:nres-2)
     dwchir = 0.5*(wid(1:nres-1) + wid(0:nres-2))

     chir = dwchir/dpchir
  endif else begin
     wid=0
     chir=0
     mres=0
     pnres=0
  endelse

  if keyword_set(plot) then begin
     z=mycolors(/badgreen,/load)

     if not keyword_set(makeps) then window,0
     cmin =0.0
     if keyword_set(yrange) then begin
        junk=min(p.pn-yrange(0),/abs,iy)
        cmax = max(1e4*br(iy:nsurf-1,*))
     endif else cmax = max(1e4*br)
     if keyword_set(brmax) then cmax=brmax

     levels = cmax*dindgen(50)/49

     if not keyword_set(pnwant) then begin
        !p.position = [0.1,0.1,0.80,0.90]
        contour,1e4*transpose(br),marr,p.pn,levels=levels,/fill,yrange=yrange,xtitle='m',ytitle=textoidl('\psi_N'),xrange=mplot,xstyle=1
;    contour,1e4*transpose(br),marr,p.pn,nlevels=50,/fill,yrange=yrange,xtitle='m',ytitle=textoidl('\psi_N')
        oplot,nn*p.q,p.pn,color=z.white,thick=8,linestyle=2
        colorbar,/vertical,minrange=0.0,maxrange=cmax,title=textoidl('b_r^{n='+strtrim(nn,2)+'} (G)'),charsize=1.1,format='(F0.1)'
     endif

     !p.position=0
     if not keyword_set(pnwant) then begin
        if not keyword_set(makeps) then window,1
        plot,p.pn,p.q,xtitle=textoidl('\psi_N'),ytitle='q',xrange=yrange
        for i=0,nres-1 do begin
           qtmp = interpol(p.q,p.pn,pnres(i))
           oplot,[pnres(i)-0.5*wid(i),pnres(i)+0.5*wid(i)],[qtmp,qtmp],color=z.brick,thick=5
        endfor

        if not keyword_set(makeps) then window,2
        plot,pnchir,chir,psym=5,xtitle=textoidl('\psi_N'),ytitle=textoidl('\sigma_{chir}'),xrange=yrange
        oplot,pnchir,0*pnchir+1.0,linestyle=2,color=z.blue
     endif

     if not keyword_set(makeps) then window,3
     !p.position = [0.1,0.1,0.80,0.90]
     contour,1e4*transpose(bnorm),phi/!dpi,theta/!dpi,nlevels=50,/fill,xtitle=textoidl('\phi'),ytitle=textoidl('\theta')
;    colorbar,/vertical,minrange=0.0,maxrange=cmax,title=textoidl('b_r^{n='+strtrim(nn,2)+'} (G)'),charsize=1.1,format='(F0.1)'
     !p.position=0
     if keyword_set(makeps) then begin
        device,/close
        set_plot,'x'
     endif 
  endif


  print,'elapsed time: ',systime(1)-tstart,' sec'

  return,{m:marr,n:nn,br:br,br_c:br_c,br_s:br_s,pn:p.pn,pnres:pnres,wid:wid,mres:mres,chir:chir,bres:bres}
end

pro compare_b_spectrum,a,b,ps=ps,xrange=xrange,brange=brange,ratrange=ratrange,legstr=legstr,mrange=mrange

if keyword_set(ps) then begin
   set_plot,'ps'
   !p.font=0
   device,file='compare_spectra.ps',/color,/inches,xsize=7,ysize=7,xoff=.5,yoff=.5,/helvetica,/bold
endif

z=mycolors()
make_nice_plots,ps=ps

pos = plot_positions(nrow=2,ncol=2)
!p.position=pos(*,0)
if not keyword_set(brange) then brange=[0,max([a.bres*1e4,b.bres*1e4])]
plot,a.pnres,a.bres*1e4,thick=4,xrange=xrange,yrange=brange,xtickname=replicate(' ',60),ytitle=textoidl('b_r')
oplot,b.pnres,b.bres*1e4,thick=4,color=z.brick
if keyword_set(legstr) then legend,legstr,pos(*,0)

!p.position=pos(*,1)
plot,a.pnres,b.bres/a.bres,thick=4,xtitle=textoidl('\psi_N'),ytitle=textoidl('f_{scr}'),/noerase
oplot,[-1e30,1e30],[1.0,1.0],linestyle=2,thick=2

!p.position=pos(*,2)
plot,a.mres,a.bres*1e4,thick=4,xrange=mrange,yrange=brange,xtickname=replicate(' ',60),ytitle=textoidl('b_r'),/noerase
oplot,b.mres,b.bres*1e4,thick=4,color=z.brick
if keyword_set(legstr) then legend,legstr,pos(*,2)

!p.position=pos(*,3)
plot,a.mres,b.bres/a.bres,thick=4,xtitle=textoidl('m'),ytitle=textoidl('f_{scr}'),/noerase
oplot,[-1e30,1e30],[1.0,1.0],linestyle=2,thick=2

if keyword_set(ps) then begin
   device,/close
   set_plot,'x'
endif

!p.position=0

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function find_pnrat,g,mm,nn,dpn=dpn

if not keyword_set(dpn) then dpn=0.01

pn0 = interpol(g.pn,g.qpsi,float(mm)/nn)

rmp=build_nstx_rwmcoils(0d3)

f=follow_fieldlines_phi(g,rmp,pn0-dpn,pn0+dpn,500,1,0.0d,!dpi/180,360*mm,/quiet)

pn = f.psin(*,0)
z0 = f.z(*,0)
z1 = f.z(*,360*mm)

pnres = interpol(pn,z1-z0,0.0d)

;f2=follow_fieldlines_phi(g,rmp,pnres,pnres,1,1,0.0d,!dpi/180,360*mm)

;stop

return,pnres
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_jparrat,g,mm,nn,phidec=phidec,thmult=thmult

if not keyword_set(thmult) then thmult=2
if not keyword_set(phidec) then phidec=5


nth = thmult*mm
nphi = 360L
nphi1 = nphi+1

pn = find_pnrat(g,mm,nn)
pest = get_pest_coords(g,pnwant=pn,ntheta=nth)

nfine = 500
ptmp = get_pest_coords(g,pnwant=pn,ntheta=nfine)
thfine = 2*!dpi*dindgen(nfine)/nfine

dr = shift(ptmp.r,0,-1)-shift(ptmp.r,0,1)
dz = shift(ptmp.z,0,-1)-shift(ptmp.z,0,1)

dfine = sqrt(dr^2+dz^2)

dfine = dfine/mean(dfine)

; b = bfield_geq_bicub(g,ptmp.r,ptmp.z)
; bfine = sqrt(b.br^2+b.bz^2+b.bphi^2)

; alffine = 2*!dpi*dindgen(nfine)/nfine

; Ifine = alffine*dfine*bfine

th0 = 2*!dpi*dindgen(nth)/nth

if thmult gt 2 then phshift=(!dpi/thmult)/mm else phshift=0.0d
;phshift=0.0d

alf0 = cos(mm*(th0+phshift))

dx0 = interpol(dfine,thfine,th0)

b = bfield_geq_bicub(g,pest.r,pest.z)
b0 = sqrt(b.br^2+b.bz^2+b.bphi^2)

;I0 = alf0*dx0*b0

I0 = alf0

print,I0
 
;arc = [0,total(dfine,/cumulative)]
;arcmax = max(arc)

;arceven = arcmax*dindgen(nth)/

;stop

rmp=build_nstx_rwmcoils(0d3)
f=follow_fieldlines_phi(g,rmp,0.5,0.5,nth,1,0.0d,!dpi/180,nphi,/quiet,rstart=pest.r,zstart=pest.z)

r1s = shift(f.r(*,0),0)
z1s = shift(f.z(*,0),0)

;print,f.r(*,nphi)
;print,shift(r1s,-thmult*nn)
;print,f.r(*,nphi)-shift(r1s,-thmult*nn)

;print,f.z(*,nphi)
;print,shift(z1s,-thmult*nn)
;print,f.z(*,nphi)-shift(z1s,-thmult*nn)

;stop

f.r(*,nphi) = shift(r1s,-thmult*nn)
f.z(*,nphi) = shift(z1s,-thmult*nn)

x = f.r*cos(f.phi)
y = f.r*sin(f.phi)
z = f.z

if phidec gt 1 then begin
    xend=x(*,nphi)
    yend=y(*,nphi)
    zend=z(*,nphi)
    xtmp=x(*,0:nphi:phidec)
    ytmp=y(*,0:nphi:phidec)
    ztmp=z(*,0:nphi:phidec)
    nphi1 = n_elements(xtmp(0,*))
    if total(abs(xtmp(*,nphi1-1)-xend)) gt 1d-3 then begin
        print,'Having to fix this up'
        nphi1+=1
        x = dblarr(nth,nphi1)
        y = x & z=x
        x(*,0:nphi1-2)=xtmp
        y(*,0:nphi1-2)=ytmp
        z(*,0:nphi1-2)=ztmp
        x(*,nphi1-1) = xend
        y(*,nphi1-1) = yend
        z(*,nphi1-1) = zend       
    endif else begin
        x=xtmp
        y=ytmp
        z=ztmp
    endelse
endif
    

coil = dblarr(nth*nphi1,3)
current = dblarr(nth*nphi1)
for i=0,nth-1 do begin
    coil(i*nphi1:(i+1)*nphi1-1,0) = x(i,*)
    coil(i*nphi1:(i+1)*nphi1-1,1) = y(i,*)
    coil(i*nphi1:(i+1)*nphi1-1,2) = z(i,*)
    b=bfield_geq_bicub(g,sqrt(x(i,*)^2+y(i,*)^2),z(i,*))
    bref = sqrt(b.br^2+b.bz^2+b.bphi^2)
;    current(i*nphi1:(i+1)*nphi1-1) = 1.0d3*alf0(i);*bref
    current(i*nphi1:(i+1)*nphi1-1) = 1.0d3*I0(i)
    current((i+1)*nphi1-1) = 0.0d
endfor
;stop
return,{coil:coil,current:current,pn:pn}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_jparrat_singleclose,g,mm,nn,phidec=phidec,thmult=thmult

if not keyword_set(thmult) then thmult=2
if not keyword_set(phidec) then phidec=5


nth = thmult*mm
nphi = 360L*mm
nphi1 = nphi+1

nfine = 500
thfine = 2*!dpi*dindgen(nfine)/nfine

pn = find_pnrat(g,mm,nn)
pest = get_pest_coords(g,pnwant=pn,ntheta=nth)

th0 = 2*!dpi*dindgen(nth)/nth
alf0 = cos(mm*th0)

pfine = get_pest_coords(g,pnwant=pn,ntheta=nfine)

;dr = shift(pfine.r,0,-1)-shift(pfine.r,0,1)
;dz = shift(pfine.z,0,-1)-shift(pfine.z,0,1)

dr = pfine.r(1:nfine-1) - pfine.r(0:nfine-2)
dz = pfine.z(1:nfine-1) - pfine.z(0:nfine-2)

dfine = reform(sqrt(dr^2+dz^2))

arc = [0,total(dfine,/cumulative)]
arcmax = max(arc)

arceven = arcmax*dindgen(nth)/nth

theven = interpol(thfine,arc,arceven)

reven = interpol(pfine.r,thfine,theven)
zeven = interpol(pfine.z,thfine,theven)

z=mycolors()
plot,pfine.r,pfine.z,psym=1,/iso
oplot,reven,zeven,psym=5,color=z.brick,thick=5.0

;stop

rmp=build_nstx_rwmcoils(0d3)
f=follow_fieldlines_phi(g,rmp,0.5,0.5,nth,1,0.0d,!dpi/180,nphi,/quiet,rstart=pest.r,zstart=pest.z)

r1s = shift(f.r(*,0),0)
z1s = shift(f.z(*,0),0)

;print,f.r(*,nphi)
;print,r1s

;print,f.z(*,nphi)
;print,z1s

f.r(*,nphi) = r1s
f.z(*,nphi) = z1s

x = f.r*cos(f.phi)
y = f.r*sin(f.phi)
z = f.z

if phidec gt 1 then begin
    xend=x(*,nphi)
    yend=y(*,nphi)
    zend=z(*,nphi)
    xtmp=x(*,0:nphi:phidec)
    ytmp=y(*,0:nphi:phidec)
    ztmp=z(*,0:nphi:phidec)
    nphi1 = n_elements(xtmp(0,*))
    if total(abs(xtmp(*,nphi1-1)-xend)) gt 1d-3 then begin
        print,'Having to fix this up'
        nphi1+=1
        x = dblarr(nth,nphi1)
        y = x & z=x
        x(*,0:nphi1-2)=xtmp
        y(*,0:nphi1-2)=ytmp
        z(*,0:nphi1-2)=ztmp
        x(*,nphi1-1) = xend
        y(*,nphi1-1) = yend
        z(*,nphi1-1) = zend       
    endif else begin
        x=xtmp
        y=ytmp
        z=ztmp
    endelse
endif
    

coil = dblarr(nth*nphi1,3)
current = dblarr(nth*nphi1)
for i=0,nth-1 do begin
    coil(i*nphi1:(i+1)*nphi1-1,0) = x(i,*)
    coil(i*nphi1:(i+1)*nphi1-1,1) = y(i,*)
    coil(i*nphi1:(i+1)*nphi1-1,2) = z(i,*)
;    b=bfield_geq_bicub(g,sqrt(x(i,*)^2+y(i,*)^2),z(i,*))
;    bref = sqrt(b.br^2+b.bz^2+b.bphi^2)
    current(i*nphi1:(i+1)*nphi1-1) = 1.0d3*alf0(i)
    current((i+1)*nphi1-1) = 0.0d
endfor
;stop
return,{coil:coil,current:current,pn:pn}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function scale_jparrat_fem,g,trmp,tjpar,mm,nn,maxerr=maxerr,maxiter=maxiter,nphi=nphi

if not keyword_set(maxerr) then maxerr= 1d-3
if not keyword_set(maxiter) then maxiter=10
if not keyword_set(nphi) then nphi=256

rmp=build_nstx_rwmcoils(0d3)

pres = find_pnrat(g,mm,nn)

sv = get_b_spectrum(g,rmp,pnwant=pres+0.001,fem=trmp,ntheta=2000,nphi=nphi)

im = where(sv.m eq mm)
brv = sv.br(im)

phv = atan(sv.br_s(im),sv.br_c(im))

tout = tjpar

err = 100.d
iter = 0.d

so = get_b_spectrum(g,rmp,pnwant=pres+0.001,fem=tout,ntheta=2000,nphi=nphi)

bro = so.br(im)
pho = atan(so.br_s(im),so.br_c(im))

mult = brv/bro
dphi = phv-!dpi - pho

err = sqrt((so.br_c(im)+sv.br_c(im))^2 + (so.br_s(im)+sv.br_s(im))^2)/bro

print,'Amp= ',brv,bro
print,'Phase= ',phv,pho+!dpi
print,'Error= ',err

;mult = 1.0d
;dphi = !dpi


while err gt maxerr and iter lt maxiter do begin

    tout.f = tout.f*mult(0)
    tout.phi = tout.phi + dphi(0)/nn

    so = get_b_spectrum(g,rmp,pnwant=pres+0.001,fem=tout,ntheta=2000,nphi=nphi)

    bro = so.br(im)
    pho = atan(so.br_s(im),so.br_c(im))

    mult = brv/bro
    dphi = phv-!dpi - pho

    err = sqrt((so.br_c(im)+sv.br_c(im))^2 + (so.br_s(im)+sv.br_s(im))^2)/bro
    iter += 1

    print,'Amp= ',brv,bro
    print,'Phase= ',phv,pho+!dpi
    print,'Error= ',err

endwhile

return,tout
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function change_jpar_phase,j,dphi,nn
jout=j
phic = atan(jout.coil(*,1),jout.coil(*,0))
phinew = phic+dphi(00)/nn
rc = sqrt(jout.coil(*,0)^2 + jout.coil(*,1)^2)
jout.coil(*,0) = rc*cos(phinew)
jout.coil(*,1) = rc*sin(phinew)
return,jout
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_jpar_phase,g,j,mm,nn,ntheta=ntheta,nphi=nphi

;if not keyword_set(nphi) then nphi=129
if not keyword_set(ntheta) then ntheta=2000

pres = find_pnrat(g,mm,nn)

sv = get_b_spectrum(g,j,pnwant=pres+0.001,ntheta=ntheta,nphi=nphi)

im = where(sv.m eq mm)
brv = sv.br(im)
phv = atan(sv.br_s(im),sv.br_c(im))

return,phv
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function check_jpar_phase,g,j,mm,nn,ntheta=ntheta,nphi=nphi
ndphi = 10

dphi = 0.1*dindgen(ndphi)/(ndphi-1)

phase = 0*dphi
for i=0,ndphi-1 do begin
   jtmp=change_jpar_phase(j,dphi(i),3)
   phase(i)=get_jpar_phase(g,jtmp,mm,nn,ntheta=ntheta,nphi=nphi)
endfor

return,{dphi:dphi,phase:phase}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function scale_jparrat_bs,g,rmp,jpar,mm,nn,maxerr=maxerr,maxiter=maxiter,nphi=nphi,ntheta=ntheta

if not keyword_set(maxerr) then maxerr= 1d-3
if not keyword_set(maxiter) then maxiter=10
if not keyword_set(nphi) then nphi=129
if not keyword_set(ntheta) then ntheta=2000

pres = find_pnrat(g,mm,nn)

sv = get_b_spectrum(g,rmp,pnwant=pres+0.001,ntheta=ntheta,nphi=nphi)

im = where(sv.m eq mm)
brv = sv.br(im)
phv = atan(sv.br_s(im),sv.br_c(im))

jout = jpar

so = get_b_spectrum(g,jout,pnwant=pres+0.001,ntheta=ntheta,nphi=nphi)
bro = so.br(im)
pho = atan(so.br_s(im),so.br_c(im))

mult = brv/bro
dphi = phv-!dpi - pho

err = sqrt((so.br_c(im)+sv.br_c(im))^2 + (so.br_s(im)+sv.br_s(im))^2)/bro

print,'Amp= ',brv,bro
print,'Phase= ',phv,pho+!dpi
print,'Error= ',err

;mult = 1.0d
;dphi = 0.0d
iter = 0.d

while err gt maxerr and iter lt maxiter do begin

    jout.current = jout.current*mult(0)
    phic = atan(jout.coil(*,1),jout.coil(*,0))
    phinew = phic+dphi(0)/nn
    rc = sqrt(jout.coil(*,0)^2 + jout.coil(*,1)^2)
    jout.coil(*,0) = rc*cos(phinew)
    jout.coil(*,1) = rc*sin(phinew)
;stop
;    tout.f = tout.f*mult(0)
;    tout.phi = tout.phi + dphi(0)/nn

    so = get_b_spectrum(g,jout,pnwant=pres+0.001,ntheta=ntheta,nphi=nphi)

    bro = so.br(im)
    pho = atan(so.br_s(im),so.br_c(im))

    mult = brv/bro
    dphi = phv-!dpi - pho

    err = sqrt((so.br_c(im)+sv.br_c(im))^2 + (so.br_s(im)+sv.br_s(im))^2)/bro
    iter += 1

    print,'Amp= ',brv,bro
    print,'Phase= ',phv,pho+!dpi
    print,'Error= ',err

endwhile

return,jout
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function screening_current,g,rmp,mm,nn,phidec=phidec,thmult=thmult

j=get_jparrat(g,mm,nn,phidec=phidec,thmult=thmult)

sv = get_b_spectrum(g,rmp,pnwant=j.pn)
sj = get_b_spectrum(g,j,pnwant=j.pn+0.001)
;sjp = get_b_spectrum(g,j,pnwant=j.pn+0.001)
;sjp2 = get_b_spectrum(g,j,pnwant=j.pn+0.01)
;sjm = get_b_spectrum(g,j,pnwant=j.pn-0.001)

im = where(sj.m eq mm)
brv = sv.br(im)
brj = sj.br(im)

mult = -brv/brj

jout=j
jout.current = jout.current*mult(0)

;stop

return,jout
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_ext_spectrum,g,rmp,pnstart,pnend,nsurfs
tstart = systime(1)
ntheta=511L
nphi=127L
r=dblarr(nsurfs,ntheta+1)
z=r & dsdr=r & dsdz=r
br = dblarr(nphi+1,ntheta+1) & bz=br & bphi=br & bnorm=br
spectrum = dblarr(nsurfs,nphi+1,ntheta+1)

theta = 2*!dpi*dindgen(ntheta+1)/ntheta
phi = 2*!dpi*dindgen(nphi+1)/nphi
a=findgen(nphi/2+2)
b=-reverse(findgen(nphi/2)+1)
nn = [a,b]
a=findgen(ntheta/2+2)
b=-reverse(findgen(ntheta/2)+1)
mm = [a,b]
;subn=sort(nn)
;nn = nn(subn)
;subm = sort(mm)
;mm = mm(subm)


if nsurfs gt 1 then $
  pnvals = (pnend-pnstart)*dindgen(nsurfs)/(nsurfs-1) + pnstart $
else pnvals = pnstart

q = interpol(g.qpsi,g.pn,pnvals)
normp = rmp
normp.current = 0.0d*normp.current
for i=0,nsurfs-1 do begin
    f=follow_fieldlines_phi(g,normp,pnvals(i),pnvals(i),1,2,0.,2*!dpi*q(i)/(100*ntheta),100*ntheta)
    r(i,*) = f.r(1,0:100*ntheta:100)
    z(i,*) = f.z(1,0:100*ntheta:100)
    p = get_psi_bicub(g,r(i,*),z(i,*),/derivs)
    normr = p.dsdr/(sqrt(p.dsdr^2 + p.dsdz^2))
    normz = p.dsdz/(sqrt(p.dsdr^2 + p.dsdz^2))
    for j=0,nphi do begin
        b = bfield_bs_cyl(r(i,*),phi(j),z(i,*),rmp.coil,rmp.current)
        br(j,*) = b.br
        bphi(j,*) = b.bphi
        bz(j,*) = b.bz
        bnorm(j,*) = b.br*normr + b.bz*normz
    endfor
    C=fft(bnorm,-1)
    spectrum(i,*,*) = C
;    for m=0,ntheta do begin
;        for n=0,nphi do begin
;            spectrum(i,subn(n),subm(m)) = C(n,m)
;        endfor
;    endfor
endfor

print,'elapsed time: ',systime(1)-tstart,' sec'
return,{pn:pnvals,q:q,mm:mm,nn:nn,spectrum:spectrum,bnorm:bnorm,r:r,z:z}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function follow_fieldlines_ripple,g,rmp,pnstart,pnend,nsurfs,ntheta,phistart,ds,nsteps, $
                                  rstart=rstart,zero=zero,zstart=zstart
tstart = systime(1)
x = dblarr(nsurfs*ntheta,nsteps+1)
z=x & y=x & s=x
px=x & py=x & pz=x & kg=x & magB=x

if nsurfs gt 1 then $
  pnvals = (pnend-pnstart)*dindgen(nsurfs)/(nsurfs-1) + pnstart $
else pnvals = pnstart
s(*,0) = 0


if ntheta gt 1 then theta = 2*!dpi*dindgen(ntheta)/ntheta else theta=0.0
nline = 1e4
for i=0,ntheta-1 do begin
    Lmax = max([g.r(g.mw-1)-g.r(0),g.z(g.mh-1) - g.z(0)])
    rline = g.rmaxis+Lmax*cos(theta(i))*dindgen(nline)/(nline-1)
    zline = g.zmaxis+Lmax*sin(theta(i))*dindgen(nline)/(nline-1)
    psiline = (-g.cpasma/abs(g.cpasma)*get_psi_bicub(g,rline,zline)-g.ssimag)/(g.ssibry-g.ssimag)
    ipn1 = where(psiline(0:nline-2) le 1.01d and psiline(1:nline-1) gt 1.01d,count)
    if count ge 1 then ipn1 = ipn1(0)
    iline = interpol(dindgen(ipn1+1),psiline(0:ipn1),pnvals)
    x(ntheta*dindgen(nsurfs)+i,0)=interpolate(rline,iline)*cos(phistart)
    y(ntheta*dindgen(nsurfs)+i,0)=interpolate(rline,iline)*sin(phistart)
    z(ntheta*dindgen(nsurfs)+i,0)=interpolate(zline,iline)
endfor

if keyword_set(rstart) then begin
    x(*,0) = pnvals*cos(phistart)
    y(*,0) = pnvals*sin(phistart)
    z(*,0) = 0.
endif

if keyword_set(zstart) then begin
    x(*,0) = zstart*cos(phistart)
    y(*,0) = zstart*sin(phistart)
    z(*,0) = pnvals
endif

if keyword_set(zstart) and keyword_set(rstart) then begin
    r = dblarr(n_elements(rstart),nsteps+1)
    z=r & phi=r
    x(*,0) = rstart*cos(phistart)
    y(*,0) = rstart*sin(phistart)
    z(*,0) = zstart
endif

psidat = get_psi_bicub(g,sqrt(x(*,0)^2+y(*,0)^2),z(*,0),/derivs)
dsdr=-g.cpasma/abs(g.cpasma)*psidat.dsdr/(g.ssibry-g.ssimag)
dsdz=-g.cpasma/abs(g.cpasma)*psidat.dsdz/(g.ssibry-g.ssimag)

px(*,0) = dsdr*x(*,0)/sqrt(x(*,0)^2+y(*,0)^2)
py(*,0) = dsdr*y(*,0)/sqrt(x(*,0)^2+y(*,0)^2)
pz(*,0) = dsdz

px_efit=px
py_efit=py
pz_efit=pz

magP = sqrt(px(*,0)^2 + py(*,0)^2 + pz(*,0)^2) ;
b=bfield_geq_cartes(g,x(*,0),y(*,0),z(*,0),/derivs,zero=zero)
if abs(max(rmp.current,/abs)) gt 1d-8 then begin
    b2=bfield_bs(x(*,0),y(*,0),z(*,0),rmp.coil,rmp.current,/derivs)
    b.bx += b2.bx
    b.by += b2.by
    b.bz += b2.bz
    b.dbxdx += b2.dbxdx
    b.dbxdy += b2.dbxdy
    b.dbxdz += b2.dbxdz
    b.dbydx += b2.dbydx
    b.dbydy += b2.dbydy
    b.dbydz += b2.dbydz
    b.dbzdx += b2.dbzdx
    b.dbzdy += b2.dbzdy
    b.dbzdz += b2.dbzdz
endif

bb=sqrt(b.bx^2 + b.by^2 + b.bz^2)
magB(*,0) = bb

hx = b.bx/bb                
hy = b.by/bb                
hz = b.bz/bb                

dBdx = (1/bb)*(b.bx*b.dbxdx+b.by*b.dbydx+b.bz*b.dbzdx) ;
dBdy = (1/bb)*(b.bx*b.dbxdy+b.by*b.dbydy+b.bz*b.dbzdy) ;
dBdz = (1/bb)*(b.bx*b.dbxdz+b.by*b.dbydz+b.bz*b.dbzdz) ;

dhxdx = (1/bb)*b.dbxdx - (b.bx/(bb^2))*dBdx ;
dhydx = (1/bb)*b.dbydx - (b.by/(bb^2))*dBdx ;
dhzdx = (1/bb)*b.dbzdx - (b.bz/(bb^2))*dBdx ;

dhxdy = (1/bb)*b.dbxdy - (b.bx/(bb^2))*dBdy ;
dhydy = (1/bb)*b.dbydy - (b.by/(bb^2))*dBdy ;
dhzdy = (1/bb)*b.dbzdy - (b.bz/(bb^2))*dBdy ;

dhxdz = (1/bb)*b.dbxdz - (b.bx/(bb^2))*dBdz ;
dhydz = (1/bb)*b.dbydz - (b.by/(bb^2))*dBdz ;
dhzdz = (1/bb)*b.dbzdz - (b.bz/(bb^2))*dBdz ;

hdotdelhx = hx*dhxdx + hy*dhxdy + hz*dhxdz ;
hdotdelhy = hx*dhydx + hy*dhydy + hz*dhydz ;
hdotdelhz = hx*dhzdx + hy*dhzdy + hz*dhzdz ;

hcrossx = hy*hdotdelhz - hz*hdotdelhy ;
hcrossy = hz*hdotdelhx - hx*hdotdelhz ;
hcrossz = hx*hdotdelhy - hy*hdotdelhx ;

kg(*,0) = hcrossx*(px(*,0)/magP) + hcrossy*(py(*,0)/magP) + hcrossz*(pz(*,0)/magP) ;
kg_efit = kg

if 1 then begin
    for i=1L,nsteps do begin
        s(*,i)=s(*,i-1)+ds

        b=bfield_geq_cartes(g,x(*,i-1),y(*,i-1),z(*,i-1),/derivs,zero=zero)
        if abs(max(rmp.current,/abs)) gt 1d-8 then begin
            b2=bfield_bs(x(*,i-1),y(*,i-1),z(*,i-1),rmp.coil,rmp.current,/derivs)
            b.bx += b2.bx
            b.by += b2.by
            b.bz += b2.bz
            b.dbxdx += b2.dbxdx
            b.dbxdy += b2.dbxdy
            b.dbxdz += b2.dbxdz
            b.dbydx += b2.dbydx
            b.dbydy += b2.dbydy
            b.dbydz += b2.dbydz
            b.dbzdx += b2.dbzdx
            b.dbzdy += b2.dbzdy
            b.dbzdz += b2.dbzdz
        endif
        btot = sqrt(b.bx^2+b.by^2+b.bz^2)
        k1x = ds*b.bx/btot
        k1y = ds*b.by/btot
        k1z = ds*b.bz/btot
        k1px = -ds*(b.dbxdx*px(*,i-1) + b.dbydx*py(*,i-1) + b.dbzdx*pz(*,i-1))/btot
        k1py = -ds*(b.dbxdy*px(*,i-1) + b.dbydy*py(*,i-1) + b.dbzdy*pz(*,i-1))/btot
        k1pz = -ds*(b.dbxdz*px(*,i-1) + b.dbydz*py(*,i-1) + b.dbzdz*pz(*,i-1))/btot

        b=bfield_geq_cartes(g,x(*,i-1)+0.5*k1x,y(*,i-1)+0.5*k1y,z(*,i-1)+0.5*k1z,/derivs,zero=zero)
        if abs(max(rmp.current,/abs)) gt 1d-8 then begin
            b2=bfield_bs(x(*,i-1)+0.5*k1x,y(*,i-1)+0.5*k1y,z(*,i-1)+0.5*k1z,rmp.coil,rmp.current,/derivs)
            b.bx += b2.bx
            b.by += b2.by
            b.bz += b2.bz
            b.dbxdx += b2.dbxdx
            b.dbxdy += b2.dbxdy
            b.dbxdz += b2.dbxdz
            b.dbydx += b2.dbydx
            b.dbydy += b2.dbydy
            b.dbydz += b2.dbydz
            b.dbzdx += b2.dbzdx
            b.dbzdy += b2.dbzdy
            b.dbzdz += b2.dbzdz
        endif
        btot = sqrt(b.bx^2+b.by^2+b.bz^2)
        k2x = ds*b.bx/btot
        k2y = ds*b.by/btot
        k2z = ds*b.bz/btot
        k2px = -ds*(b.dbxdx*(px(*,i-1)+0.5*k1px) + b.dbydx*(py(*,i-1)+0.5*k1py) + b.dbzdx*(pz(*,i-1)+0.5*k1pz))/btot
        k2py = -ds*(b.dbxdy*(px(*,i-1)+0.5*k1px) + b.dbydy*(py(*,i-1)+0.5*k1py) + b.dbzdy*(pz(*,i-1)+0.5*k1pz))/btot
        k2pz = -ds*(b.dbxdz*(px(*,i-1)+0.5*k1px) + b.dbydz*(py(*,i-1)+0.5*k1py) + b.dbzdz*(pz(*,i-1)+0.5*k1pz))/btot

        b=bfield_geq_cartes(g,x(*,i-1)+0.5*k2x,y(*,i-1)+0.5*k2y,z(*,i-1)+0.5*k2z,/derivs,zero=zero)
        if abs(max(rmp.current,/abs)) gt 1d-8 then begin
            b2=bfield_bs(x(*,i-1)+0.5*k2x,y(*,i-1)+0.5*k2y,z(*,i-1)+0.5*k2z,rmp.coil,rmp.current,/derivs)
            b.bx += b2.bx
            b.by += b2.by
            b.bz += b2.bz
            b.dbxdx += b2.dbxdx
            b.dbxdy += b2.dbxdy
            b.dbxdz += b2.dbxdz
            b.dbydx += b2.dbydx
            b.dbydy += b2.dbydy
            b.dbydz += b2.dbydz
            b.dbzdx += b2.dbzdx
            b.dbzdy += b2.dbzdy
            b.dbzdz += b2.dbzdz
        endif
        btot = sqrt(b.bx^2+b.by^2+b.bz^2)
        k3x = ds*b.bx/btot
        k3y = ds*b.by/btot
        k3z = ds*b.bz/btot
        k3px = -ds*(b.dbxdx*(px(*,i-1)+0.5*k2px) + b.dbydx*(py(*,i-1)+0.5*k2py) + b.dbzdx*(pz(*,i-1)+0.5*k2pz))/btot
        k3py = -ds*(b.dbxdy*(px(*,i-1)+0.5*k2px) + b.dbydy*(py(*,i-1)+0.5*k2py) + b.dbzdy*(pz(*,i-1)+0.5*k2pz))/btot
        k3pz = -ds*(b.dbxdz*(px(*,i-1)+0.5*k2px) + b.dbydz*(py(*,i-1)+0.5*k2py) + b.dbzdz*(pz(*,i-1)+0.5*k2pz))/btot

        b=bfield_geq_cartes(g,x(*,i-1)+k3x,y(*,i-1)+k3y,z(*,i-1)+k3z,/derivs,zero=zero)
        if abs(max(rmp.current,/abs)) gt 1d-8 then begin
            b2=bfield_bs(x(*,i-1)+k3x,y(*,i-1)+k3y,z(*,i-1)+k3z,rmp.coil,rmp.current,/derivs)
            b.bx += b2.bx
            b.by += b2.by
            b.bz += b2.bz
            b.dbxdx += b2.dbxdx
            b.dbxdy += b2.dbxdy
            b.dbxdz += b2.dbxdz
            b.dbydx += b2.dbydx
            b.dbydy += b2.dbydy
            b.dbydz += b2.dbydz
            b.dbzdx += b2.dbzdx
            b.dbzdy += b2.dbzdy
            b.dbzdz += b2.dbzdz
        endif
        btot = sqrt(b.bx^2+b.by^2+b.bz^2)
        k4x = ds*b.bx/btot
        k4y = ds*b.by/btot
        k4z = ds*b.bz/btot
        k4px = -ds*(b.dbxdx*(px(*,i-1)+k3px) + b.dbydx*(py(*,i-1)+k3py) + b.dbzdx*(pz(*,i-1)+k3pz))/btot
        k4py = -ds*(b.dbxdy*(px(*,i-1)+k3px) + b.dbydy*(py(*,i-1)+k3py) + b.dbzdy*(pz(*,i-1)+k3pz))/btot
        k4pz = -ds*(b.dbxdz*(px(*,i-1)+k3px) + b.dbydz*(py(*,i-1)+k3py) + b.dbzdz*(pz(*,i-1)+k3pz))/btot

        x(*,i) = x(*,i-1) + k1x/6 + k2x/3 + k3x/3 + k4x/6
        y(*,i) = y(*,i-1) + k1y/6 + k2y/3 + k3y/3 + k4y/6
        z(*,i) = z(*,i-1) + k1z/6 + k2z/3 + k3z/3 + k4z/6
        px(*,i) = px(*,i-1) + k1px/6 + k2px/3 + k3px/3 + k4px/6
        py(*,i) = py(*,i-1) + k1py/6 + k2py/3 + k3py/3 + k4py/6
        pz(*,i) = pz(*,i-1) + k1pz/6 + k2pz/3 + k3pz/3 + k4pz/6

        psidat = get_psi_bicub(g,sqrt(x(*,i)^2+y(*,i)^2),z(*,i),/derivs)
        dsdr=-g.cpasma/abs(g.cpasma)*psidat.dsdr/(g.ssibry-g.ssimag)
        dsdz=-g.cpasma/abs(g.cpasma)*psidat.dsdz/(g.ssibry-g.ssimag)
        px_efit(*,i) = dsdr*x(*,i)/sqrt(x(*,i)^2+y(*,i)^2)
        py_efit(*,i) = dsdr*y(*,i)/sqrt(x(*,i)^2+y(*,i)^2)
        pz_efit(*,i) = dsdz

        magP = sqrt(px(*,i)^2 + py(*,i)^2 + pz(*,i)^2) 
        magP_efit = sqrt(px_efit(*,i)^2 + py_efit(*,i)^2 + pz_efit(*,i)^2) 
        b=bfield_geq_cartes(g,x(*,i),y(*,i),z(*,i),/derivs,zero=zero)
        if abs(max(rmp.current,/abs)) gt 1d-8 then begin
            b2=bfield_bs(x(*,i),y(*,i),z(*,i)+k3z,rmp.coil,rmp.current,/derivs)
            b.bx += b2.bx
            b.by += b2.by
            b.bz += b2.bz
            b.dbxdx += b2.dbxdx
            b.dbxdy += b2.dbxdy
            b.dbxdz += b2.dbxdz
            b.dbydx += b2.dbydx
            b.dbydy += b2.dbydy
            b.dbydz += b2.dbydz
            b.dbzdx += b2.dbzdx
            b.dbzdy += b2.dbzdy
            b.dbzdz += b2.dbzdz
        endif

        bb=sqrt(b.bx^2 + b.by^2 + b.bz^2)
        magB(*,i) = bb

        hx = b.bx/bb                
        hy = b.by/bb                
        hz = b.bz/bb                

        dBdx = (1/bb)*(b.bx*b.dbxdx+b.by*b.dbydx+b.bz*b.dbzdx) ;
        dBdy = (1/bb)*(b.bx*b.dbxdy+b.by*b.dbydy+b.bz*b.dbzdy) ;
        dBdz = (1/bb)*(b.bx*b.dbxdz+b.by*b.dbydz+b.bz*b.dbzdz) ;

        dhxdx = (1/bb)*b.dbxdx - (b.bx/(bb^2))*dBdx ;
        dhydx = (1/bb)*b.dbydx - (b.by/(bb^2))*dBdx ;
        dhzdx = (1/bb)*b.dbzdx - (b.bz/(bb^2))*dBdx ;

        dhxdy = (1/bb)*b.dbxdy - (b.bx/(bb^2))*dBdy ;
        dhydy = (1/bb)*b.dbydy - (b.by/(bb^2))*dBdy ;
        dhzdy = (1/bb)*b.dbzdy - (b.bz/(bb^2))*dBdy ;
        
        dhxdz = (1/bb)*b.dbxdz - (b.bx/(bb^2))*dBdz ;
        dhydz = (1/bb)*b.dbydz - (b.by/(bb^2))*dBdz ;
        dhzdz = (1/bb)*b.dbzdz - (b.bz/(bb^2))*dBdz ;

        hdotdelhx = hx*dhxdx + hy*dhxdy + hz*dhxdz ;
        hdotdelhy = hx*dhydx + hy*dhydy + hz*dhydz ;
        hdotdelhz = hx*dhzdx + hy*dhzdy + hz*dhzdz ;

        hcrossx = hy*hdotdelhz - hz*hdotdelhy ;
        hcrossy = hz*hdotdelhx - hx*hdotdelhz ;
        hcrossz = hx*hdotdelhy - hy*hdotdelhx ;

        kg(*,i) = hcrossx*(px(*,i)/magP) + hcrossy*(py(*,i)/magP) + hcrossz*(pz(*,i)/magP) 
        kg_efit(*,i) = hcrossx*(px_efit(*,i)/magP_efit) + hcrossy*(py_efit(*,i)/magP_efit) + hcrossz*(pz_efit(*,i)/magP_efit) 

    endfor

endif
print,'elapsed time: ',systime(1)-tstart,' sec'
return,{x:x,y:y,z:z,s:s,magB:magB,magP:sqrt(px^2+py^2+pz^2),px:px,py:py,pz:pz,$
        magPe:sqrt(px_efit^2+py_efit^2+pz_efit^2),pxe:px_efit,pye:py_efit,pze:pz_efit,kg:kg,kge:kg_efit}
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_ripple,g,f,plotbp=plotbp,plotrip=plotrip

numbp = 500

s = size(f.x,/dim)
nlines = s(0)
nsteps = s(1)

R0 = g.rzero
B0 = abs(g.bcentr)

bhalf = 0.5d*(f.magB(*,0:nsteps-2)+f.magB(*,1:nsteps-1))
phalf = 0.5d*(f.magP(*,0:nsteps-2)+f.magP(*,1:nsteps-1))

ds = double(f.s(0,1)-f.s(0,0))

eps_eff_3_2 = dblarr(nlines)

nwin=0
z=mycolors()

for i=0,nlines-1 do begin
    I1 = total(ds/bhalf(i,*))
    I2 = total(ds*phalf(i,*)/bhalf(i,*))^(-2.0d)

    bmin = min(f.magB(i,*))
    bmax = max(f.magB(i,*))

    bp = bmin/B0 + (bmax/B0-bmin/B0)*dindgen(numbp)/(numbp-1)
    dbp = bp(1)-bp(0)
    
    sum_H2_I = dblarr(numbp)
    for ibp = 1,numbp-2 do begin
        ileft = where(f.magB(i,0:nsteps-2)/B0 ge bp(ibp) and f.magB(i,1:nsteps-1)/B0 lt bp(ibp),nleft)
        ileft = ileft+1
        iright = where(f.magB(i,0:nsteps-2)/B0 le bp(ibp) and f.magB(i,1:nsteps-1)/B0 gt bp(ibp),nright)
        if nleft eq 0 or nright eq 0 then continue
        if ileft(0) gt iright(0) then begin
            if nleft eq 1 and nright eq 1 then continue
            iright = iright(1:nright-1)
            nright = nright-1
        endif
        if ileft(nleft-1) gt iright(nright-1) then begin
            ileft = ileft(0:nleft-2)
            nleft = nleft-1
        endif
        if n_elements(ileft) ne n_elements(iright) then $
          print,'Bounce points are wrong: nleft,nright= ',nleft,nright
        if keyword_set(plotbp) then begin  
;            print,ileft,iright
            sleft = f.s(ileft)
            sright = f.s(iright)
            bleft = f.magB(i,ileft)
            bright = f.magB(i,iright)
            window,nwin
            nwin += 1
            plot,f.s,f.magB(i,*)
            oplot,sleft,bleft,psym=4,color=z.blue
            oplot,sright,bright,psym=6,color=z.red
        endif
        capH = dblarr(nleft) & capI=capH
        for irip = 0,nleft-1 do begin
            Bint = f.magB(i,ileft(irip):iright(irip))
            Pint = f.magP(i,ileft(irip):iright(irip))
            kgint = f.kg(i,ileft(irip):iright(irip))
            if keyword_set(plotrip) then begin
                window,nwin
                nwin += 1
                plot,f.s,f.magB(i,*)
                oplot,f.s(ileft(irip):iright(irip)),Bint,color=z.red,thick=2
            endif
            capH(irip) = (1/bp(ibp))*total((ds/Bint)*sqrt(bp(ibp)-Bint/B0)* $
                                           (4*(B0/Bint)-1/bp(ibp))*Pint*kgint)
            capI(irip) = total((ds/Bint)*sqrt(1-(Bint/(B0*bp(ibp)))))
        endfor
        sum_H2_I(ibp) = total(capH^2/capI)
    endfor
    I3 = total(dbp*sum_H2_I)
    eps_eff_3_2(i) = (!dpi*R0^2./(8*sqrt(2)))*(I1*I2*I3)
endfor

return,eps_eff_3_2^(2./3.)
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro plot_fieldlines_heatflux,g,s,filename,ps=ps,toffset=toffset,temp=temp,scale=scale,xrange=xrange,yrange=yrange,symsize=symsize

if keyword_set(ps) then begin
    dold = !d.name
    set_plot,'ps'
    !p.font=0
    device,/color,xsize=7,ysize=9,/inches,file='heatflux_fieldlines.ps',/helvetica,/bold
    charsize=2.0
    !p.charthick=5.0
    !p.thick=5.0
    !x.thick=5.0
    !y.thick=5.0
endif else charsize=2.

if not keyword_set(toffset) then toffset=0
if not keyword_set(scale) then scale=1.0
if not keyword_set(xrange) then xrange=[0.,1.5]
if not keyword_set(yrange) then yrange=[-1.8,-0.5]
if not keyword_set(symsize) then symsize = 0.001
tang =  360.+s.phi(0,s.ip0(0)+toffset)*180/!dpi
tang_nstx = 360 - (tang-90)
print,'Toroidal angle is: ',strtrim(tang,2)
print,'In NSTX system: ',strtrim(tang_nstx,2)

d=read_ascii(filename)
tmp=size(d.field1,/dimension)

rir = .01*d.field1(0,1:tmp(1)-1)
qir = d.field1(3,1:tmp(1)-1)
i60 = where(rir le 0.6)
rir = rir(i60)
qie = qir(i60)

z=mycolors()
inz = where(g.lim(0,*) ne 0.0 or g.lim(1,*) ne 0.)
rv = g.lim(0,inz)
zv = g.lim(1,inz)
ss=size(s.r(*,s.ip0+toffset),/dimension)
;print,ss
rplot = reform(s.r(*,s.ip0+toffset),ss(0)*ss(1))
zplot = reform(s.z(*,s.ip0+toffset),ss(0)*ss(1))
vesobj = obj_new('IDLanROI',rv,zv)
ins = vesobj->containspoints(rplot,zplot)
ins2 = where(ins gt 0.5)
;help,ins2
;help,rplot
rplot=rplot(ins2)
zplot=zplot(ins2)
;help,rplot
plot,rv,zv,/iso,charsize=charsize,xtitle='R (m)',ytitle='Z (m)',yrange=yrange,xrange=xrange
oplot,rv,zv,thick=4,color=z.blue
;oplot,s.r(*,s.ip0+toffset),s.z(*,s.ip0+toffset),psym=3
oplot,rplot,zplot,psym=1,symsize=symsize
if not keyword_set(temp)then oplot,rir,scale*qir*0.1 - 1.6,thick=8,color=z.red $
else oplot,rir,scale*qir*0.03 - 1.6,thick=6,color=z.red

if keyword_set(ps) then begin
    device,/close
    set_plot,dold
end
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro plot_footprint_heatflux,f,filename,ps=ps,scale=scale

if keyword_set(ps) then begin
    devold=!d.name
    set_plot,'ps'
    device,file='heatflux_footprint.ps',/inches,xsize=4,ysize=4,/color
    !p.thick=2.
    !x.thick=2.
    !y.thick=2.
    !p.charthick=2.
    !p.charsize=1.1
endif

d=read_ascii(filename)
tmp=size(d.field1,/dimension)
if not keyword_set(scale) then scale=10.
rir = .01*d.field1(0,1:tmp(1)-1)
qir = d.field1(1,1:tmp(1)-1)

z=mycolors()
plot,rir,qir,thick=0.1,xrange=[0.3,0.6],xtitle='R (m)',ytitle=textoidl('q_{IR} (MW/m^2)')
oplot,rir,qir,thick=2.0,color=z.brick
;oplot,f.rr,f.lconn*max(qir)/max(f.lconn),thick=2,color=z.blue
oplot,f.rr,scale*(f.psistart-f.psimin),color=z.blue,thick=2.0

if keyword_set(ps) then begin
    device,/close
    set_plot,devold
endif

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function find_manifolds,g,rmp,second=second,npts=npts,ntransits=ntransits,phioff=phioff

if not keyword_set(npts) then npts=300
if not keyword_set(ntransits) then ntransits=40
if not keyword_set(phioff) then phioff=0.0d

x=find_xpt(g,/refine,second=second)

eigp=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,!dpi/180,720,rstart=x.rx,zstart=x.zx)
eign=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,-!dpi/180,720,rstart=x.rx,zstart=x.zx)
rp = eigp.r(*,eigp.ip0)
zp = eigp.z(*,eigp.ip0)
rn = eign.r(*,eign.ip0)
zn = eign.z(*,eign.ip0)

rstartp = rp(0)+(dindgen(npts)/(npts-1))*(rp(1)-rp(0))
zstartp = zp(0)+(dindgen(npts)/(npts-1))*(zp(1)-zp(0))

rstartn = rn(0)+(dindgen(npts)/(npts-1))*(rn(1)-rn(0))
zstartn = zn(0)+(dindgen(npts)/(npts-1))*(zn(1)-zn(0))

pos=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,!dpi/180,ntransits*360,rstart=rstartp,zstart=zstartp)
neg=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,-!dpi/180,ntransits*360,rstart=rstartn,zstart=zstartn)

if keyword_set(second) then begin
    eigp=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,!dpi/180,720,rstart=x.rx2,zstart=x.zx2)
    eign=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,-!dpi/180,720,rstart=x.rx2,zstart=x.zx2)
    rp = eigp.r(*,eigp.ip0)
    zp = eigp.z(*,eigp.ip0)
    rn = eign.r(*,eign.ip0)
    zn = eign.z(*,eign.ip0)

    rstartp = rp(0)+(dindgen(npts)/(npts-1))*(rp(1)-rp(0))
    zstartp = zp(0)+(dindgen(npts)/(npts-1))*(zp(1)-zp(0))

    rstartn = rn(0)+(dindgen(npts)/(npts-1))*(rn(1)-rn(0))
    zstartn = zn(0)+(dindgen(npts)/(npts-1))*(zn(1)-zn(0))

    pos2=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,!dpi/180,ntransits*360,rstart=rstartp,zstart=zstartp)
    neg2=follow_fieldlines_phi(g,rmp,0.1,0.1,1,1,phioff,-!dpi/180,ntransits*360,rstart=rstartn,zstart=zstartn)
endif

if not keyword_set(second) then return,{pos:pos,neg:neg} $
  else return, {pos:pos,neg:neg,pos2:pos2,neg2:neg2}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function process_manifolds,man



return,0
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro plot_manifolds,man,overplot=overplot,second=second,indend=indend,allred=allred,psym=psym,toffset=toffset

z=mycolors()

if not keyword_set(toffset) then toffset=0
if not keyword_set(indend) then begin
    ipos = man.pos.ip0+toffset
    ineg = man.neg.ip0+toffset
    if keyword_set(second) then begin
        ipos2 = man.pos2.ip0
        ineg2 = man.neg2.ip0
    endif
endif else begin
    ipos = man.pos.ip0(0:indend(0))+toffset
    ineg = man.neg.ip0(0:indend(1))+toffset
    if keyword_set(second) then begin
        ipos2 = man.pos2.ip0(0:indend(2))
        ineg2 = man.neg2.ip0(0:indend(3))
    endif
endelse

if not keyword_set(allred) then $
  cols = [z.brick,z.blue,z.cyan,z.green] $
else cols = [z.brick,z.brick,z.brick,z.brick] 

if keyword_set(overplot) then oplot,man.neg.r(*,ineg),man.neg.z(*,ineg),thick=2,color=cols(0) $
else plot,man.neg.r(*,ineg),man.neg.z(*,ineg),thick=2,color=cols(0),psym=psym
oplot,man.pos.r(*,ipos),man.pos.z(*,ipos),thick=2,color=cols(1),psym=psym

if keyword_set(second) then begin
    oplot,man.pos2.r(*,ipos2),man.pos2.z(*,ipos2),thick=2,color=cols(2),psym=psym
    oplot,man.neg2.r(*,ineg2),man.neg2.z(*,ineg2),thick=2,color=cols(3),psym=psym
endif

end

pro compare_ipec_fields,g,rmp,ip,ps=ps,btype=btype

if not keyword_set(btype) then btype='vacuum'
if keyword_set(ps) then begin
    devold=!d.name
    set_plot,'ps'
    device,file='db_bs_ipec_'+btype+'.ps',/color
    !p.charthick=2.0
    !p.charsize=1.4
    !p.thick=2.0
    !x.thick=2.0
    !y.thick=2.0
endif

r=dindgen(100)/99 + 0.5
z=0.*r

phi = 0.*r
bc0 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi0 = bfield_ipec(ip,r,phi,z,btype=btype)

phi = 0.*r + 0.5*!dpi/3
bc1 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi1 = bfield_ipec(ip,r,phi,z,btype=btype)

phi = 0.*r + 1.0*!dpi/3
bc2 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi2 = bfield_ipec(ip,r,phi,z,btype=btype)

zz=mycolors()
if not keyword_set(ps) then window,0
plot,r,bc0.br,yrange=[-0.0025,0.001],xtitle='R (m)',ytitle=textoidl('B_R (T)')
oplot,r,bi0.br,linestyle=2
oplot,r,bc1.br,color=zz.brick
oplot,r,bi1.br,linestyle=2,color=zz.brick
oplot,r,bc2.br,color=zz.blue
oplot,r,bi2.br,linestyle=2,color=zz.blue
legend,[textoidl('\phi=0'),textoidl('\phi=30'),textoidl('\phi=60')],colors=[zz.black,zz.brick,zz.blue]

if not keyword_set(ps) then window,1
plot,r,bc0.bz,yrange=[-0.0001,0.0001],xtitle='R (m)',ytitle=textoidl('B_Z (T)')
oplot,r,bi0.bz,linestyle=2
oplot,r,bc1.bz,color=zz.brick
oplot,r,bi1.bz,linestyle=2,color=zz.brick
oplot,r,bc2.bz,color=zz.blue
oplot,r,bi2.bz,linestyle=2,color=zz.blue
legend,[textoidl('\phi=0'),textoidl('\phi=30'),textoidl('\phi=60')],colors=[zz.black,zz.brick,zz.blue]

if not keyword_set(ps) then window,2
plot,r,bc0.bphi,yrange=[-0.003,0.002],xtitle='R (m)',ytitle=textoidl('B_\phi (T)')
oplot,r,bi0.bphi,linestyle=2
oplot,r,bc1.bphi,color=zz.brick
oplot,r,bi1.bphi,linestyle=2,color=zz.brick
oplot,r,bc2.bphi,color=zz.blue
oplot,r,bi2.bphi,linestyle=2,color=zz.blue
legend,[textoidl('\phi=0'),textoidl('\phi=30'),textoidl('\phi=60')],colors=[zz.black,zz.brick,zz.blue]

z=3*dindgen(100)/99 -1.5
r=0.*z + 1.4

phi = 0.*r
bc0 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi0 = bfield_ipec(ip,r,phi,z,btype=btype)

phi = 0.*r + 0.5*!dpi/3
bc1 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi1 = bfield_ipec(ip,r,phi,z,btype=btype)

phi = 0.*r + 1.0*!dpi/3
bc2 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi2 = bfield_ipec(ip,r,phi,z,btype=btype)

if not keyword_set(ps) then window,3
plot,z,bc0.br,yrange=[-0.002,0.0005],xtitle='Z (m)',ytitle=textoidl('B_R (T)')
oplot,z,bi0.br,linestyle=2
oplot,z,bc1.br,color=zz.brick
oplot,z,bi1.br,linestyle=2,color=zz.brick
oplot,z,bc2.br,color=zz.blue
oplot,z,bi2.br,linestyle=2,color=zz.blue
legend,[textoidl('\phi=0'),textoidl('\phi=30'),textoidl('\phi=60')],colors=[zz.black,zz.brick,zz.blue]

if not keyword_set(ps) then window,4
plot,z,bc0.bz,yrange=[-0.001,0.001],xtitle='Z (m)',ytitle=textoidl('B_Z (T)')
oplot,z,bi0.bz,linestyle=2
oplot,z,bc1.bz,color=zz.brick
oplot,z,bi1.bz,linestyle=2,color=zz.brick
oplot,z,bc2.bz,color=zz.blue
oplot,z,bi2.bz,linestyle=2,color=zz.blue
legend,[textoidl('\phi=0'),textoidl('\phi=30'),textoidl('\phi=60')],colors=[zz.black,zz.brick,zz.blue]

if not keyword_set(ps) then window,5
plot,z,bc0.bphi,yrange=[-0.0025,0.0025],xtitle='Z (m)',ytitle=textoidl('B_\phi (T)')
oplot,z,bi0.bphi,linestyle=2
oplot,z,bc1.bphi,color=zz.brick
oplot,z,bi1.bphi,linestyle=2,color=zz.brick
oplot,z,bc2.bphi,color=zz.blue
oplot,z,bi2.bphi,linestyle=2,color=zz.blue
legend,[textoidl('\phi=0'),textoidl('\phi=30'),textoidl('\phi=60')],colors=[zz.black,zz.brick,zz.blue]

phi=2*!dpi*dindgen(100)/99

r=0.*phi + 1.4
z=0.*phi
bc0 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi0 = bfield_ipec(ip,r,phi,z,btype=btype)

r=0.*phi + 1.2
z=0.*phi
bc1 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi1 = bfield_ipec(ip,r,phi,z,btype=btype)

r=0.*phi + 1.2
z=0.*phi + 0.5
bc2 = bfield_bs_cyl(r,phi,z,rmp.coil,rmp.current)
bi2 = bfield_ipec(ip,r,phi,z,btype=btype)

if not keyword_set(ps) then window,6
plot,phi,bc0.br,yrange=[-0.002,0.002],xtitle='Phi',ytitle=textoidl('B_R (T)')
oplot,phi,bi0.br,linestyle=2
oplot,phi,bc1.br,color=zz.brick
oplot,phi,bi1.br,linestyle=2,color=zz.brick
oplot,phi,bc2.br,color=zz.blue
oplot,phi,bi2.br,linestyle=2,color=zz.blue
legend,['R=1.4, Z=0.0','R=1.2, Z=0.0','R=1.2, Z=0.5'],colors=[zz.black,zz.brick,zz.blue]

if not keyword_set(ps) then window,7
plot,phi,bc0.bz,yrange=[-0.0006,0.0006],xtitle='Phi',ytitle=textoidl('B_Z (T)')
oplot,phi,bi0.bz,linestyle=2
oplot,phi,bc1.bz,color=zz.brick
oplot,phi,bi1.bz,linestyle=2,color=zz.brick
oplot,phi,bc2.bz,color=zz.blue
oplot,phi,bi2.bz,linestyle=2,color=zz.blue
legend,['R=1.4, Z=0.0','R=1.2, Z=0.0','R=1.2, Z=0.5'],colors=[zz.black,zz.brick,zz.blue]

if not keyword_set(ps) then window,8
plot,phi,bc0.bphi,yrange=[-0.0025,0.0025],xtitle='Phi',ytitle=textoidl('B_\phi (T)')
oplot,phi,bi0.bphi,linestyle=2
oplot,phi,bc1.bphi,color=zz.brick
oplot,phi,bi1.bphi,linestyle=2,color=zz.brick
oplot,phi,bc2.bphi,color=zz.blue
oplot,phi,bi2.bphi,linestyle=2,color=zz.blue
legend,['R=1.4, Z=0.0','R=1.2, Z=0.0','R=1.2, Z=0.5'],colors=[zz.black,zz.brick,zz.blue]

if keyword_set(ps) then begin
    device,/close
    set_plot,devold
endif


end





function get_bfield_singleslice_old,s,r,z,iphi

inds = lonarr(n_elements(r),2)
sr = reform(s.r(*,*,iphi))
sz = reform(s.z(*,*,iphi))
for i=0,n_elements(r)-1 do begin
    dist = sqrt((s.r(*,*,iphi(i))-r(i))^2+(s.z(*,*,iphi(i))-z(i))^2)
;    dist = sqrt((sr-r(i))^2+(sz-z(i))^2)
;    help,dist
    junk = min(dist,imin)
    inds(i,*) = array_indices(dist,imin)
endfor

mp = inds(*,1)-1
iw = where(mp lt 0,count)
if count gt 0 then mp(iw) = n_elements(s.r(0,*,0))-2
;if mp lt 0 then mp = n_elements(s.r(0,*,0))-2
pp = inds(*,1)+1
iw = where(pp gt n_elements(s.r(0,*,0))-1,count)
if count gt 0 then pp(iw) = 1
;if pp gt n_elements(s.r(0,*,0))-1 then pp = 1
;pr = inds(0)+1
;if pr gt n_elements(s.r


r00 = s.r(inds(*,0),inds(*,1),iphi)
r10 = s.r(inds(*,0)+1,inds(*,1),iphi)
r01 = s.r(inds(*,0),pp,iphi)
r11 = s.r(inds(*,0)+1,pp,iphi)
rm0 = s.r(inds(*,0)-1,inds(*,1),iphi)
rm1 = s.r(inds(*,0)-1,pp,iphi)
rmm = s.r(inds(*,0)-1,mp,iphi)
r1m = s.r(inds(*,0)+1,mp,iphi)
r0m = s.r(inds(*,0),mp,iphi)

z00 = s.z(inds(*,0),inds(*,1),iphi)
z10 = s.z(inds(*,0)+1,inds(*,1),iphi)
z01 = s.z(inds(*,0),pp,iphi)
z11 = s.z(inds(*,0)+1,pp,iphi)
zm0 = s.z(inds(*,0)-1,inds(*,1),iphi)
zm1 = s.z(inds(*,0)-1,pp,iphi)
zmm = s.z(inds(*,0)-1,mp,iphi)
z1m = s.z(inds(*,0)+1,mp,iphi)
z0m = s.z(inds(*,0),mp,iphi)

br00 = s.br(inds(*,0),inds(*,1),iphi)
br10 = s.br(inds(*,0)+1,inds(*,1),iphi)
br01 = s.br(inds(*,0),pp,iphi)
br11 = s.br(inds(*,0)+1,pp,iphi)
brm0 = s.br(inds(*,0)-1,inds(*,1),iphi)
brm1 = s.br(inds(*,0)-1,pp,iphi)
brmm = s.br(inds(*,0)-1,mp,iphi)
br1m = s.br(inds(*,0)+1,mp,iphi)
br0m = s.br(inds(*,0),mp,iphi)

bz00 = s.bz(inds(*,0),inds(*,1),iphi)
bz10 = s.bz(inds(*,0)+1,inds(*,1),iphi)
bz01 = s.bz(inds(*,0),pp,iphi)
bz11 = s.bz(inds(*,0)+1,pp,iphi)
bzm0 = s.bz(inds(*,0)-1,inds(*,1),iphi)
bzm1 = s.bz(inds(*,0)-1,pp,iphi)
bzmm = s.bz(inds(*,0)-1,mp,iphi)
bz1m = s.bz(inds(*,0)+1,mp,iphi)
bz0m = s.bz(inds(*,0),mp,iphi)

bphi00 = s.bphi(inds(*,0),inds(*,1),iphi)
bphi10 = s.bphi(inds(*,0)+1,inds(*,1),iphi)
bphi01 = s.bphi(inds(*,0),pp,iphi)
bphi11 = s.bphi(inds(*,0)+1,pp,iphi)
bphim0 = s.bphi(inds(*,0)-1,inds(*,1),iphi)
bphim1 = s.bphi(inds(*,0)-1,pp,iphi)
bphimm = s.bphi(inds(*,0)-1,mp,iphi)
bphi1m = s.bphi(inds(*,0)+1,mp,iphi)
bphi0m = s.bphi(inds(*,0),mp,iphi)

;print,[r00,z00],[r01,z01],[r11,z11],[r,z]
;print,[r00,z00],[rm1,zm1],[r01,z01],[r,z]
;print,[r00,z00],[r11,z11],[r10,z10],[r,z]
;print,[r00,z00],[r10,z10],[r0m,z0m],[r,z]
;print,[r00,z00],[r0m,z0m],[rmm,zmm],[r,z]
;print,[r00,z00],[rmm,zmm],[rm0,zm0],[r,z]

;print,in_tri([r00,z00],[r01,z01],[r11,z11],[r,z])
;print,in_tri([r00,z00],[rm1,zm1],[r01,z01],[r,z])
;print,in_tri([r00,z00],[r11,z11],[r10,z10],[r,z])
;print,in_tri([r00,z00],[r10,z10],[r0m,z0m],[r,z])
;print,in_tri([r00,z00],[r0m,z0m],[rmm,zmm],[r,z])
;print,in_tri([r00,z00],[rmm,zmm],[rm0,zm0],[r,z])

bar1=in_tri([[r00],[z00]],[[r01],[z01]],[[r11],[z11]],[[r],[z]])
bar2=in_tri([[r00],[z00]],[[r11],[z11]],[[r10],[z10]],[[r],[z]])
bar3=in_tri([[r00],[z00]],[[r10],[z10]],[[r0m],[z0m]],[[r],[z]])
bar4=in_tri([[r00],[z00]],[[r0m],[z0m]],[[rmm],[zmm]],[[r],[z]])
bar5=in_tri([[r00],[z00]],[[rmm],[zmm]],[[rm0],[zm0]],[[r],[z]])
bar6=in_tri([[r00],[z00]],[[rm0],[zm0]],[[r01],[z01]],[[r],[z]])
bar7=in_tri([[r0m],[z0m]],[[r10],[z10]],[[r1m],[z1m]],[[r],[z]])
bar8=in_tri([[rm0],[zm0]],[[rm1],[zm1]],[[r01],[z01]],[[r],[z]])

ins = bar1.ins + bar2.ins +bar3.ins +bar4.ins +bar5.ins +bar6.ins + bar7.ins+ bar8.ins
;ins = [[bar1.ins],[bar2.ins],[bar3.ins],[bar4.ins],[bar5.ins],[bar6.ins],[bar7.ins],[bar8.ins]]
;stop

iprob = where(ins ne 1,count)
if count gt 0 then begin
    print,'Problem with finding triangle surrounding points'
    print,'r= ',r(iprob)
    print,'z= ',z(iprob)
    print,'ins= ',ins(iprob)
    print,'u = ',bar1.u(iprob),bar2.u(iprob),bar3.u(iprob),bar4.u(iprob),bar5.u(iprob),bar6.u(iprob),bar7.u(iprob),bar8.u(iprob)
    print,'v = ',bar1.v(iprob),bar2.v(iprob),bar3.v(iprob),bar4.v(iprob),bar5.v(iprob),bar6.v(iprob),bar7.v(iprob),bar8.v(iprob)
;    print,'u+v = ',bar1.u+bar1.v,bar2.u+bar2.v,bar3.u+bar3.v,bar4.u+bar4.v,bar5.u+bar5.v,bar6.u+bar6.v
;    print,'tri1: ',[r00,z00],[r01,z01],[r11,z11],[r,z]
;    print,'tri2: ',[r00,z00],[rm1,zm1],[r01,z01],[r,z]
;    print,'tri3: ',[r00,z00],[r11,z11],[r10,z10],[r,z]
;    print,'tri4: ',[r00,z00],[r10,z10],[r0m,z0m],[r,z]
;    print,'tri5: ',[r00,z00],[r0m,z0m],[rmm,zmm],[r,z]
;    print,'tri6: ',[r00,z00],[rmm,zmm],[rm0,zm0],[r,z]
    col=mycolors()
    plot,[r00,r01,r11,r00,r11,r10,r00,r10,r0m,r00,r0m,rmm,r00,rmm,rm0,r00,rm0,r01],[z00,z01,z11,z00,z11,z10,z00,z10,z0m,z00,z0m,zmm,z00,zmm,zm0,z00,zm0,z01],xrange=[min([rm1,rmm,rm0]),max([r1m,r10,r11])],yrange=[min([zmm,z0m,z1m]),max([zm1,z01,z11])],/iso
;    plot,[r00,r01,r11,r00,rm1,r01,r00,r11,r10,r00,r10,r0m,r00,r0m,rmm,r00,rmm,rm0],[z00,z01,z11,z00,zm1,z01,z00,z11,z10,z00,z10,z0m,z00,z0m,zmm,z00,zmm,zm0]
    oplot,[r0m,r10,r1m,r0m],[z0m,z10,z1m,z0m],color=col.blue
    oplot,[rm0,rm1,r01,rm0],[zm0,zm1,z01,zm0],color=col.red
    oplot,[r,r],[z,z],psym=1,color=col.brick
    stop
endif

br1 = br00 + bar1.u*(br01-br00) + bar1.v(br11-br00)
br2 = br00 + bar2.u*(brm1-br00) + bar2.v(br01-br00)
br3 = br00 + bar3.u*(br11-br00) + bar3.v(br10-br00)
br4 = br00 + bar4.u*(br10-br00) + bar4.v(br0m-br00)
br5 = br00 + bar5.u*(br0m-br00) + bar5.v(brmm-br00)
br6 = br00 + bar6.u*(brmm-br00) + bar6.v(brm0-br00)
br7 = br00 + bar7.u*(brmm-br00) + bar7.v(brm0-br00)
br8 = br00 + bar8.u*(brmm-br00) + bar8.v(brm0-br00)

bz1 = bz00 + bar1.u*(bz01-bz00) + bar1.v(bz11-bz00)
bz2 = bz00 + bar2.u*(bzm1-bz00) + bar2.v(bz01-bz00)
bz3 = bz00 + bar3.u*(bz11-bz00) + bar3.v(bz10-bz00)
bz4 = bz00 + bar4.u*(bz10-bz00) + bar4.v(bz0m-bz00)
bz5 = bz00 + bar5.u*(bz0m-bz00) + bar5.v(bzmm-bz00)
bz6 = bz00 + bar6.u*(bzmm-bz00) + bar6.v(bzm0-bz00)
bz7 = bz00 + bar7.u*(bzmm-bz00) + bar7.v(bzm0-bz00)
bz8 = bz00 + bar8.u*(bzmm-bz00) + bar8.v(bzm0-bz00)

bphi1 = bphi00 + bar1.u*(bphi01-bphi00) + bar1.v(bphi11-bphi00)
bphi2 = bphi00 + bar2.u*(bphim1-bphi00) + bar2.v(bphi01-bphi00)
bphi3 = bphi00 + bar3.u*(bphi11-bphi00) + bar3.v(bphi10-bphi00)
bphi4 = bphi00 + bar4.u*(bphi10-bphi00) + bar4.v(bphi0m-bphi00)
bphi5 = bphi00 + bar5.u*(bphi0m-bphi00) + bar5.v(bphimm-bphi00)
bphi6 = bphi00 + bar6.u*(bphimm-bphi00) + bar6.v(bphim0-bphi00)
bphi7 = bphi00 + bar7.u*(bphimm-bphi00) + bar7.v(bphim0-bphi00)
bphi8 = bphi00 + bar8.u*(bphimm-bphi00) + bar8.v(bphim0-bphi00)

br = br1*bar1.ins + br2*bar2.ins +br3*bar3.ins +br4*bar4.ins +br5*bar5.ins +br6*bar6.ins+br7*bar7.ins+br8*bar8.ins
bz = bz1*bar1.ins + bz2*bar2.ins +bz3*bar3.ins +bz4*bar4.ins +bz5*bar5.ins +bz6*bar6.ins+bz7*bar7.ins+bz8*bar8.ins
bphi = bphi1*bar1.ins + bphi2*bar2.ins +bphi3*bar3.ins +bphi4*bar4.ins +bphi5*bar5.ins +bphi6*bar6.ins+bphi7*bar7.ins+bphi8*bar8.ins

return,{br:br,bphi:bphi,bz:bz}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_rflm_bicub_coeffs,s

nr = n_elements(s(*,0))
nz = n_elements(s(0,*))

dsdr = (shift(s,-1,0)-shift(s,1,0))/2
dsdz = (shift(s,0,-1)-shift(s,0,1))/2
d2sdrdz = (shift(s,-1,-1) - shift(s,1,-1) - shift(s,-1,1) + shift(s,1,1))/4

a = get_bicub_mat()

c=dblarr(nr*nz,4,4)
for ir=0,nr-2 do begin
    for iz=0,nz-2 do begin
        index = iz + nz*ir
        b = [ s(ir,iz),        s(ir+1,iz),        s(ir,iz+1),        s(ir+1,iz+1), $
              dsdr(ir,iz),      dsdr(ir+1,iz),      dsdr(ir,iz+1),      dsdr(ir+1,iz+1), $
              dsdz(ir,iz),      dsdz(ir+1,iz),      dsdz(ir,iz+1),      dsdz(ir+1,iz+1), $
              d2sdrdz(ir,iz),d2sdrdz(ir+1,iz),d2sdrdz(ir,iz+1),d2sdrdz(ir+1,iz+1)]
        coeff = transpose(invert(a))##b
        c(index,*,0) = coeff(0:3)
        c(index,*,1) = coeff(4:7)
        c(index,*,2) = coeff(8:11)
        c(index,*,3) = coeff(12:15)
    endfor
endfor

return,c
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_rflm_bicub,cin,is,iu,nu,derivs=derivs

ir = floor(is)
iz = floor(iu)
index = iz + nu*ir

dir = is-ir
diz = iu-iz
c = cin(index,*,*)
;stop
dir2 = dir^2
dir3 = dir^3
diz2 = diz^2
diz3 = diz^3

psi1 = c(*,0,0)       + c(*,1,0)*dir       + c(*,2,0)*dir2       + c(*,3,0)*dir3 + $
  c(*,0,1)*diz   + c(*,1,1)*dir*diz   + c(*,2,1)*dir2*diz   + c(*,3,1)*dir3*diz + $
  c(*,0,2)*diz2 + c(*,1,2)*dir*diz2 + c(*,2,2)*dir2*diz2 + c(*,3,2)*dir3*diz2 + $
  c(*,0,3)*diz3 + c(*,1,3)*dir*diz3 + c(*,2,3)*dir2*diz3 + c(*,3,3)*dir3*diz3

return,psi1

end

function rflm,a,ns_start=ns_start,nu_start=nu_start,ntransits=ntransits,s_start=s_start,u_start=u_start
if not keyword_set(ns_start) then ns_start=1
if not keyword_set(nu_start) then nu_start=1
if not keyword_set(ntransits) then ntransits=2
if keyword_set(s_start) then ns_start=n_elements(s_start)
if keyword_set(u_start) then nu_start=n_elements(u_start)
sp = dblarr(ns_start,nu_start,ntransits)
up = sp & rp=sp & zp=sp

ntot = n_elements(a.s)

s1tmp = a.s(0:ntot/2-1)
u1tmp = a.theta(0:ntot/2-1)
r1tmp = a.r(0:ntot/2-1)
z1tmp = a.z(0:ntot/2-1)

s2tmp = a.s(ntot/2:ntot-1)
u2tmp = a.theta(ntot/2:ntot-1)
r2tmp = a.r(ntot/2:ntot-1)
z2tmp = a.z(ntot/2:ntot-1)

is1 = where(s1tmp-s1tmp(0) ge 1e-3)
nu = is1(0)
ns = ntot/nu/2
nunew = nu+1
s1=dblarr(ns,nunew)
u1=s1 & r1=s1 & z1=S1 & s2=s1 & u2=s1 & r2=s1 & z2=s1

s1tmp = transpose(reform(s1tmp,[nu,ns]))
u1tmp = transpose(reform(u1tmp,[nu,ns]))
r1tmp = transpose(reform(r1tmp,[nu,ns]))
z1tmp = transpose(reform(z1tmp,[nu,ns]))

s2tmp = transpose(reform(s2tmp,[nu,ns]))
u2tmp = transpose(reform(u2tmp,[nu,ns]))
r2tmp = transpose(reform(r2tmp,[nu,ns]))
z2tmp = transpose(reform(z2tmp,[nu,ns]))

s1(*,0:nu-1) = s1tmp
s1(*,nunew-1) = s1tmp(*,0)
u1(*,0:nu-1) = u1tmp
u1(*,nunew-1) = u1tmp(*,0) + 2*!dpi
r1(*,0:nu-1) = r1tmp
r1(*,nunew-1) = r1tmp(*,0)
z1(*,0:nu-1) = z1tmp
z1(*,nunew-1) = z1tmp(*,0)

s2(*,0:nu-1) = s2tmp
s2(*,nunew-1) = s2tmp(*,0)
u2(*,0:nu-1) = u2tmp 
u2(*,nunew-1) = u2tmp(*,0)
r2(*,0:nu-1) = r2tmp
r2(*,nunew-1) = r2tmp(*,0)
z2(*,0:nu-1) = z2tmp
z2(*,nunew-1) = z2tmp(*,0)

ig = where(u2 lt u1,count)
;if count gt 0 then u2(ig) += 2*!dpi

;s2coeff = get_rflm_bicub_coeffs(s2)
;u2coeff = get_rflm_bicub_coeffs(u2)
;;print,get_rflm_bicub(s2coeff,12.5,13,nunew),s2(12,13)

;stop
smin = min(s1)
smax = max(s1)
umin = min(u1)
umax = max(u1)

if not keyword_set(s_start) then begin
    s_start = (smax-smin)*dindgen(ns_start+1)/(ns_start+1) + smin
    s_start = s_start(1:ns_start)
endif 

if not keyword_set(u_start) then begin
    u_start = (max(u1)-min(u1))*dindgen(nu_start+1)/(nu_start+1) + min(u1)
    u_start = u_start(1:nu_start)
endif

for i=0,ns_start-1 do begin
    sp(i,*,0) = s_start(i)
    up(i,*,0) = u_start
endfor

du = u1(0,1)-u1(0,0)
;ds = s1(1,0)-s1(0,0)
ds = (s1(ns-1,0)-s1(0,0))/(ns-1)

is0 = (sp(*,*,0)-smin)/ds
iu0 = (up(*,*,0)-umin)/du
rp(*,*,0) = interpolate(r1,is0,iu0)
zp(*,*,0) = interpolate(z1,is0,iu0)

; plot,u2,s2,psym=3;,xrange=[2.5,3.5],yrange=[0.4,0.6]
; for i = 0,ns-2 do begin
;     for j = 0,nu-2 do begin
;         oplot,[u2(i,j),u2(i+1,j),u2(i+1,j+1),u2(i,j+1),u2(i,j)],[s2(i,j),s2(i+1,j),s2(i+1,j+1),s2(i,j+1),s2(i,j)]
;     endfor
; endfor

;stop
for i=1L,ntransits-1 do begin
;   is = (sp(*,*,i-1)-smin)/ds
;   iu = (up(*,*,i-1)-umin)/du
    is = interpol(dindgen(ns),s1(*,0),sp(*,*,i-1))
    iu = interpol(dindgen(nunew),reform(u1(0,*)),up(*,*,i-1))
;    print,is,is2,iu,iu2
    fis = floor(is)
    fiu = floor(iu)
;   ts = (sp(*,*,i-1) - s1(fis,fiu))/ds
;   tu = (up(*,*,i-1) - u1(fis,fiu))/du
    ts = is - floor(is)
    tu = iu - floor(iu)
;    print,ts,ts2,tu,tu2
;    print,'Error 1 s: ',interpolate(s1,is,iu)-sp(*,*,i-1)
;    print,'Error 2 s: ',s1(fis,fiu)*(1-ts)*(1-tu) + s1(fis+1,fiu)*ts*(1-tu) + s1(fis,fiu+1)*(1-ts)*tu + s1(fis+1,fiu+1)*ts*tu - sp(*,*,i-1)
;   print,'Error 3 s: ',s1(fis,fiu)*(1-ts) + s1(fis+1,fiu)*ts - sp(*,*,i-1)
;   print,'Error 3 u: ',u1(fis,fiu)*(1-tu) + u1(fis,fiu+1)*tu - up(*,*,i-1)
    stmp = s1(fis,fiu)*(1-ts) + s1(fis+1,fiu)*ts
    tserr = (stmp-sp(*,*,i-1))/ds
    ts = ts-tserr
    utmp = u1(fis,fiu)*(1-tu) + u1(fis,fiu+1)*tu
    tuerr = (utmp-up(*,*,i-1))/du
    tu = tu-tuerr
;   print,'Error corrected s: ',s1(fis,fiu)*(1-ts) + s1(fis+1,fiu)*ts - sp(*,*,i-1)
;   print,'Error corrected u: ',u1(fis,fiu)*(1-tu) + u1(fis,fiu+1)*tu - up(*,*,i-1)
;   sp(*,*,i) = s2(fis,fiu)*(1-ts)*(1-tu) + s2(fis+1,fiu)*ts*(1-tu) + s2(fis,fiu+1)*(1-ts)*tu + s2(fis+1,fiu+1)*ts*tu
;   up(*,*,i) = (u2(fis,fiu)*(1-ts)*(1-tu) + u2(fis+1,fiu)*ts*(1-tu) + u2(fis,fiu+1)*(1-ts)*tu + u2(fis+1,fiu+1)*ts*tu) mod (2*!dpi)
;   rp(*,*,i) = r2(fis,fiu)*(1-ts)*(1-tu) + r2(fis+1,fiu)*ts*(1-tu) + r2(fis,fiu+1)*(1-ts)*tu + r2(fis+1,fiu+1)*ts*tu
;   zp(*,*,i) = z2(fis,fiu)*(1-ts)*(1-tu) + z2(fis+1,fiu)*ts*(1-tu) + z2(fis,fiu+1)*(1-ts)*tu + z2(fis+1,fiu+1)*ts*tu
    sp(*,*,i) = interpolate(s2,is,iu,cubic=-0.5)
    up(*,*,i) = interpolate(u2,is,iu,cubic=-0.5) mod (2*!dpi)
;   tmp = get_rflm_bicub(s2coeff,is,iu,nunew)
;   sp(*,*,i) = reform(tmp,[ns_start,nu_start])
;   tmp = get_rflm_bicub(u2coeff,is,iu,nunew)
;   up(*,*,i) = reform(tmp,[ns_start,nu_start]) mod (2*!dpi)
;    stop
    rp(*,*,i) = interpolate(r2,is,iu)
    zp(*,*,i) = interpolate(z2,is,iu)
;    print,'Next s ',reform(sp(*,*,i))
;    print,'Next u ',reform(up(*,*,i))
endfor
;stop
;help
return,{s:sp,u:up,r:rp,z:zp}
end

function convert_poincare,s,ns,nu,nt

sg = (reform(s.s,[nu*ns,nt]))
ug = (reform(s.theta,[nu*ns,nt]))
rg = (reform(s.r,[nu*ns,nt]))
zg = (reform(s.z,[nu*ns,nt]))
Lg = 2*!dpi*rg
Lg = total(Lg,2,/cumulative)
ip0 = dindgen(nt)
;stop

return,{r:rg,z:zg,theta:ug,psin:sg,L:Lg,ip0:ip0}
end

function convert_rflm,s

ns = n_elements(s.s(*,0,0))
nu = n_elements(s.s(0,*,0))
nt = n_elements(s.s(0,0,*))
sg = dblarr(nu*ns,nt)
ug = sg & rg=sg& zg=sg

for i=0,ns-1 do begin
    sg(i*nu:(i+1)*nu-1,*) = s.s(i,*,*)
    ug(i*nu:(i+1)*nu-1,*) = s.u(i,*,*)
    rg(i*nu:(i+1)*nu-1,*) = s.r(i,*,*)
    zg(i*nu:(i+1)*nu-1,*) = s.z(i,*,*)
endfor

;sg = (reform(s.s,[nu*ns,nt]))
;ug = (reform(s.u,[nu*ns,nt]))
;rg = (reform(s.r,[nu*ns,nt]))
;zg = (reform(s.z,[nu*ns,nt]))
Lg = 2*!dpi*rg
Lg = total(Lg,2,/cumulative)
ip0 = dindgen(nt)
;stop

return,{r:rg,z:zg,theta:ug,psin:sg,L:Lg,ip0:ip0}
end



pro test_rflm,ps=ps

;restore,'siesta73pert_p05top95_36surfs_10pol_150transits.sav'
;restore,'poincare_127317_7_3perturb_p05p95_36surfs_10pol_200transits.sav'
;restore,'poincare_135183_M31N1_p05p95_36surfs_10pol_200transits.sav'
restore,'poincare_135183_M31N1_p05p95_36surfs_10pol_200transits.sav'
;ntot = n_elements(s.s)
;is1 = where(s.s-s.s(0) ge 1e-3)
nu = 10
nt = 201
ns = 36
sg = (reform(s.s,[nu,ns,nt]))
ug = (reform(s.theta,[nu,ns,nt]))

m44=read_poincare('nstx_135183-bfield_tracing-M31N1_rflm_s50u100.dat')
m14=read_poincare('nstx_135183-bfield_tracing-M31N1_rflm_s50u100.dat')
m51=read_poincare('nstx_135183-bfield_tracing-M31N1_rflm_s50u100.dat')

;m44=read_poincare('nstx-bfield_tracing-M31N1_largepert_rflm_s40u40.dat')
;m14=read_poincare('nstx-bfield_tracing-M31N1_largepert_rflm_s100u40.dat')
;m51=read_poincare('nstx-bfield_tracing-M31N1_largepert_rflm_s50u100.dat')

;m44=read_poincare('nstx-bfield_tracing-127317_7_3perturb_rflm_s40u40.dat')
;m14=read_poincare('nstx-bfield_tracing-127317_7_3perturb_rflm_s100u40.dat')
;m41=read_poincare('nstx-bfield_tracing-127317_7_3perturb_rflm_s40u100.dat')
;m11=read_poincare('nstx-bfield_tracing-127317_7_3perturb_rflm_s100u100.dat')

b44=rflm(m44,ntransits=nt,s_start=reform(sg(0,*,0)),u_start=reform(ug(*,0,0)))
b14=rflm(m14,ntransits=nt,s_start=reform(sg(0,*,0)),u_start=reform(ug(*,0,0)))
b51=rflm(m51,ntransits=nt,s_start=reform(sg(0,*,0)),u_start=reform(ug(*,0,0)))
;b41=rflm(m41,ntransits=nt,s_start=reform(sg(0,*,0)),u_start=reform(ug(*,0,0)))
;b11=rflm(m11,ntransits=nt,s_start=reform(sg(0,*,0)),u_start=reform(ug(*,0,0)))
;b44=rflm(m44,ntransits=200,ns=35,nu=10)
;b41=rflm(m41,ntransits=200,ns=35,nu=10)
;b11=rflm(m11,ntransits=200,ns=35,nu=10)

if keyword_set(ps) then begin
    set_plot,'ps'
    !p.font=0
    device,file='check_rflm.ps',/color,/inches,xsize=7.5,ysize=10,xoffset=0.5,yoffset=0.5,/helvetica,/bold
endif

pos = plot_positions(nrow=3,ncol=2,vertspac=0.1)
z=mycolors()
cf = color_flipper()
make_nice_plots,ps=ps
!p.position = pos(*,0)
nsurf=ns
nsurf=ns
nlines=nu
plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),title=textoidl('Poincare')
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[sg(0,i,0),sg(0,i,0)],color=cf(i)
    oplot,ug(*,i,*),sg(*,i,*),psym=3,color=cf(i)
;    oplot,[0,0.1*2*!dpi],[s.psiN(i*nlines,0),s.psiN(i*nlines,0)],color=cf(i)
;    oplot,s.theta(i*nlines:(i+1)*nlines-1,s.ip0),s.psiN(i*nlines:(i+1)*nlines-1,s.ip0),psym=3,color=cf(i)
endfor

!p.position = pos(*,1)
nsurf=36
nlines=10
plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),/noerase,title=textoidl('Poincare')
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[sg(0,i,0),sg(0,i,0)],color=cf(i)
    oplot,ug(*,i,*),sg(*,i,*),psym=3,color=cf(i)
;    oplot,[0,0.1*2*!dpi],[s.psiN(i*nlines,0),s.psiN(i*nlines,0)],color=cf(i)
;    oplot,s.theta(i*nlines:(i+1)*nlines-1,s.ip0),s.psiN(i*nlines:(i+1)*nlines-1,s.ip0),psym=3,color=cf(i)
endfor

!p.position = pos(*,2)
nsurf=36
nlines=10
plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),/noerase,title=textoidl('Poincare')
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[sg(0,i,0),sg(0,i,0)],color=cf(i)
    oplot,ug(*,i,*),sg(*,i,*),psym=3,color=cf(i)
;    oplot,[0,0.1*2*!dpi],[s.psiN(i*nlines,0),s.psiN(i*nlines,0)],color=cf(i)
;    oplot,s.theta(i*nlines:(i+1)*nlines-1,s.ip0),s.psiN(i*nlines:(i+1)*nlines-1,s.ip0),psym=3,color=cf(i)
endfor

!p.position = pos(*,3)
nsurf=n_elements(b44.s(*,0,0))
nlines=n_elements(b44.s(0,*,0))

plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),title=textoidl('RFLM: s40,u40'),/noerase
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[b44.s(i,0,0),b44.s(i,0,0)],color=cf(i)
    oplot,b44.u(i,*,*),b44.s(i,*,*),psym=3,color=cf(i)
endfor

!p.position = pos(*,4)
nsurf=n_elements(b14.s(*,0,0))
nlines=n_elements(b14.s(0,*,0))

plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),title=textoidl('RFLM: s100,u40'),/noerase
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[b14.s(i,0,0),b14.s(i,0,0)],color=cf(i)
    oplot,b14.u(i,*,*),b14.s(i,*,*),psym=3,color=cf(i)
endfor

!p.position = pos(*,5)
nsurf=n_elements(b51.s(*,0,0))
nlines=n_elements(b51.s(0,*,0))

plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),title=textoidl('RFLM: s50,u100'),/noerase
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[b51.s(i,0,0),b51.s(i,0,0)],color=cf(i)
    oplot,b51.u(i,*,*),b51.s(i,*,*),psym=3,color=cf(i)
endfor

if keyword_set(ps) then begin
    device,/close
    set_plot,'x'
endif

!p.position=0

end


pro compare_rflm,file1,file2,ps=ps

m1 = read_poincare(file1)
m2 = read_poincare(file2)

b1=rflm(m1,ntransits=200,ns=60,nu=10)
b2=rflm(m2,ntransits=200,ns=60,nu=10)

if keyword_set(ps) then begin
    set_plot,'ps'
    !p.font=0
    device,file='compare_rflm.ps',/color,/inches,xsize=7.5,ysize=3,xoffset=0.5,yoffset=0.5,/helvetica,/bold
endif

pos = plot_positions(nrow=1,ncol=2)
z=mycolors()
cf = color_flipper()
make_nice_plots,ps=ps
!p.position = pos(*,0)
nsurf=n_elements(b1.s(*,0,0))
nlines=n_elements(b1.s(0,*,0))

plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),title=file1,charsize=0.5,xcharsize=2.0,ycharsize=2.0
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[b1.s(i,0,0),b1.s(i,0,0)],color=cf(i)
    oplot,b1.u(i,*,*),b1.s(i,*,*),psym=3,color=cf(i)
endfor

!p.position = pos(*,1)
nsurf=n_elements(b2.s(*,0,0))
nlines=n_elements(b2.s(0,*,0))

plot,[0,0],[0,0],xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,ytitle=textoidl('\psi_{tor}^{1/2}'),xtitle=textoidl('\theta'),title=file2,/noerase,charsize=0.5,xcharsize=2.0,ycharsize=2.0
for i = 0,nsurf-1 do begin
    oplot,[0,0.1*2*!dpi],[b2.s(i,0,0),b2.s(i,0,0)],color=cf(i)
    oplot,b2.u(i,*,*),b2.s(i,*,*),psym=3,color=cf(i)
endfor

if keyword_set(ps) then begin
    device,/close
    set_plot,'x'
endif

!p.position=0

end


function spiral_2d,g,rmp,n=n,plot=plot
if not keyword_set(n) then n=3

ps = 1.001
pe = 1.075
nlines=500

rhit = dblarr(nlines)
phit = rhit
zhit=rhit

il = where(g.lim(0,*) gt 1e-4)
rlim = g.lim(0,il)
zlim = g.lim(1,il)
vesobj = obj_new('IDLanROI',rlim,zlim)

phistart=0.0
s=follow_fieldlines_phi(g,rmp,ps,pe,nlines,1,phistart,-!dpi/180.,3600)
for i=0,nlines-1 do begin
    ins = vesobj->containspoints(s.r(i,*),s.z(i,*))
    iw = where(ins eq 0,count)
    if count ne 0 then begin
        rhit(i)=interpol(s.r(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
        phit(i)=interpol(s.phi(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
        zhit(i)=interpol(s.z(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
    endif else print,'Line ',i,'didn''t hit wall'
endfor

phit = phit mod (2*!dpi) + 2*!dpi

rn = dblarr(nlines,n)
phin = rn
zn = rn

for i=0,n-1 do begin
    rn(*,i) = rhit
    zn(*,i) = zhit
    phin(*,i) = (phit + i*2*!dpi/n) mod (2*!dpi)
endfor

if keyword_set(plot) then begin
    plot,phin,rn,xrange=[0,2*!dpi],yrange=[0.3,0.6],xstyle=1,ystyle=1,psym=1
endif

return,{r:rn,phi:phin,z:zn}
end


function transport_spiral_torweight,g,rmp,n=n,plot=plot,eps=eps,nlines=nlines,yrange=yrange,ps=ps
if not keyword_set(n) then n=3
if not keyword_set(eps) then eps=0.2

ps = 1.001
pe = 1.075
if not keyword_set(nlines) then nlines=500
ntor = 360

rhit = dblarr(nlines)
phit = rhit
zhit=rhit

il = where(g.lim(0,*) gt 1e-4)
rlim = g.lim(0,il)
zlim = g.lim(1,il)
vesobj = obj_new('IDLanROI',rlim,zlim)

phistart=0.0
s=follow_fieldlines_phi(g,rmp,ps,pe,nlines,1,phistart,-!dpi/180.,3600)
for i=0,nlines-1 do begin
    ins = vesobj->containspoints(s.r(i,*),s.z(i,*))
    iw = where(ins eq 0,count)
    if count ne 0 then begin
        rhit(i)=interpol(s.r(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
        phit(i)=interpol(s.phi(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
        zhit(i)=interpol(s.z(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
    endif else print,'Line ',i,'didn''t hit wall'
endfor

phit = phit mod (2*!dpi) + 2*!dpi

rn = dblarr(nlines,ntor)
phin = rn
zn = rn
wn = rn

for i=0,ntor-1 do begin
    phitmp = i*2*!dpi/ntor
    fac = 1 + eps*cos(n*phitmp)
    wn(*,i) = fac
    rn(*,i) = rhit
    zn(*,i) = zhit
    phin(*,i) = (phit + phitmp) mod (2*!dpi)
endfor



if keyword_set(plot) then begin
   if keyword_set(ps) then begin
      set_plot,'ps'
      !p.font=0
      device,file='transport_spirals.ps',/inches,xsize=4,ysize=3.5,xoff=.5,yoff=.5,/helvetica,/bold,/color,bits_per_pixel=8
   endif
   z=mycolors()
   make_nice_plots,ps=ps
   levs = (max(wn)-min(wn))*dindgen(50)/49 + min(wn)
    contour,wn,phin,rn,/irregular,levels=levs,/fill,xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,xtitle=textoidl('\phi'),ytitle='R (m)'
   if keyword_set(ps) then begin
      device,/close
      set_plot,'x'
   endif
endif

return,{r:rn,phi:phin,z:zn,w:wn}
end

function transport_spiral_torpos,g,rmp,n=n,plot=plot,eps=eps,nlines=nlines,yrange=yrange,ps=ps
if not keyword_set(n) then n=3
if not keyword_set(eps) then eps=0.2

ps = 1.001
pe = 1.075
if not keyword_set(nlines) then nlines = 50
ntor = 360
lam_nat = 0.01*(0.1/0.02)
dlam = 0.003*(0.1/0.02)

rhit = dblarr(nlines)
phit = rhit
zhit=rhit

il = where(g.lim(0,*) gt 1e-4)
rlim = g.lim(0,il)
zlim = g.lim(1,il)
vesobj = obj_new('IDLanROI',rlim,zlim)

phistart=0.0
s=follow_fieldlines_phi(g,rmp,ps,pe,nlines,1,phistart,-!dpi/180.,3600)
pnstart = s.psin(*,0)
rstart = s.r(*,0)
zstart = s.z(*,0)

for i=0,nlines-1 do begin
    ins = vesobj->containspoints(s.r(i,*),s.z(i,*))
    iw = where(ins eq 0,count)
    if count ne 0 then begin
        rhit(i)=interpol(s.r(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
        phit(i)=interpol(s.phi(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
        zhit(i)=interpol(s.z(i,0:iw(0)),s.z(i,0:iw(0)),-1.603) 
    endif else print,'Line ',i,'didn''t hit wall'
endfor

phit = phit mod (2*!dpi) + 2*!dpi

rn = dblarr(nlines,ntor)
phin = rn
zn = rn
wn = rn
x0 = dblarr(ntor)

for i=0,ntor-1 do begin
    phitmp = i*2*!dpi/ntor
    x0(i) = 1.0 + dlam + dlam*cos(n*phitmp)
    junk = min(pnstart - x0(i),/abs,im)
;    wn(0:im,i) = 0.0d
;    wn(im:nlines-1,i) = exp(-(pnstart(im:nlines-1)-x0(i))/lam_nat)
    wn(*,i) = exp(-(pnstart-x0(i))/lam_nat)
;    stop
    rn(*,i) = rhit
    zn(*,i) = zhit
    phin(*,i) = (phit + phitmp) mod (2*!dpi)
endfor

if keyword_set(plot) then begin
   if keyword_set(ps) then begin
      set_plot,'ps'
      !p.font=0
      device,file='transport_spirals_pos.ps',/inches,xsize=4,ysize=3.5,xoff=.5,yoff=.5,/helvetica,/bold
   endif
   levs = (max(wn)-min(wn))*dindgen(50)/49 + min(wn)
   contour,wn,phin,rn,/irregular,levels=levs,/fill,xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,xtitle=textoidl('\phi'),ytitle='R (m)'
   
   if keyword_set(ps) then begin
      device,/close
      set_plot,'x'
   endif
endif

return,{r:rn,phi:phin,z:zn,w:wn,rstart:rstart,zstart:zstart,psistart:pnstart,x0:x0}
end



function transport_spiral_phiplane,g,rmp,n=n,plot=plot,eps=eps,nlines=nlines,xrange=xrange,yrange=yrange,phiwant=phiwant,iso=iso,lambda=lambda, $
                                   psiplot=psiplot,radplot=radplot,makeps=makeps,nturn=nturn,nphi=nphi,npol=npol,pstart=pstart,pend=pend
if not keyword_set(n) then n=3
if not keyword_set(eps) then eps=0.2
if not keyword_set(phiwant) then phiwant = 0.0d
if eps lt 1e-5 then eps = 0.0d

if keyword_set(lambda) then begin
   rline = (max(g.bdry(0,*))+0.1-g.rmaxis)*dindgen(1000)/999 + g.rmaxis
   ptmp = get_psi_bicub(g,rline,0*rline+g.zmaxis)
   pnline = (-ptmp-g.ssimag)/(g.ssibry-g.ssimag)
   rsep = interpol(rline,pnline,1.0d)
   lambda_psi = interpol(pnline,rline,rsep+lambda) - 1.0d
endif

if not keyword_set(pstart) then pstart = 1.0001
if not keyword_set(pend) then pend = 1.5
if not keyword_set(nlines) then nlines=500
ntor = 360

il = where(g.lim(0,*) gt 1e-4)
rlim = g.lim(0,il)
zlim = g.lim(1,il)
vesobj = obj_new('IDLanROI',rlim,zlim)

phistart=0.0
if not keyword_set(nphi) then nphi = 3600
s=follow_fieldlines_phi(g,rmp,pstart,pend,nlines,1,phistart,-!dpi/180.,nphi)
s2=follow_fieldlines_phi(g,rmp,pstart,pend,nlines,1,phistart,!dpi/180.,nphi)

isol = where(s.psin(*,0) ge 1.0d)
icore = where(s.psin(*,0) lt 1.0d,ncore)

if keyword_set(npol) and ncore gt 0 then begin
   for i=0,max(icore) do begin
;      if s.psin(i,0) lt 1.0d then begin
         if s.theta(i,1) lt !dpi then ipi = where((s.theta(i,0:nphi-1) lt !dpi and s.theta(i,1:nphi) ge !dpi) or $
                          (abs(s.theta(i,1:nphi)-s.theta(i,0:nphi-1)) gt !dpi and s.theta(i,0:nphi-1) gt !dpi),count) $
         else ipi = where((s.theta(i,0:nphi-1) gt !dpi and s.theta(i,1:nphi) le !dpi) or $
                          (abs(s.theta(i,1:nphi)-s.theta(i,0:nphi-1)) gt !dpi and s.theta(i,0:nphi-1) gt 0.0d),count)
         if count ge npol then begin
            s.r(i,ipi(npol-1)+1:nphi) *= 0.0d
            s.z(i,ipi(npol-1)+1:nphi) *= 0.0d
;            s.theta(i,ipi(npol-1):nphi) *= 0.0d
         endif else print,'Line ',strtrim(i,2),' only goes ',strtrim(count,2),' poloidal transits'
 ;     endif
 ;     if s2.psin(i,0) lt 1.0d then begin
         if s2.theta(i,1) lt !dpi then ipi = where((s2.theta(i,0:nphi-1) lt !dpi and s2.theta(i,1:nphi) ge !dpi) or $
                          (abs(s2.theta(i,1:nphi)-s2.theta(i,0:nphi-1)) gt !dpi and s2.theta(i,0:nphi-1) gt !dpi),count) $
         else ipi = where((s2.theta(i,0:nphi-1) gt !dpi and s2.theta(i,1:nphi) le !dpi) or $
                          (abs(s2.theta(i,1:nphi)-s2.theta(i,0:nphi-1)) gt !dpi and s2.theta(i,0:nphi-1) gt 0.0d),count)
         if count ge npol then begin
            s2.r(i,ipi(npol-1)+1:nphi) *= 0.0d
            s2.z(i,ipi(npol-1)+1:nphi) *= 0.0d
;            s2.theta(i,ipi(npol-1):nphi) *= 0.0d
         endif else print,'Line ',strtrim(i,2),' only goes ',strtrim(count,2),' poloidal transits'
;      endif
   endfor
endif


nip0 = n_elements(s.ip0)
for j=0,ntor-1 do begin
;   for i=0,nlines-1 do begin
      phitmp = j*2*!dpi/ntor
      fac = 1 + eps*cos(n*phitmp)
      junk = min(phitmp-abs(s.phi(0,s.ip0(0):s.ip0(1))),ioff,/abs)
      if n_elements(rip0) eq 0 then rip0 = reform(s.r(*,s.ip0+ioff),nlines*nip0) else rip0 = [rip0, reform(s.r(*,s.ip0+ioff),nlines*nip0)]
      if n_elements(zip0) eq 0 then zip0 = reform(s.z(*,s.ip0+ioff),nlines*nip0) else zip0 = [zip0, reform(s.z(*,s.ip0+ioff),nlines*nip0)]
      if n_elements(pnip0) eq 0 then pnip0 = reform(s.psin(*,s.ip0+ioff),nlines*nip0) else pnip0 = [pnip0, reform(s.psin(*,s.ip0+ioff),nlines*nip0)]
;      if n_elements(wip0) eq 0 then wip0 = 0.0*dblarr(nlines*nip0) + fac else wip0 = [wip0,0.0*dblarr(nlines*nip0) + fac
      if keyword_set(lambda) then radfac = exp(-(reform(s.psin(*,s.ip0+ioff),nlines*nip0)-1.0d)/lambda_psi) else radfac = 1.0 + 0.0*dblarr(nlines*nip0)
      ig1 = where(radfac gt 1.0,count)
      if count gt 0 then radfac(ig1) = 1.0
      if n_elements(wip0) eq 0 then wip0 = radfac*fac else wip0 = [wip0,radfac*fac]
      if n_elements(radip0) eq 0 then radip0 = radfac else radip0 = [radip0,radfac]
;      stop
endfor

nip0 = n_elements(s2.ip0)
for j=0,ntor-1 do begin
;   for i=0,nlines-1 do begin
      phitmp = j*2*!dpi/ntor
      fac = 1 + eps*cos(n*phitmp)
      junk = min(phitmp-abs(s2.phi(0,s2.ip0(0):s2.ip0(1))),ioff,/abs)
      if keyword_set(lambda) then radfac = exp(-(reform(s2.psin(*,s2.ip0+ioff),nlines*nip0)-1.0d)/lambda_psi) else radfac = 1.0 + 0.0*dblarr(nlines*nip0)
       ig1 = where(radfac gt 1.0,count)
      if count gt 0 then radfac(ig1) = 1.0
      if n_elements(rip0) eq 0 then rip0 = reform(s2.r(*,s2.ip0+ioff),nlines*nip0) else rip0 = [rip0, reform(s2.r(*,s2.ip0+ioff),nlines*nip0)]
      if n_elements(zip0) eq 0 then zip0 = reform(s2.z(*,s2.ip0+ioff),nlines*nip0) else zip0 = [zip0, reform(s2.z(*,s2.ip0+ioff),nlines*nip0)]
      if n_elements(pnip0) eq 0 then pnip0 = reform(s2.psin(*,s2.ip0+ioff),nlines*nip0) else pnip0 = [pnip0, reform(s2.psin(*,s2.ip0+ioff),nlines*nip0)]
;      if n_elements(wip0) eq 0 then wip0 = 0.0*dblarr(nlines*nip0) + fac else wip0 = [wip0,0.0*dblarr(nlines*nip0) + fac
      if n_elements(wip0) eq 0 then wip0 = radfac*fac else wip0 = [wip0,radfac*fac]
      if n_elements(radip0) eq 0 then radip0 = radfac else radip0 = [radip0,radfac]
endfor

;ins = vesobj->containspoints(rip0,zip0)
;iw = where(ins gt 0)
;rip0 = rip0(iw)
;zip0 = zip0(iw)
;wip0 = wip0(iw)
;pnip0 = pnip0(iw)
;radip0 = radip0(iw)

;radip0 = exp(-(pnip0-1.0d)/lambda)

;;rr = g.r(0) + (g.r(g.mw-1)-g.r(0))*dindgen(nr)/(nr-1)
;;zz = g.z(0) + (g.z(g.mw-1)-g.z(0))*dindgen(nz)/(nz-1)
;rr = min(rip0) + (max(rip0)+1e-3-min(rip0))*dindgen(nr)/(nr-1)
;zz = min(zip0) + (max(zip0)+1e-3-min(zip0))*dindgen(nz)/(nz-1)
;nhit = dblarr(nr,nz)
;for ir = 0,nr-2 do begin
;   for iz=0,nz-2 do begin
;      junk = where(rip0 ge rr(ir) and rip0 lt rr(ir+1) and zip0 ge zz(iz) and zip0 lt zz(iz+1),count)
;      nhit(ir,iz) = count
;   endfor
;endfor



if keyword_set(plot) then begin
   if keyword_set(makeps) then begin
      set_plot,'ps'
      !p.font=0
      device,file='transport_lobes.ps',/color,/inches,xsize=5,ysize=8,xoff=.75,yoff=.75,/helvetica,/bold
   endif 

;   ncolor = !d.table_size
;   maxc = max(wip0)
;   minc = min(wip0)
;   ib = where(g.bdry(0,*) gt 1e-4)
;   rb = g.bdry(0,il)
;   zb = g.bdry(1,il)
;   z=mycolors()
;   plot,rlim,zlim,xrange=xrange,yrange=yrange,xtitle='R (m)',ytitle='Z (m)',iso=iso,color=z.black,thick=8
;   for i=0L,n_elements(wip0)-1 do begin
;      oplot,[rip0(i),rip0(i)],[zip0(i),zip0(i)],psym=3,color=(ncolor-2)*(wip0(i)-minc)/(maxc-minc)
;   endfor
;;   oplot,rb,zb,color=z.white,thick=2.0
;   
   ;r1,z1,w1 are inputs, grid_input reduces number of points to make sure triangulate works
   grid_input,rip0,zip0,wip0,r2,z2,w2,duplicates='Avg',epsilon = 5e-4	;r2,z2,w2 are outputs
   triangulate,r2,z2,triangles                                  ;triangles is an output

;setup R & Z for grid output:
   rgrid = findgen(280)*0.005 + 1.0
   zgrid = findgen(540)*0.005 - 1.4

   grid = trigrid(r2,z2,w2,triangles,xout=rgrid,yout=zgrid) ;still messy inside separatrix

;zero out inside separatrix…
;   ob=obj_new('IDLanROI',s.r(0,*),s.z(0,*))
   icross = where(abs(s.theta(0,1:nphi)-s.theta(0,0:nphi-1)) gt !dpi,ncross)
   if ncross gt 0 then begin
      if s.theta(0,1) lt 0.0 then izero = icross(0) 
      if s.theta(0,1) gt 0.0 and ncross gt 1 then izero = icross(1)
   endif else izero = nphi
;   if izero(0) eq 0 then izero = izero(1) else izero = izero(0)
   if s.psin(0,0) lt 1.0 then ob=obj_new('IDLanROI',s.r(0,0:izero),s.z(0,0:izero)) $
   else ob=obj_new('IDLanROI',g.bdry)
;stop
   nrgrid = n_elements(rgrid)
   nzgrid = n_elements(zgrid)
   rgrid2d = reform(rgrid # replicate(1,nzgrid),nrgrid*nzgrid)
   zgrid2d = reform(zgrid ## replicate(1,nrgrid),nrgrid*nzgrid)
   zeropts = ob->containspoints(rgrid2d,zgrid2d)
   obj_destroy,ob
   grid[where(zeropts eq 1)]*=0
   
   sepobj = obj_new('IDLanROI',g.bdry)
   psin2d = dblarr(n_elements(rgrid),n_elements(zgrid))
   for iz = 0,n_elements(zgrid)-1 do begin
      ptmp = get_psi_bicub(g,rgrid,0.0*rgrid+zgrid(iz))
      psin2d(*,iz) = (-ptmp-g.ssimag)/(g.ssibry-g.ssimag)
   endfor
   ipfr = where((psin2d lt 1.0d) and (sepobj->containspoints(rgrid2d,zgrid2d) eq 0))
   grid(ipfr)*=0
   ives = where(vesobj->containspoints(rgrid2d,zgrid2d) eq 0)
   grid(ives)*=0
   grid(where(psin2d gt max(s.psin)))*=0

;   levels=(max(grid)-min(grid))*dindgen(50)/49 + min(grid)
   if keyword_set(lambda) then levels=(1.0+eps-0.1)*dindgen(50)/49+0.1 else levels = 2.0d*eps*dindgen(50)/49 + 1.0d - eps

   make_nice_plots,ps=makeps
   z=mycolors()

   plot,rlim,zlim,xrange=xrange,yrange=yrange,xtitle='R (m)',ytitle='Z (m)',iso=iso,color=z.black,thick=8
   contour,grid,rgrid,zgrid,/fill,levels=levels,xrange=xrange,yrange=yrange,/overplot
   ib = where(g.bdry(0,*) gt 1e-4)
   rb = g.bdry(0,il)
   zb = g.bdry(1,il)
   oplot,rb,zb,color=z.white
   
   pn = dindgen(g.mw)/(g.mw-1)
   pn123=interpol(pn,g.qpsi,12.0/3.0)
   pn113=interpol(pn,g.qpsi,11.0/3.0)
   pn103=interpol(pn,g.qpsi,10.0/3.0)

;   contour,psin2d,rgrid,zgrid,levels=[pn103,pn113,pn123],color=z.black,/overplot

   if keyword_set(makeps) then begin
      device,/close
      set_plot,'x'
   endif

;    contour,wn,phin,rn,/irregular,nlevels=50,/fill,xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,xtitle=textoidl('\phi'),ytitle='R (m)'
endif


if keyword_set(psiplot) then begin
   ncolor = !d.table_size
   maxc = max(pnip0)
   minc = min(pnip0)
   ib = where(g.bdry(0,*) gt 1e-4)
   rb = g.bdry(0,il)
   zb = g.bdry(1,il)
   z=mycolors()
   plot,rb,zb,xrange=xrange,yrange=yrange,xtitle='R (m)',ytitle='Z (m)',iso=iso
   for i=0L,n_elements(pnip0)-1 do begin
      oplot,[rip0(i),rip0(i)],[zip0(i),zip0(i)],psym=3,color=(ncolor-2)*(pnip0(i)-minc)/(maxc-minc)
   endfor

   oplot,rb,zb,color=z.white,thick=2.0

;    contour,wn,phin,rn,/irregular,nlevels=50,/fill,xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,xtitle=textoidl('\phi'),ytitle='R (m)'
endif

if keyword_set(radplot) then begin
   ncolor = !d.table_size
   maxc = max(radip0)
   minc = min(radip0)
   ib = where(g.bdry(0,*) gt 1e-4)
   rb = g.bdry(0,il)
   zb = g.bdry(1,il)
   z=mycolors()
   plot,rb,zb,xrange=xrange,yrange=yrange,xtitle='R (m)',ytitle='Z (m)',iso=iso
   for i=0L,n_elements(radip0)-1 do begin
      oplot,[rip0(i),rip0(i)],[zip0(i),zip0(i)],psym=3,color=(ncolor-2)*(radip0(i)-minc)/(maxc-minc)
   endfor
   oplot,rb,zb,color=z.white,thick=2.0

;    contour,wn,phin,rn,/irregular,nlevels=50,/fill,xrange=[0,2*!dpi],yrange=yrange,xstyle=1,ystyle=1,xtitle=textoidl('\phi'),ytitle='R (m)'
endif

return,{r:rip0,z:zip0,w:wip0}
end


function transport_spiral_diffuse,g,rmp,pnstart=pnstart,nsurfs=nsurfs,ntheta=ntheta,ntor=ntor,diffuse=diffuse,nsteps=nsteps,$
                                  nbins=nbins,sigtheta=sigtheta,sigphi=sigphi

if not keyword_set(nsurfs) then nsurfs=10L
if not keyword_set(ntheta) then ntheta=10L
if not keyword_set(pnstart) then pnstart = 0.98
if not keyword_set(diffuse) then diffuse=1d-4
if not keyword_set(nsteps) then nsteps=3600L
if not keyword_set(ntor) then ntor = 10L

phistart = 2*!dpi*dindgen(ntor)/ntor

s=follow_fieldlines_phi(g,rmp,pnstart,pnstart,nsurfs,ntheta,0.0d,-!dpi/180,nsteps,diffuse=diffuse,sigtheta=sigtheta,sigphi=sigphi)
;s = replicate(s,ntor)

;for j = 1,ntor-1 do begin
;   s(j)=follow_fieldlines_phi(g,rmp,pnstart,pnstart,nsurfs,ntheta,phistart(j),-!dpi/180,nsteps,diffuse=diffuse,sigtheta=sigtheta,sigphi=sigphi)
;endfor

il = where(g.lim(0,*) gt 1e-4)
rlim = g.lim(0,il)
zlim = g.lim(1,il)
vesobj = obj_new('IDLanROI',rlim,zlim)

for i=0,nsurfs*ntheta-1 do begin
   ins = vesobj->containspoints(s.r(i,*),s.z(i,*))
   iw = where(ins eq 0,count)
   if count ne 0 then begin
      if n_elements(rhit) eq 0 then begin
         rhit = s.r(i,iw(0))
         zhit = s.z(i,iw(0))
         phit = s.phi(i,iw(0))
         pnhit = s.psin(i,iw(0))
      endif else begin
         rhit = [rhit,s.r(i,iw(0))]
         zhit = [zhit,s.z(i,iw(0))]
         phit = [phit,s.phi(i,iw(0))]
         pnhit = [pnhit,s.psin(i,iw(0))]
      endelse
   endif
;      endif else print,'Line ',i,' didn''t hit wall'
endfor

for j=1,ntor-1 do begin
   s=0.0d
   s=follow_fieldlines_phi(g,rmp,pnstart,pnstart,nsurfs,ntheta,phistart(j),-!dpi/180,nsteps,diffuse=diffuse,sigtheta=sigtheta,sigphi=sigphi)
   for i=0,nsurfs*ntheta-1 do begin
      ins = vesobj->containspoints(s.r(i,*),s.z(i,*))
      iw = where(ins eq 0,count)
      if count ne 0 then begin
         if n_elements(rhit) eq 0 then begin
            rhit = s.r(i,iw(0))
            zhit = s.z(i,iw(0))
            phit = s.phi(i,iw(0))
            pnhit = s.psin(i,iw(0))
         endif else begin
            rhit = [rhit,s.r(i,iw(0))]
            zhit = [zhit,s.z(i,iw(0))]
            phit = [phit,s.phi(i,iw(0))]
            pnhit = [pnhit,s.psin(i,iw(0))]
         endelse
      endif
;      endif else print,'Line ',i,' didn''t hit wall'
   endfor
endfor

print,n_elements(rhit),' out of ',ntor*nsurfs*ntheta,' lines hit wall'

;psi2d = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag)
rtest = dindgen(1000)/990 + g.rmaxis
ztest = 0.0*rtest + g.zmaxis
psitest = get_psi_bicub(g,rtest,ztest)
pntest = (-psitest-g.ssimag)/(g.ssibry-g.ssimag)
dpdr = deriv(rtest,pntest)
ip=where(dpdr lt 0.0 and pntest gt 1.0,count)
if count gt 0 then begin 
   rtest = rtest(0:ip(0))
   pntest = pntest(0:ip(0))
endif
romp = interpol(rtest,pntest,pnhit) - interpol(rtest,pntest,1.0d)

if not keyword_set(nbins) then nbins = 10
dr = (max(romp))/(nbins-1)
rbin = dindgen(nbins)*dr
rbin = [rbin(0)-2*dr,rbin(0)-dr,rbin]
nbins = nbins+2
hbin = dblarr(nbins)
for i=0,nbins-1 do begin
   junk = where(romp ge rbin(i) and romp lt rbin(i)+dr,count)
   if count gt 0 then hbin(i)=count else hbin(i)=0
endfor

;stop


return,{rbin:rbin+0.5*dr,hbin:hbin,rhit:rhit,phit:phit,zhit:zhit,romp:romp,pnhit:pnhit}
end

function process_spiral_diffuse,d,nphibin=nphibin,nrinbin=nrinbin,rinstart=rinstart,rinend=rinend,period=period,plot=plot,ps=ps

if not keyword_set(nphibin) then nphibin=10
if not keyword_set(nrinbin) then nrinbin=10
if not keyword_set(rinstart) then rinstart=1.02d
if not keyword_set(rinend) then rinend=1.13d
if not keyword_set(period) then period=1.0d

rinbin = (rinend-rinstart)*dindgen(nrinbin)/(nrinbin-1) + rinstart
phibin = (2*!dpi/period)*dindgen(nphibin)/(nphibin-1)

hinbin=dblarr(nphibin-1,nrinbin-1)
rincen = 0.5*(rinbin(0:nrinbin-2)+rinbin(1:nrinbin-1))
phicen = 0.5*(phibin(0:nphibin-2)+phibin(1:nphibin-1))

for ir = 0,nrinbin-2 do begin
   for ip = 0,nphibin-2 do begin
      ptmp = d.phit mod (2*!dpi/period)
      ptmp = (ptmp + 2*!dpi) mod (2*!dpi/period)
      junk = where((d.rhit ge rinbin(ir)) and (d.rhit lt rinbin(ir+1)) and (ptmp ge phibin(ip)) and (ptmp lt phibin(ip+1)),count)
      if count gt 0 then hinbin(ip,ir)=count else hinbin(ip,ir) = 0
   endfor
endfor

if keyword_set(plot) then begin

   if keyword_set(ps) then begin
      set_plot,'ps'
      !p.font=0
      device,file='plot_spiral_diffuse.ps',/color,/inches,xsize=5,ysize=4,xoff=.5,yoff=.5,/helvetica,/bold
   endif
   make_nice_plots,ps=ps

   contour,hinbin,phicen/!dpi,rincen,nlevel=50,/fill,xrange=[0,2/period],yrange=[rinstart,rinend],xstyle=1,ystyle=1,xtitle=textoidl('\phi/\pi'),ytitle='R (m)'

   if keyword_set(ps) then begin
      device,/close
      set_plot,'x'
   endif
   !p.position=0
endif


return,{rinbin:rinbin,rincen:rincen,phibin:phibin,phicen:phicen,hinbin:hinbin}
end


pro check_n_converge,g,rmp

ntest = [3,6,12,24,48,96,194]
n_ns = n_elements(ntest)
mmax = 60

btest = dblarr(n_ns,2*mmax+1)
cf=color_flipper()
make_nice_plots

pos = plot_positions(nrow=2,ncol=1,vertspace=0.1)

!p.position=pos(*,0)
plot,[0,0],[0,0],/nodata,xrange=[-mmax,mmax],yrange=[0,3.e-4],xtitle='m'

for i=0,n_ns-1 do begin
    s=get_b_spectrum(g,rmp,nn=3,nphi=ntest(i),nsurf=10,mmax=mmax)
    btest(i,*) = s.br(8,*)
    oplot,s.m,btest(i,*),color=cf(i)
endfor

!p.position=pos(*,1)
plot,ntest,btest(*,mmax-5),xtitle='nphi',/noerase


!p.position=0

end


pro check_jrat_converge,g,mm,nn

;mtest=[50,75,200,350]
mtest=[10,50,200]
n_m = n_elements(mtest)

mmax=60

plot,[0,0],[0,0],xrange=[-mmax,mmax],yrange=[0,1],xtitle='m'

make_nice_plots
cf=color_flipper()

for i=0,n_m-1 do begin
    j=get_jparrat(g,mm,nn,thmult=mtest(i))
    s=get_b_spectrum(g,j,pnwant=j.pn+0.001,ntheta=4*mmax)
    ind = where(s.m eq mm)
    bmax = s.br(ind)
    oplot,s.m,s.br/bmax(0),color=cf(i)
endfor

end


pro check_nphi_converge,g,t,mm,nn

ntest = [8,16,32,64,128,256]*nn
n_n = n_elements(ntest)

bs=dblarr(n_n)
bc = bs

rmp=build_nstx_rwmcoils(0d3)
pres = find_pnrat(g,mm,nn)

for i=0,n_n-1 do begin
    s=get_b_spectrum(g,rmp,pnwant=pres+0.001,ntheta=1000,nphi=ntest(i))
    im = where(s.m eq mm)
    bs(i) = s.br_s(im)
    bc(i) = s.br_c(im)
endfor

stop


end


function flux_ave_jtor,g,ntheta=ntheta

if not keyword_set(ntheta) then ntheta = 100
nline = 1e4
rs = dblarr(g.mw,ntheta)
zs = rs
theta = 2*!dpi*dindgen(ntheta)/(ntheta-1)
bp = rs
jtor = rs

for i=0,ntheta-1 do begin
    Lmax = max([g.r(g.mw-1)-g.r(0),g.z(g.mh-1) - g.z(0)])
    rline = g.rmaxis+Lmax*cos(theta(i))*dindgen(nline)/(nline-1)
    zline = g.zmaxis+Lmax*sin(theta(i))*dindgen(nline)/(nline-1)
    psiline = (-g.cpasma/abs(g.cpasma)*get_psi_bicub(g,rline,zline)-g.ssimag)/(g.ssibry-g.ssimag)
    ipn1 = where(psiline(0:nline-2) le 1.001d and psiline(1:nline-1) gt 1.001d,count)
    if count ge 1 then ipn1 = ipn1(0)
    iline = interpol(dindgen(ipn1+1),psiline(0:ipn1),g.pn)
    rs(*,i) = interpolate(rline,iline)
    zs(*,i) = interpolate(zline,iline)
    b=bfield_geq_bicub(g,rs(*,i),zs(*,i))
    bp(*,i) = sqrt(b.br^2+b.bz^2)
    jtor(*,i) = rs(*,i)*g.pprime + 4*!dpi*1d-7*g.ffprim/rs(*,i)
endfor

dl = sqrt((rs(*,1:ntheta-1)-rs(*,0:ntheta-2))^2 + (zs(*,1:ntheta-1)-zs(*,0:ntheta-2))^2)

dlob = total(dl/bp,2)
dlob(0) = 0.0
jdlob = total(jtor*dl/bp,2)
jdlob(0) = 0.0

jave = jdlob/dlob

dvdpsi = dlob/(2*!dpi)
dpsi = (g.pn(1)-g.pn(0))*(g.ssibry-g.ssimag)
vol = total(dvdpsi*dpsi,/cumulative)

stop



return,{rs:rs,zs:zs}
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function get_punct_midplane,g,s,sign=sign,plot=plot,ps=ps
if not keyword_set(sign) then sign=-1

nfl = n_elements(s.r(*,0))
npts = n_elements(s.r(0,*))

for i=0,nfl-1 do begin
   if sign eq 1 then iomp = where(s.z(i,0:npts-2) lt 0.0d and s.z(i,1:npts-1) gt 0.0d,count) else iomp = where(s.z(i,0:npts-2) gt 0.0d and s.z(i,1:npts-1) lt 0.0d,count)
   for j=0,count-1 do begin
      rtmp = interpol([s.r(i,iomp(j)),s.r(i,iomp(j)+1)],[s.z(i,iomp(j)),s.z(i,iomp(j)+1)],0.0d)
      if n_elements(romp) eq 0 then romp=rtmp else romp=[romp,rtmp]
      phitmp = interpol([s.phi(i,iomp(j)),s.phi(i,iomp(j)+1)],[s.z(i,iomp(j)),s.z(i,iomp(j)+1)],0.0d)
      if n_elements(phiomp) eq 0 then phiomp=phitmp else phiomp=[phiomp,phitmp]
   endfor
endfor

if keyword_set(plot) then begin
   make_nice_plots,ps=ps
   z=mycolors()
   plot,phiomp,romp,psym=3
endif

return,{r:romp,phi:phiomp}
end


function get_mercier_luc_u,g,r,z,delta=delta

if not keyword_set(delta) then delta = sqrt((g.r(1)-g.r(0))^2 + (g.z(1)-g.z(0))^2)

b1=bfield_geq_bicub(g,r,z)

norm_R = -b1.bz
norm_Z = b1.br
norm_M = sqrt(norm_R^2 + norm_Z^2)

;u = acos(norm_R/norm_M)
u = atan(norm_Z,norm_R)

bp=bfield_geq_bicub(g,r-delta*sin(u),z+delta*cos(u))
bpol_p = sqrt(bp.br^2 + bp.bz^2)

np_R = -bp.bz
np_Z = bp.br
np_M = sqrt(np_R^2 + np_Z^2)

;up = acos(np_R/np_M)
up = atan(np_Z,np_R)

bm=bfield_geq_bicub(g,r+delta*sin(u),z-delta*cos(u))
bpolm = sqrt(bm.br^2+bm.bz^2)

nm_R = -bm.bz
nm_Z = bm.br
nm_M = sqrt(nm_R^2 + nm_Z^2)

;um = acos(nm_R/nm_M)
um = atan(nm_Z,nm_R)

ip = where(abs(up-u) gt !dpi,count)
if count gt 0 then up(ip) = up(ip) - 2*!dpi*(up(ip)-u(ip))/abs(up(ip)-u(ip))
im = where(abs(um-u) gt !dpi,count)
if count gt 0 then um(im) = um(im) - 2*!dpi*(um(im)-u(im))/abs(um(im)-u(im))

;print,'Br,Bz,u',b1.br,b1.bz,u
;print,'rp,zp,up',r-delta*sin(u),z+delta*cos(u),up
;print,'rp,zp,up',r+delta*sin(u),z-delta*cos(u),um

dudl = 0.5d*(up-um)/delta
dbpdl = 0.5d*(bpol_p-bpol_m)/delta

drhodr = norm_R/norm_M
drhodz = norm_Z/norm_M

;print,'Rc= ',dudl

return,{u:u,rc:dudl,drhodr:drhodr,drhodz:drhodz,bpol:sqrt(b1.br^2+b1.bz^2),bphi:b1.bphi,dbpdl:dbpdl}
end



function bfield_mercluc,g,r,z,delta=delta

u = get_mercier_luc_u(g,r,z,delta=delta)

psi = get_psi_bicub(g,r,z)
psin = (-psi-g.ssimag)/(g.ssibry-g.ssimag)
pprime = interpol(g.pprime,g.pn,psin)
;fpol = interpol(g.fpol,g.pn,psin)
ffprime = interpol(g.ffprim,g.pn,psin)

psi1 = r*u.bpol
psi2 = 0.5d*(u.bpol*cos(u.u) - u.bpol*r/u.rc - 4*!dpi*r^2*pprime - ffprime)

;dp1dr = 

return,0
end


function culssg3d,g,ps=ps

sloc = dblarr(g.mw,g.mh)
cks = sloc
ckth = sloc
psi2d = sloc

;b = bfield_geq_bicub(g,g.rmaxis,g.zmaxis)
;b=replicate(b,[g.mw,g.mh])

for iz=0,g.mh-1 do begin
   b = bfield_geq_bicub(g,g.r,0*g.r + g.z(iz))
   bth = sqrt(b.br^2+b.bz^2)
   bb = sqrt(b.br^2+b.bz^2+b.bphi^2)
   psi = get_psi_bicub(g,g.r,0*g.r + g.z(iz),/secderiv)
   psi.psi *= -1
   psi.dsdr *= -1
   psi.dsdz *= -1
   psi.d2sdr2 *= -1
   psi.d2sdz2 *= -1
   psi.d2sdrdz *= -1
   psin = (psi.psi - g.ssimag)/(g.ssibry-g.ssimag)
   psi2d(*,iz) = psin
   fpol = interpol(g.fpol,g.pn,psin)
   ffp = interpol(g.ffprim,g.pn,psin)
   fprim = ffp/fpol

   sloc(*,iz) = (fpol/(g.r^4*bth^2*bb^2))*((psi.d2sdr2-psi.d2sdz2)*(psi.dsdr^2-psi.dsdz^2)+4*psi.dsdr*psi.dsdz*psi.d2sdrdz) $
          + fpol*psi.dsdr/(g.r^3*bb^2) - fprim*bth^2/(bb^2)
   cks(*,iz) = (fpol^2*psi.dsdr + g.r*(psi.d2sdr2*psi.dsdz^2 + psi.dsdr^2*psi.d2sdz2 - 2*psi.dsdr*psi.d2sdrdz*psi.dsdz))/(g.r^4*bth*bb^2)
   ckth(*,iz) = -fpol*(g.r*(psi.d2sdrdz*(psi.d2sdr2-psi.d2sdz2)-psi.dsdr*psi.dsdz*(psi.d2sdr2-psi.d2sdz2))+psi.dsdz*g.r^2*bb^2)/(g.r^5*bth*bb^3)
endfor

rgrid2d = reform(g.r # replicate(1,g.mh),g.mw*g.mh)
zgrid2d = reform(g.z ## replicate(1,g.mw),g.mw*g.mh)
ig = where(g.bdry(0,*) gt 1e-4)
vesobj = obj_new('IDLanROI',g.bdry(*,ig))
ibad = where(vesobj->containspoints(rgrid2d,zgrid2d) eq 0)

sloc(ibad) = 0
cks(ibad) = 0
ckth(ibad) = 0

if keyword_set(ps) then begin
   set_plot,'ps'
   !p.font=0
   device,file='culssg3d.ps',/inches,/color,xsize=4,ysize=7,xoff=.5,yoff=.5,/helvetica,/bold,bits_per_pixel=8
endif

make_nice_plots,ps=ps

pos = plot_positions(nrow=3,ncol=1,vertspace=0.1)
pos(2,*) = 0.8
;print,pos(*,1)
cpos = pos
cpos(0,*) = 0.9
cpos(2,*) = 0.99

;r = g.r + randomn(seed,n_elements(is))*1d-8
;z = g.z + randomn(seed,n_elements(is))*1d-8

!p.position=pos(*,0)
levs = (max(sloc)-min(sloc))*dindgen(50)/49 + min(sloc)
contour,sloc,g.r,g.z,levels=levs,/fill,xtitle=textoidl('R'),ytitle=textoidl('Z'),$
        title=textoidl('Local shear'),charsize=0.8,xcharsize=1.1,ycharsize=1.1,/iso
colorbar,range=[min(levs),max(levs)],/vertical,format='(F5.2)',position=cpos(*,0),title=textoidl('s_{loc}'),charsize=0.7

!p.position=pos(*,1)
levs = (max(cks)-min(cks))*dindgen(50)/49 + min(cks)
contour,cks,g.r,g.z,levels=levs,/fill,xtitle=textoidl('R'),ytitle=textoidl('Z'),$
        title=textoidl('Normal curvature'),/noerase,charsize=0.8,xcharsize=1.1,ycharsize=1.1,/iso
colorbar,range=[min(levs),max(levs)],/vertical,format='(F5.2)',position=cpos(*,1),title=textoidl('\kappa_s'),charsize=0.7

!p.position=pos(*,2)
levs = (max(ckth)-min(ckth))*dindgen(50)/49 + min(ckth)
contour,ckth,g.r,g.z,levels=levs,/fill,xtitle=textoidl('R'),ytitle=textoidl('Z'),$
        title=textoidl('Geodesic curvature'),/noerase,charsize=0.8,xcharsize=1.1,ycharsize=1.1,/iso
colorbar,range=[min(levs),max(levs)],/vertical,format='(F5.2)',position=cpos(*,2),title=textoidl('\kappa_{\theta}'),charsize=0.7

if keyword_set(ps) then begin
   device,/close
   set_plot,'x'
endif


return,{r:g.r,z:g.z,sloc:sloc,cks:cks,ckth:ckth}
end


pro plot_toroidal_current,g,xrange=xrange,yrange=yrange,iso=iso,ps=ps,fill=fill,neg=neg,wout=wout,diffscale=diffscale,xstyle=xstyle,ystyle=ystyle

r2d = dblarr(g.mw,g.mh)
z2d = r2d
bb2d=r2d
for i =0,g.mh-1 do r2d(*,i) = g.r
for j= 0,g.mw-1 do z2d(j,*) = g.z

if not keyword_set(xrange) then xrange=[1.0,2.4]
if not keyword_set(yrange) then yrange=[-1.5,1.0]

il = where(g.bdry(0,*) gt 1e-4)
sepobj = obj_new('IDLanROI',g.bdry(0,il),g.bdry(1,il))
ins = sepobj->containspoints(r2d,z2d)
ins = reform(ins,[g.mw,g.mh])

for i=0,g.mh-1 do begin
   b=bfield_geq_bicub(g,r2d(*,i),z2d(*,i))
   bb=sqrt(b.br^2+b.bz^2+b.bphi^2)
   bb2d(*,i)=bb
endfor

bb2d *= ins

psi2d = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag)

pp2d = ins*interpol(g.pprime,g.pn,psi2d)
ffp2d = ins*interpol(g.ffprim,g.pn,psi2d)
mu0 = 4*!dpi*1d-7

jt2d = r2d*pp2d + ffp2d/(r2d*mu0)
if keyword_set(neg) then jt2d *= -1.0d
jt2d_not0 = jt2d(where(abs(jt2d) gt 1d-6))
jtmin=min(jt2d_not0)
jtmax=max(jt2d_not0)
;jt2d *=sepobj->containspoints(r2d,z2d)
;jt2d = reform(jt2d,[g.mw,g.mh])
jttot = total(jt2d*(g.r(1)-g.r(0))*(g.z(1)-g.z(0)))
print,'Jtot_2D, Jtot_efit'
print,jttot,g.cpasma
print,'Min, Max Jtor'
print,jtmin,jtmax

if keyword_set(ps) then begin
   set_plot,'ps'
   !p.font=0
   device,file='jtor_2d.ps',/inches,/color,xsize=7,ysize=4,/helvetica,/bold,xoffset=0.5,yoffset=0.5
endif
nlev = 100
jlevs = (jtmax-jtmin)*dindgen(nlev)/(nlev-1) + jtmin
jcol = !d.table_size*dindgen(nlev)/(nlev-1)
make_nice_plots,ps=ps
z=mycolors()
pos = plot_positions(nrow=1,ncol=2)
;if not keyword_set(ps) then window,0
!p.position=pos(*,0)
contour,jt2d,g.r,g.z,level=jlevs,c_color=jcol,xrange=xrange,yrange=yrange,iso=iso,xtitle='R (m)',ytitle='Z (m)', $
        title=textoidl('j_{tor}: EFIT'),thick=0.5,fill=fill,xstyle=xstyle,ystyle=ystyle
legend,['j^{MIN}='+string(jtmin,format='(E9.2)'),'j^{MAX}='+string(jtmax,format='(E9.2)')],pos(*,0),color=[z.black,z.black],corner='br',spacing=0.05,offsetx=-0.055,offsety=0.01
;oplot,g.bdry(0,*),g.bdry(1,*),color=z.brick,thick=4

;if not keyword_set(ps) then window,1
if not keyword_set(wout) then begin
   !p.position=pos(*,1)
   contour,bb2d,g.r,g.z,nlevel=50,xrange=xrange,yrange=yrange,iso=iso,xtitle='R (m)',ytitle='Z (m)', $
           title=textoidl('|B|'),thick=0.5,/noerase,xstyle=xstyle,ystyle=ystyle
endif else begin
   jt_vmec = get_vmec_jtor(wout)  
   if keyword_set(diffscale) then jlevs = (jt_vmec.jtmax-jt_vmec.jtmin)*dindgen(nlev)/(nlev-1) + jt_vmec.jtmin
   !p.position=pos(*,1)
   contour,jt_vmec.jtor,jt_vmec.r,jt_vmec.z,level=jlevs,c_color=jcol,xrange=xrange,yrange=yrange,iso=iso, $
           xtitle='R (m)',ytitle='Z (m)', title=textoidl('j_{tor}: VMEC'),thick=0.5,fill=fill,/noerase,xstyle=xstyle,ystyle=ystyle
   legend,['j^{MIN}='+string(jt_vmec.jtmin,format='(E9.2)'),'j^{MAX}='+string(jt_vmec.jtmax,format='(E9.2)')],pos(*,1),color=[z.black,z.black],corner='br',spacing=0.05,offsetx=-0.055,offsety=0.01
;   if keyword_set(diffscale) then begin
;      legend,['DIFFERENT CONTOUR LEVELS'],pos(*,1),color=z.brick,corner='bl'
;   endif else begin
;      legend,['SAME CONTOUR LEVELS AS EFIT'],pos(*,1),color=z.brick,corner='bl'
;   endelse
endelse
!p.position=0

if keyword_set(ps) then begin
   device,/close
   set_plot,'x'
endif

end


pro plot_punct,g,s,ps=ps,poincare=poincare,yrange=yrange

if not keyword_set(yrange) then yrange=[0,1]
if not keyword_set(xrange) then xrange=[0,1]

if keyword_set(ps) then begin
   set_plot,'ps'
   !p.font=0
   device,file='punct.ps',xsize=5,ysize=5,xoff=.5,yoff=.5,/helvetica,/bold,/inches
endif
   
if not keyword_set(poincare) then sqrts = sqrt(interpol(g.psitor,g.pn,s.psin)) else sqrts=s.s

make_nice_plots,ps=ps

if keyword_set(poincare) then begin
   plot,s.theta/(2*!dpi),s.s,psym=1,symsize=0.3,xtitle=textoidl('\theta/2\pi'),ytitle=textoidl('\psi_t^{1/2}'),xrange=xrange,yrange=yrange,xstyle=1,ystyle=1
endif else begin
   plot,s.theta(*,s.ip0)/(2*!dpi),sqrts(*,s.ip0),psym=1,symsize=0.3,xtitle=textoidl('\theta/2\pi'),ytitle=textoidl('\psi_t^{1/2}'),xrange=xrange,yrange=[0.5,1.],xstyle=1,ystyle=1
endelse


if keyword_set(ps) then begin
   device,/close
   set_plot,'x'
endif


end
