pro rr,name
;procedure to call resolve_routine and return to calling program in event of error
catch, ierror
if ierror ne 0 then return
resolve_routine,name
return
end

function ri,name,source=source

;procedure to call routine_info for a known function
;return to calling program in event of error
a={path:'Function not found',num_args:0, num_kw_args:0, args:0, kw_args:0}
catch, ierror
if ierror ne 0 then return,a
if not keyword_set(source) then a=routine_info(name,/param,/func) else $
  a=routine_info(name,/source,/func)
return,a
end

function rip,name,source=source

;procedure to call routine_info for a known procedure
;return to calling program in event of error
a={path:'Function not found',num_args:0, num_kw_args:0, args:0, kw_args:0}
catch, ierror
if ierror ne 0 then return,a
if not keyword_set(source) then a=routine_info(name,/param) else $
  a=routine_info(name,/source)
return,a
end

pro get_info,name,help=help

; Purpose - retrieve various information on a routine
; Written ~ 1997 by R. Maingi and modified several times

name=strlowcase(name)
title='PROCEDURE: '
if name eq '' or keyword_set(help) then begin
  print,'program get_info(name,/help),
  print,' This routine returns the arguments and keywords of '
  print,' argument NAME. The help keyword prints out this help'
  print,' page. Written by R. Maingi 2/18/98'
  print,'Keywords: func, if the routine is a function'
  return
endif
uname=strupcase(name)
if not keyword_set(func) then func=0
if max(strpos(routine_info(),uname) eq 0) eq 0 and $ 
   max(strpos(routine_info(/func),uname) eq 0) eq 0 then rr,name
if max(strpos(routine_info(),uname) eq 0) eq 0 and $
   max(strpos(routine_info(/func),uname) eq 0) eq 0 then $
  print,'get_info: could not find '+name+' in IDL path, returning' $
else begin
  a=rip(name) ;check list of compiled procedures
  if max(a.num_args) eq 0 and max(a.num_kw_args) eq 0 then begin 
    a=ri(name) ;probably a function
    b=ri(name,/source)
    title='FUNCTION: '
  endif else b=routine_info(name,/source)
  print,title+b.path
  if max(tag_names(a) eq 'ARGS') ge 1 then print,'Arguments: ',a.args else print,'Arguments: None'
  if max(tag_names(a) eq 'KW_ARGS') ge 1 then begin
    print,'Keywords: ',a.kw_args 
    if max(a.kw_args eq 'HELP') ge 1 then print,'Note: HELP keyword is available for this routine'
  endif else print,'Keywords: None'
endelse

return
end
