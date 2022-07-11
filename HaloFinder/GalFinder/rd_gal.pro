function getcarnumlong,nframes
ndigitmx = floor( alog10( 9999990 ) ) + 1
car = strarr(nframes)
for i = 1L, nframes do begin
       a = string(i)
       if (i gt 0) then $
         ndigit =  floor( alog10( float(i) ) ) + 1 $
       else $
          ndigit = 1
       for j = 1L, ndigitmx - ndigit do begin
           a = '0' + a
       endfor
       car(i-1) = strcompress(a, /remove_all)
;       print,car(i-1)
   endfor
return,car
end

pro rd_gal,gal,nout=nout,dir=dir,idgal=idgal,metal=metal,nchem=nchem,long=long

if not keyword_set(dir)then begin
   if not keyword_set(nout)then stop
   car=getcarnum(nout)
   dir='GAL_'+car(nout-1L)
endif
if not keyword_set(long)then begin
    cargal=getcarnum(idgal)
endif else begin
    cargal=getcarnumlong(idgal)
endelse
openr,9,dir+'/gal_stars_'+cargal(idgal-1L),/f77_unformatted
mynumber=0L&level=0L&mg=0d0&xg=0d0&yg=0d0&zg=0d0&vxg=0d0&vyg=0d0&vzg=0d0
Lxg=0d0&Lyg=0d0&Lzg=0d0&nlist=0L
readu,9,my_number
readu,9,level
readu,9,mg
readu,9,xg,yg,zg
readu,9,vxg,vyg,vzg
readu,9,Lxg,Lyg,Lzg
readu,9,nlist
mg=mg*1d11
xx     =dblarr(nlist)
xpart  =dblarr(nlist,3)
vpart  =dblarr(nlist,3)
mpart  =dblarr(nlist)
idpart =lonarr(nlist)
agepart=dblarr(nlist)
metpart=dblarr(nlist)
if keyword_set(nchem)then chempart=dblarr(nlist,nchem)
for idim=1,3 do begin
   readu,9,xx
   xpart(0L:nlist-1L,idim-1L)=xx
endfor
for idim=1,3 do begin
   readu,9,xx
   vpart(0L:nlist-1L,idim-1L)=xx
endfor
readu,9,mpart
readu,9,idpart
readu,9,agepart
if keyword_set(metal)then begin
   readu,9,metpart
   if keyword_set(nchem)then begin
      for ichem=1L,nchem do begin
         readu,9,xx
         chempart(0L:nlist-1L,ichem-1L)=xx
      endfor
   endif
endif

close,9
if keyword_set(metal)then begin
   if keyword_set(nchem)then begin
      gal={ mg:mg, xg:xg, yg:yg, zg:zg, vxg:vxg, vyg:vyg, vzg:vzg $
            ,Lxg:Lxg, Lyg:Lyg, Lzg:Lzg, level:level $
            ,npart:nlist $
            ,xp:fltarr(nlist,3) $
            ,vp:fltarr(nlist,3) $
            ,id:lonarr(nlist) $
            ,mp:fltarr(nlist) $
            ,zp:fltarr(nlist) $
            ,cp:fltarr(nlist,nchem) $
            ,ap:fltarr(nlist) }
   endif else begin
      gal={ mg:mg, xg:xg, yg:yg, zg:zg, vxg:vxg, vyg:vyg, vzg:vzg $
            ,Lxg:Lxg, Lyg:Lyg, Lzg:Lzg, level:level $
            ,npart:nlist $
            ,xp:fltarr(nlist,3) $
            ,vp:fltarr(nlist,3) $
            ,id:lonarr(nlist) $
            ,mp:fltarr(nlist) $
            ,zp:fltarr(nlist) $
            ,ap:fltarr(nlist) }      
   endelse
endif else begin
   gal={ mg:mg, xg:xg, yg:yg, zg:zg, vxg:vxg, vyg:vyg, vzg:vzg $
         ,Lxg:Lxg, Lyg:Lyg, Lzg:Lzg, level:level $
         ,npart:nlist $
         ,xp:fltarr(nlist,3) $
         ,vp:fltarr(nlist,3) $
         ,id:lonarr(nlist) $
         ,mp:fltarr(nlist) $
         ,ap:fltarr(nlist) }
endelse
gal.xp(0L:nlist-1L,0)=xpart(0L:nlist-1L,0)
gal.xp(0L:nlist-1L,1)=xpart(0L:nlist-1L,1)
gal.xp(0L:nlist-1L,2)=xpart(0L:nlist-1L,2)
gal.vp(0L:nlist-1L,0)=vpart(0L:nlist-1L,0)
gal.vp(0L:nlist-1L,1)=vpart(0L:nlist-1L,1)
gal.vp(0L:nlist-1L,2)=vpart(0L:nlist-1L,2)
gal.mp(0L:nlist-1L)=mpart(0L:nlist-1L)
gal.ap(0L:nlist-1L)=agepart(0L:nlist-1L)
gal.id(0L:nlist-1L)=idpart(0L:nlist-1L)
if keyword_set(metal)then begin
   gal.zp(0L:nlist-1L)=metpart(0L:nlist-1L)
   if keyword_set(nchem)then begin
      for ichem=1L,nchem do begin
         gal.cp(0L:nlist-1L,ichem-1L)=chempart(0L:nlist-1L,ichem-1L)
      endfor
   endif
endif

end
