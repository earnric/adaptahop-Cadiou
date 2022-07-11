pro read_halo,filegal=filegal,fileDM=fileDM,fileim=fileim,nout=nout,mass=mass,mvirial=mvirial,radius=radius,zoom=zoom,x0=x0,y0=y0,z0=z0,th=th,fileinfo=fileinfo,min=min,max=max,log=log,temp=temp,stars=stars,mmin=mmin,dir=dir

if not keyword_set(dir)then dir='x'
if not keyword_set(mmin)then mmin=1d9
if not keyword_set(zoom)then zoom=1d0
if not keyword_set(x0)then x0=0.5d0
if not keyword_set(y0)then y0=0.5d0
if not keyword_set(z0)then z0=0.5d0

if not keyword_set(th)then th=1

nbodies=0L & nb_of_halos=0L & nb_of_subhalos=0L 
nb_of_parts=0L
my_number=0L & my_timestep=0L & nbsub=0L
hosthalo=0L & hostsub=0L & level=0L
nextsub=0L

if keyword_set(fileinfo)then begin
    rd_info,info,file=fileinfo
endif else begin
    rd_info,info,nout=nout
endelse

xmin=(x0-0.5/zoom-0.5d0)*(info.unit_l/3.08d24)
xmax=(x0+0.5/zoom-0.5d0)*(info.unit_l/3.08d24)
ymin=(y0-0.5/zoom-0.5d0)*(info.unit_l/3.08d24)
ymax=(y0+0.5/zoom-0.5d0)*(info.unit_l/3.08d24)
zmin=(z0-0.5/zoom-0.5d0)*(info.unit_l/3.08d24)
zmax=(z0+0.5/zoom-0.5d0)*(info.unit_l/3.08d24)
x00=(x0-0.5d0)*(info.unit_l/3.08d24)
y00=(y0-0.5d0)*(info.unit_l/3.08d24)
z00=(z0-0.5d0)*(info.unit_l/3.08d24)
print,xmin,xmax,ymin,ymax,zmin,zmax
print,x00,y00,z00
rd_im,d,file=fileim
mH=1.66d-24
kB=1.38d-16
lbox=info.unit_l/3.08d24/zoom
if not keyword_set(min)then min=-5.0
if not keyword_set(max)then max=0.001
if keyword_set(temp)then begin
    d=d*mH/kB*info.unit_l^2d0/info.unit_t^2d0*5d-8
endif else if not keyword_set(stars)then begin
    d=d*info.unit_d/mH
endif
if keyword_set(log)then d=alog10(d)
mycontour,d,/table,ncont=50,lbox=lbox,unit='(Mpc)',title=' ',min=min,max=max

if keyword_set(filegal)then begin
openr,1,filegal,/f77_unformatted
readu,1,nbodies
readu,1,massp
readu,1,aexp
readu,1,omega_t
readu,1,age_univ
readu,1,nb_of_halos,nb_of_subhalos
print,'nbodies =',nbodies
print,'massp   =',massp
print,'aexp    =',aexp
print,'omega_t =',omega_t
print,'age_univ=',age_univ
print,'Nombre de galaxies     =',nb_of_halos
print,'Nombre de sous-galaxies=',nb_of_subhalos
nmax=nb_of_halos+nb_of_subhalos

for i=0L,nmax-1 do begin
    readu,1,nb_of_parts
    members=lonarr(nb_of_parts)
    readu,1,members

    ; Read properties of each halo
    readu,1,my_number
    readu,1,my_timestep
    readu,1,level,hosthalo,hostsub,nbsub,nextsub
    readu,1,m
    readu,1,px,py,pz
    readu,1,vx,vy,vz
    readu,1,Lx,Ly,Lz
    readu,1,r,sha,shb,shc
    readu,1,ek,ep,et
    readu,1,spin
    readu,1,sigma,sigmab,mb
    readu,1,rvir,mvir,tvir,cvel
    readu,1,rho_0,r_c

    tek_color
    if(m*1d11 gt mmin)then begin
    if (px lt xmax and py lt ymax and pz lt zmax and px gt xmin and py gt ymin and pz gt zmin) then begin
        if(dir eq 'x')then begin
            xc=py
            yc=pz
        endif 
        if(dir eq 'y')then begin
            xc=px
            yc=pz
        endif 
        if(dir eq 'z')then begin
            xc=px
            yc=py
        endif
        rc=rvir
        if keyword_set(radius)then rc=rvir
        if keyword_set(mass)then rc=(alog10(double(m*1d11   )/mmin)+1)*0.005d0*lbox
        if keyword_set(mvirial)then rc=(alog10(double(mvir*1d11)/mmin)+1)*0.005d0*lbox
        points = (2 * !PI / 99.0) * FINDGEN(100)
        x = xc - x00 + rc * COS(points )
        y = yc - y00 + rc * SIN(points )
        tek_color
        oplot,x,y,th=th,color=3
    endif
    endif

endfor
close,1
endif

if keyword_set(fileDM)then begin
openr,1,fileDM,/f77_unformatted
readu,1,nbodies
readu,1,massp
readu,1,aexp
readu,1,omega_t
readu,1,age_univ
readu,1,nb_of_halos,nb_of_subhalos
print,'nbodies =',nbodies
print,'massp   =',massp
print,'aexp    =',aexp
print,'omega_t =',omega_t
print,'age_univ=',age_univ
print,'Nombre de halos     =',nb_of_halos
print,'Nombre de sous-halos=',nb_of_subhalos
nmax=nb_of_halos+nb_of_subhalos

tek_color

for i=0L,nmax-1 do begin
    readu,1,nb_of_parts
    members=lonarr(nb_of_parts)
    readu,1,members

    ; Read properties of each halo
    readu,1,my_number
    readu,1,my_timestep
    readu,1,level,hosthalo,hostsub,nbsub,nextsub
    readu,1,m
    readu,1,px,py,pz
    readu,1,vx,vy,vz
    readu,1,Lx,Ly,Lz
    readu,1,r,sha,shb,shc
    readu,1,ek,ep,et
    readu,1,spin
    readu,1,rvir,mvir,tvir,cvel
    readu,1,rho_0,r_c

    if (px lt xmax and py lt ymax and pz lt zmax and px gt xmin and py gt ymin and pz gt zmin) then begin
        xc=px
        yc=py
        rc=rvir
        if keyword_set(radius)then rc=rvir
        if keyword_set(mass)then rc=(alog10(double(m*1d11   )/1d9)+1)*0.005d0*lbox
        if keyword_set(mvirial)then rc=(alog10(double(mvir*1d11)/1d9)+1)*0.005d0*lbox
        if(m*1d11 gt mmin)then begin
            points = (2 * !PI / 99.0) * FINDGEN(100)
            x = xc - x00 + rc * COS(points )
            y = yc - y00 + rc * SIN(points )
            oplot,x,y,th=th,color=2
            
            points = FINDGEN(100)/100d0-0.5d0
            x = xc - x00 + lbox/2d2 * points
            y = yc - y00 + lbox/2d2 * points
            oplot,x,y,th=th
            x = xc - x00 + lbox/2d2 * points
            y = yc - y00 - lbox/2d2 * points
            oplot,x,y,th=th
        endif
    endif

endfor
close,1
endif

loadct,0
end

