function make_triangle,polyarr
  nt = 0L
  ip = 0L
  while ( ip lt n_elements(polyarr)) do begin
    nv = polyarr[ip]
    nt = nt+nv-2
    ip = ip+nv+1
  endwhile

  res = lonarr(3,nt)

  ip = 0L
  it = 0L
  while ( ip lt n_elements(polyarr)) do begin
    nv   = polyarr[ip]
    poly = polyarr[ip+1:ip+nv]
    ip   = ip+nv+1
    for iv=1L,n_elements(poly)-2 do begin
      res[*,it] = [poly[0],poly[iv],poly[iv+1]]
      it = it+1
    endfor
  endwhile
 
  return,res
end   

pro swap,a,b
  tmp=a
  a=b
  b=tmp
end

function cross,a,b
  return,[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
end

function scalar,a,b
  ta = a
  tb = b
  na = n_elements(ta[0,*])
  nb = n_elements(tb[0,*])

  if (na eq nb) then begin
    res = fltarr(na)
    for i=0L,na-1 do res[i] = total(ta[*,i]*tb[*,i])
    return, res
  endif
  if (nb gt na) then begin
    swap,ta,tb
    swap,na,nb    
  endif
  res = fltarr(na)
  for i=0L,na-1 do res[i] = total(ta[*,i]*tb)
  return, res
end

function normal,vec
  nt = n_elements(vec)/3
  res = vec
  for i=0L,nt-1 do begin
    res[*,i] = res[*,i]/sqrt(total(res(*,i)^2))
  endfor
  return,res
end

function calc_normal,vert,triarr
  nt = n_elements(triarr)/3
  res = fltarr(3,nt)

  for i=0L,nt-1 do begin
    v1 = vert[*,triarr[1,i]]-vert[*,triarr[0,i]]
    v2 = vert[*,triarr[2,i]]-vert[*,triarr[1,i]]
    res[*,i] = cross(v1,v2)
    res[*,i] = res[*,i]/sqrt(total(res(*,i)^2))
  endfor

  return,res
end

function calc_pos,vert,tri
  nt = n_elements(tri)/3
  res = fltarr(3,nt)
  
  for i=0L,nt-1 do begin
    res[*,i] = (vert[*,tri[0,i]]+vert[*,tri[1,i]]+vert[*,tri[2,i]])/3
  endfor
  
  return,res
end

function sub,a,b
  ta = a
  tb = b
  na = n_elements(ta[0,*])
  nb = n_elements(tb[0,*])

  if (na eq nb) then begin
    return, ta-tb
  endif
  if (na gt nb) then begin
    res = ta
    for i=0L,na-1 do res[*,i] = ta[*,i]-tb[*,0]
    return, res
  endif
  res = tb
  for i=0L,nb-1 do res[*,i] = ta[*,0]-tb[*,i]
  return, res
end

pro mypoly,vert,poly,light=light,ambient=ambient,new=new
  common mypoly,tri,pos,norm,pol

  if (n_elements(pol) eq 0) then pol = 0
  if (not keyword_set(light)) then light=[-1,-1,1]*10000
  if (not keyword_set(ambient)) then ambient=0.15
  if (not keyword_set(new)) then new=0

  mincol = ambient*255
  maxcol = 230

  view = (!P.T[0:2,0:2])#[0,0,-1]
 
  if (new or total(poly[0:min([100,n_elements(poly)-1])]-pol) ne 0) $
  then begin
    pol  = poly[0:min([100,n_elements(poly)-1])]
    tri  = make_triangle(poly)
    pos  = calc_pos(vert,tri)
    norm = calc_normal(vert,tri)
  endif

  shade=bytscl(scalar(norm,normal(sub(pos,light))),min=0,$
               top=maxcol-mincol)+mincol

  ind = sort(scalar(pos,-view))
  ind2 = where(scalar(norm(*,ind),view) ge 0)
  ind = ind(ind2)

  for i=0L,n_elements(ind)-1 do begin 
    polyfill,vert(*,tri(*,ind(i))),/t3d,color=shade(ind(i))
  endfor
end
