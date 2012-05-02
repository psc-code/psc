import matplotlib.pyplot as plt
import numpy as np
plt.spectral()#setting the default color map to "spectral"
plt.rcParams['contour.negative_linestyle'] = 'solid'

def getfile(ts=None,dtype=None,filename=None,fmt="%s.%06d_p000000.h5"):
	'''
	Return the hdf5 file pointer. "filename" is either given explicitly, or formatted from "ts"(timestep) and "dtype"('pfd'/'tfd') according to "fmt".
	'''
	import h5py,os
	if filename==None:
		if ts!=None and dtype!=None:
			filename=fmt%(dtype,ts)
		else:
			print "ERROR: either ts & dtype together, or filename alone should be given"
			return
	try:
		f=h5py.File(filename,'r')
	except IOError:
		print "ERROR: Cannot open %s"%filename
		return
	return f

def getfld_file(fld,f,dims='xz',n=0,i1b=0,i1e=None,i2b=0,i2e=None,stride1=None,stride2=None,stride=1):
	'''
	Retrieve field named by "fld" from a hdf5 file pointer. For the first dimension(e.g., 'x' in 'xz'), indices ;i1b'/'i1e' specify a range from which data is retrieved. Data array will also be strided every 'stride1' step. If 'stride1' is not set, the value of 'stride' will be used, which is 1 (no striding) by default.	Similarly, for the second dimension, 'i2b', 'i2e', 'stride2' can be used.
	'''
	if stride1==None:
		stride1=stride
	if stride2==None:
		stride2=stride
	try:
		_field_=f[fld][fld]['p0']['3d']
	except KeyError:
		print "Error: Cannot open object. Possibly %s is not stored in the %s."%(fld,f.filename)
		return
	if dims=='xz':
		field=_field_[i2b:i2e:stride2,n,i1b:i1e:stride1]
	elif dims=='yz':
		field=_field_[i2b:i2e:stride2,i1b:i1e:stride1,n]
	elif dims=='xy':
		field=_field_[n,i2b:i2e:stride2,i1b:i1e:stride1]
	if fld=='ne':
		return -field
	return field

def getfld(fld,ts=None,dtype=None,filename=None,**kwargs):
	'''
	Retrieve field named by 'fld' from a PSC data file whose name is either given explicitly by 'filename', or formatted by 'ts'(timestep) and 'dtype'('pfd' or 'tfd').
	'**kwargs' will be passed to getfld_file, which includes: 
		'dims' specifies which two dimensions are wanted(one of 'xz', 'yz', and 'xy'; default is 'xz'). 
		For the first dimension(e.g., 'x' in 'xz'), indices ;i1b'/'i1e' specify a range from which data is retrieved. Data array will also be strided every 'stride1' step. If 'stride1' is not set, the value of 'stride' will be used, which is 1 (no striding) by default.	Similarly, for the second dimension, 'i2b', 'i2e', 'stride2' can be used.
	'''
	f = getfile(ts, dtype, filename)
	field = getfld_file(fld,f,**kwargs)
	f.close()
	return field

def v2d(field,clim=None,show=True,**kwargs):
	'''
	Visualize the 2d array using plt.imshow. '**kwargs' will be passed to plt.imshow, where most useful one is 'extent' specifies the extent of the two coordinates. 
	'clim'=[-cmin,-cmax] gives the lower and upper limit of the colormap. For example, setting clim so that zero is in the middle.
	The boolean 'show' determines if the figure is to be shown right (for example, by setting show=False you can keep the figure and plot magnetic field line on top of it, and then show/save it by plt.savefig).
	'''
	plt.imshow(field,origin='lower',**kwargs)
	plt.colorbar()
	plt.clim(clim)
	if show:
		plt.show()

def view2d(fld,ts=None,dtype=None,filename=None,clim=None,show=True,extent=None,**kwargs):
	'''Retrieve and view the data for 'fld'. '**kwargs' is passed to getfld.'''
	field=getfld(fld,ts,dtype,filename,**kwargs)
	v2d(field,clim,show,extent=extent)

def v1d(field,dim,n=None,r=None,extent=None,show=True,**kwargs):
	'''View a slice view. 'dim'=1 or 2 for the first or second dimension. 'n' gives the index of the slice in the other dimension. '**kwargs' is passed to plt.plot.'''
	if extent==None:
		plt.plot(slicing(field,dim,n,r),**kwargs)
	else:
		n2,n1=field.shape
		if dim==1:
			length=n1
		else:
			length=n2
		x=np.linspace(extent[0],extent[1],length,endpoint=False)
		plt.plot(x,slicing(field,dim,n,r),**kwargs)
	plt.axis('tight')
	if show:
		plt.show()

'''
Calculate physical quantities that are not stored directly in the PSC data files.
'''
def calc_pressure(kind,dim1,dim2,f,**kwargs):
	'''Calculate pressure tensor n<vv>-n<v><v>.'''
	fn   = getfld_file('n%c' %kind,f,**kwargs)
	fv1  = getfld_file('v%c%c'%(kind,dim1),f,**kwargs)
	fv2  = getfld_file('v%c%c'%(kind,dim2),f,**kwargs)
	fv12 = getfld_file('v%c%cv%c%c'%(kind,dim1,kind,dim2), f,**kwargs)
	return fvxy - fvx*fvy/fn

def calc_u(kind,dim,f,**kwargs):
	'''
	Calculate flow velocities. kind = 'e', 'i', or 'all'; dim='x', 'y', or 'z'
	'''
	if kind in ['e','i']:
		fn  = getfld_file('n%c'%kind,f,**kwargs)
		fnv = getfld_file('v%c%c'%(kind,dim),f,**kwargs)
		return fnv/fn
	elif kind =='all':
		fne  = getfld_file('ne',f,**kwargs)
		fnve = getfld_file('ve%c'%dim,f,**kwargs)
		fni  = getfld_file('ni',f,**kwargs)
		fnvi = getfld_file('vi%c'%dim,f,**kwargs)
		return (fnve+fnvi)/(fne+fni)

def calc_2(fld,f,**kwargs):
	'''
	Calculate square of a quantity. fld = 'e', 'h'
	'''
	if fld in ['e','h']:
		fx  = getfld_file(fld+'x',f,**kwargs)
		fy  = getfld_file(fld+'y',f,**kwargs)
		fz  = getfld_file(fld+'z',f,**kwargs)
	elif fld in ['ue','ui']:
		fx  = calcfld_file(fld+'x',f,**kwargs)
		fy  = calcfld_file(fld+'y',f,**kwargs)
		fz  = calcfld_file(fld+'z',f,**kwargs)
	else:
		print "square of %s is not supported"%fld
		return
	return (fx*fx+fy*fy+fz*fz)

def solve_ham_grad(fx,fy, dx=1, dy=1):
	'''
	Integrate dpsi = -fy dx + fx dy, psi[0,0]=0. Notes: (fx=dpsi/dy,fy=-dpsi/dx) is called hamiltonian gradient of psi, and	contours of psi give vector field (fx, fy)
	'''
	#FIXME: avoid loop for higher performance
	ny,nx=fx.shape
	psi=np.zeros((ny,nx))
	for jx in np.arange(1,nx):
		psi[0,jx]=psi[0,jx-1]-0.5*(fy[0,jx-1]+fy[0,jx])*dx
	for jy in np.arange(1,ny):
		psi[jy,:]=psi[jy-1,:]+0.5*(fx[jy-1,:]+fx[jy,:])*dy
	return psi

def calc_psi(f,vxname='hx',vyname='hz',dx=1,dy=1,**kwargs):
	'''
	Calculate magnetic flux, whose contour gives magnetic field lines.
	'''
	fx=getfld_file(vxname,f,**kwargs)
	fy=getfld_file(vyname,f,**kwargs)
	psi=solve_ham_grad(fx,fy,dx,dy)
	return psi

def calcfld_file(fld,f,**kwargs):
	if fld in ['pxxe','pyye','pzze','pxxi','pyyi','pzzi','pxye','pxze','pyze','pxyi','pxzi','pyzi','pyxe','pzxe','pzye','pyxi','pzxi','pzyi']:
		dim1=fld[1]
		dim2=fld[2]
		kind=fld[3]
		field = calc_pressure(kind,dim1,dim2,f,**kwargs)
	elif fld in ['uxe','uye','uze','uxi','uyi','uzi']:
		dim=fld[1]
		kind=fld[2]
		field = calc_u(kind,dim,f,**kwargs)
	elif fld in ['ux','uy','uz']:
		dim=fld[1]
		field = calc_u('all',dim,f,**kwargs)
	elif fld in ['jxe','jye','jze','jxi','jyi','jzi']:# assume qe=-1, qi=+1
		dim=fld[1]
		kind=fld[2]
		if kind=='e':
			q=-1
		else:
			q=1
		field = q*getfld_file("v%c%c"%(kind,dim),f,**kwargs)
	elif fld in ['b2','e2']:
		field = calc_2(fld[0],f,**kwargs)
	elif fld == 'psi':
		field = calc_psi(f,**kwargs)
	else:
		print "%s is not supported!"%fld
		return
	return field

def calcfld(fld,ts=None,dtype=None,filename=None,**kwargs):
	'''
	Calculate physical quantity "fld" that is not stored directly in the PSC data files. "filename" is either given explicitly, or formatted from "ts"(timestep) and "dtype"('pfd'/'tfd').
	'**kwargs' will be passed to getfld_file, which includes: 
		'dims' specifies which two dimensions are wanted(one of 'xz', 'yz', and 'xy'; default is 'xz'). 
		For the first dimension(e.g., 'x' in 'xz'), indices ;i1b'/'i1e' specify a range from which data is retrieved. Data array will also be strided every 'stride1' step. If 'stride1' is not set, the value of 'stride' will be used, which is 1 (no striding) by default.	Similarly, for the second dimension, 'i2b', 'i2e', 'stride2' can be used.
	Supported quantities are: [pxxe,pyye,pzze,pxxi,pyyi,pzzi,pxye,pxze,pyze,pxyi,pxzi,pyzi,pyxe,pzxe,pzye,pyxi,pzxi,pzyi,uxe,uye,uze,uxi,uyi,uzi,ux,uy,uz,jxe,jye,jze,jxi,jyi,jzi]
	'''
	f = getfile(ts, dtype, filename)
	field = calcfld_file(fld,f,**kwargs)
	f.close()
	return field

'''
utilities
'''
def slicing(field, dim, n=None, r=None):
	'''Take a slice from a 2D array. 'dim'=1 or 2 for the first and second dimension. For example, 'x' is the first dimension of 'xz'. 'n' gives the index of the slice in the other dimension. If 'n' is not given, then it will be calculated by 'r' times the length of the other dimension.'''
	n2,n1 = field.shape
	if dim==1:
		if n==None:
			n=n2*r
		_field_ = field[n, :]
	if dim==2:
		if n==None:
			n=n1*r
		_field_ = field[:, n]
	return _field_

def smooth_2dbox(field,n=2):
	'''Do a basic smoothing of using a (2*n+1)x(2*n+1) box.'''
	from numpy import roll
	if n == 0:
		return field
	a = None
	for i in range(-n, n+1):
		for j in range(-n, n+1):
			t = roll(field, i, 0)
			t = roll(t,     j, 1)
			if a == None:
				a = t
			else:
				a = a + t
	return a/np.square(2*n+1)

'''
exporting data
'''
def dump2asc(fld,ts=None,dtype=None,filename=None,dims='xz',
	reshape=True,output_filename=None,mode='w',fmt='%10g'):
	field=getfld(fld,ts,dtype,filename,dims)
	if reshape:
		field=field.reshape(field.shape[1],field.shape[0])
	if output_filename==None:
		output_filename=raw_input('enter filename for output: ')
	np.savetxt(output_filename,field,fmt=fmt)

