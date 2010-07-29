# Will Fox, Mar 2010
#
# fields2mat.py
#
# Load PSC fields data into a matlab .mat file.  It is a mirror image
# of OUT_fields.f, we just match all the write statements therein.
#
# This script is designed to be edited as needed at the top,
# but after 'GO!' you shouldn't need to change.
#

from numpy import *
from struct import unpack
from scipy.io import savemat

out = {};

# whether to use time-averaged ('t') or point ('p') field values
PT = 't';

# datadir
datadir = '../data/'

outfile = datadir  +'/' + PT + 'fields.mat'
print 'Outfile: ' + outfile

# time slices to load
ts = arange(0, 20, 10)
out['ts'] = 1.0 * ts			# 1.0 * : force to be double

# total number of processors 
npe = 2;

# decimation of data
decX = 5
decZ = 4

# whether to decimate by averaging (1) or simple choosing (0)
# for averaging to work, we require that the number of data points 
# along z be a multiple of decZ, (and analogously for decX)
meanDec = 1

out['decX'] = array( [1.0*decX] )
out['decZ'] = array( [1.0*decZ] )

# data to end up with.  1 or 0 indicates whether you want to keep it
# Order reflects order of their output in OUT_fields.f, so is important
field = [
  {'name' : 'ne',   'keep' : 1},
  {'name' : 'ni',   'keep' : 1},
  {'name' : 'nn',   'keep' : 1},
  {'name' : 'ex',   'keep' : 1},
  {'name' : 'ey',   'keep' : 1},
  {'name' : 'ez',   'keep' : 1},
  {'name' : 'bx',   'keep' : 1},
  {'name' : 'by',   'keep' : 1},
  {'name' : 'bz',   'keep' : 1},
  {'name' : 'jx',   'keep' : 1},
  {'name' : 'jy',   'keep' : 1},
  {'name' : 'jz',   'keep' : 1},
  {'name' : 'jxex', 'keep' : 0},
  {'name' : 'jyey', 'keep' : 0},
  {'name' : 'jzez', 'keep' : 0},
  {'name' : 'poyx', 'keep' : 0},
  {'name' : 'poyy', 'keep' : 0},
  {'name' : 'poyz', 'keep' : 0},
  {'name' : 'ex2',  'keep' : 0},
  {'name' : 'ey2',  'keep' : 0},
  {'name' : 'ez2',  'keep' : 0},
  {'name' : 'bx2',  'keep' : 0},
  {'name' : 'by2',  'keep' : 0},
  {'name' : 'bz2',  'keep' : 0},
] 


# GO!

for k in range(len(ts)):
    print ts[k]
    for pe in range(npe):
        filename = '%s/%cfd_%06d_%07d.psc' % ( datadir, PT, pe, ts[k] )

        # no tfield at t = 0!
        if ts[k] == 0:
           filename = '%s/pfd_%06d_%07d.psc' % ( datadir, pe, ts[k] )

#        print filename
        f = open(filename);
        
#     ...quoting OUT_fields.f
#
#         write(11) dx, dy, dz, dt
#         write(11) i1mn, i1mx, i2mn, i2mx, i3mn, i3mx
#         write(11) r1n,r1x,r2n,r2x,r3n,r3x,shift_z
#         write(11) dowrite_ne ... 


        str = f.read(4);
        if (str != "PSC "):
             raise Exception('Expected to see PSC header on data file');

#       figure out endianness
        endian_test = f.read(4);
        LB = '<';
        if ( unpack( LB + 'i', endian_test)[0] != 1061962303 ):
           LB = '>';
           if ( unpack( LB + 'i', endian_test)[0] != 1061962303 ):
              raise Exception('Couldn''t figure out byte swapping!');
 
        version = unpack( LB + 'i', f.read(4) )[0];
        if (version > 1):
            raise Exception('I only work up to version 1!');

#       get to work 
   
        out['dx'] = array( unpack( LB + 'f', f.read(4) ) )
        out['dy'] = array( unpack( LB + 'f', f.read(4) ) )
        out['dz'] = array( unpack( LB + 'f', f.read(4) ) )
        out['dt'] = array( unpack( LB + 'f', f.read(4) ) )

        # Indices calculated on local proc
        i1mn = unpack( LB + 'i', f.read(4) )[0]
        i1mx = unpack( LB + 'i', f.read(4) )[0]
        i2mn = unpack( LB + 'i', f.read(4) )[0]
        i2mx = unpack( LB + 'i', f.read(4) )[0]
        i3mn = unpack( LB + 'i', f.read(4) )[0]
        i3mx = unpack( LB + 'i', f.read(4) )[0]
        
        # These are the globally saved indices
        r1n  = unpack( LB + 'i', f.read(4) )[0]
        r1x  = unpack( LB + 'i', f.read(4) )[0]
        r2n  = unpack( LB + 'i', f.read(4) )[0]
        r2x  = unpack( LB + 'i', f.read(4) )[0]
        r3n  = unpack( LB + 'i', f.read(4) )[0]
        r3x  = unpack( LB + 'i', f.read(4) )[0]

        # indices we will eventually grab from the whole domain, 
        # including decimation.  PYTHON style 0-offset
        r1_ind = arange(r1n, r1x, decX)
        r2_ind = arange(r2n, r2x, 1)
        r3_ind = arange(r3n, r3x, decZ)

        # which fields were written to the data file... PSC dowrite_* vars
        for j in range(24):
           field[j]['wrote'] = (unpack( LB + 'c', f.read(1) )[0] == '\x01')

        # Can now initialize out field matrices 
        if k == 0 and pe == 0:
            for j in range(24):
                if field[j]['wrote'] and field[j]['keep']:
                    out[field[j]['name']] = zeros( (len(ts), len(r1_ind), 
                                                   len(r2_ind), len(r3_ind)) )


        # Temporary field matrices (pre decimation)
        if pe == 0:
           tmp = {}
           for j in range(24):
               if field[j]['wrote'] and field[j]['keep']:
                   tmp[field[j]['name']] = zeros( (r1x-r1n, 
                                               r2x-r2n, r3x-r3n) )
     
        # Index range actually written into this data file
#        i1a = unpack( LB + 'i', f.read(4) )[0]
#        i1e = unpack( LB + 'i', f.read(4) )[0]
#        i2a = unpack( LB + 'i', f.read(4) )[0]
#        i2e = unpack( LB + 'i', f.read(4) )[0]
#        i3a = unpack( LB + 'i', f.read(4) )[0]
#        i3e = unpack( LB + 'i', f.read(4) )[0]

        i1a = i1mn;
        i1e = i1mx;
        i2a = i2mn;
        i2e = i2mx;
        i3a = i3mn;
        i3e = i3mx;

        # Make sure we got past the header
        if (f.read(4) != "DATA"):
             raise Exception('Expected to be past the header!');

        # how many data points in local field matrix
        Ndata = (i1e-i1a) * (i2e-i2a) * (i3e-i3a);

        if (i1e>=i1a) and (i2e>=i2a) and (i3e>=i3a): 
        
            # indices for jamming into temporary field matrix - PYTHON style
            t_i1a = i1a - r1n;
            t_i1e = i1e - r1n;
            t_i2a = i2a - r2n;
            t_i2e = i2e - r2n;
            t_i3a = i3a - r3n;
            t_i3e = i3e - r3n;      

            for j in range(24):
                if field[j]['wrote']:

                    datachunk = array(unpack( LB + '%df' % Ndata, f.read(4*Ndata)));


#                    print datachunk.shape
                    
                    if field[j]['keep']:
                        datachunk = reshape(datachunk, (i1e-i1a, 
                                              i2e-i2a, i3e-i3a), order='F')

                        tmp[field[j]['name']][t_i1a:t_i1e, 
                           t_i2a:t_i2e, t_i3a:t_i3e] = datachunk
                else:
                    if field[j]['keep']:
                        raise Exception('HELP! Want to keep field %s \
                               but it was not saved' % field[j]['name'] )
                # end if wrote
            # end for j fields 
        # end if Ndata > 0
                        
        f.close();
    # end for each proc

    # assemble out from tmp
    for j in range(24):
       if field[j]['wrote'] and field[j]['keep']:
           if meanDec==1:
               out[field[j]['name']][k,:,:,:] =  \
                 mean( mean( mean(
                  reshape(tmp[field[j]['name']], 
                     (decX, len(r1_ind), 1, len(r2_ind), decZ, len(r3_ind)),
                     order='F'),
                 4), 2), 0)
           else:
               out[field[j]['name']][k,:,:,:] =  \
                 tmp[field[j]['name']][r1_ind, r2_ind, r3_ind]
           # if meanDec
       # if wrote and keep
    # for each field

# end for each time slice

# Compression has to wait for Scipy 0.8
#savemat(outfile, out, do_compression=True)
savemat(outfile, out)


