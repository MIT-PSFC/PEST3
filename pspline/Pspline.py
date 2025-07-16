#!/usr/bin/env python

"""
$Id: Pspline.py 1045 2011-11-08 22:19:18Z pletzer $
---------------------------------------------------------------------------
 This code was developed at Tech-X (www.txcorp.com). It is free for any one
 to use but comes with no warranty whatsoever. Use at your own risk. 
 Thanks for reporting bugs to pletzer@txcorp.com. 
---------------------------------------------------------------------------
"""

import os
import numpy
from ctypes import CDLL, POINTER, c_int, c_double, c_float, byref, \
     c_int, RTLD_GLOBAL, RTLD_LOCAL


class Pspline:
    """
    1D, 2D, 3D spline, Akima, or linear interpolation
    """

    def __init__(self, sharedObject, numUnderscores=1):
        """
        Constructor
        @param sharedObject path to libpsplineAll.{so,dylib,dll} shared object
        @param numUnderscores number of trailing underscores in symbol table
        """
        self.isReady = False
        self.isAllocated = False

        # opaque type holder
        self.holder = numpy.array([0]*12)
        self.handle = self.holder.ctypes.data_as(POINTER(c_int))
        self.object = CDLL(sharedObject, mode=RTLD_GLOBAL)
        self.underscore = '_' * numUnderscores
        self.ier = c_int()
        self.bcTypes = []
        self.bcVals  = []
        self.numDims = 0  # 1, 2, or 3
        self.type    = '' # r4 or r8
        self.cType   = ''
        self.nType   = ''
        self.interpMethod = 'linear'

    def __del__(self):
        """
        Destructor
        """
        if self.isAllocated:
            meth = 'self.object.czspline_free%d_%s%s' % \
                (self.numDims, self.type, self.underscore)
            args = '(self.handle, byref(self.ier))'
            eval(meth + args)
            assert(self.ier.value == 0)
        

    def setup(self, axes, functValues, boundTypes=[], boundValues=[]):
        """
        Allocate object and compute spline coefficients
        @param axes list of ndim axes [x1, x2, ...]
        @param functValues array of function values, *varies fastest along x1*
        @param boundTypes list of ndim tuples [(bc0, bc1), ...]
                          (-1 for periodic, 0 for not-a-knot ,...)
        @param boundValues list of ndim tuples [(bv0, bv1), ...]
                           only used if bc0, ... are > 0
        
        """
        self.numDims = len(axes)
        assert(len(functValues.shape) == self.numDims)

        if functValues.dtype == numpy.float32:
            self.type = 'r4'
            self.cType = 'c_float'
            self.nType = numpy.float32
        elif functValues.dtype == numpy.float64:
            self.type = 'r8'
            self.cType = 'c_double'
            self.nType = numpy.float64
        else:
            raise TypeError, 'functValues must be float32/float64, got "%s"'% \
                str(functValues.dtype)
        
        if boundTypes != []:
            assert(len(boundTypes) == self.numDims)
            self.bcTypes = boundTypes
        else:
            # defaults to not-a-knot
            self.bcTypes = numpy.array( [ (0,0,) for idim in range(self.numDims) ] )

        if boundValues:
            assert(len(boundValues) == self.numDims)
            self.bcVals = boundValues
        for idim in range(self.numDims):
            # order of axes is x3, x2, x1 in functValues
            n = functValues.shape[self.numDims-idim-1]
            assert(n == len(axes[idim]))

        # initialize
        if self.interpMethod.lower() == 'linear':
            meth = 'self.object.czlinear_init%d_%s%s' % \
                (self.numDims, self.type, self.underscore)
            args = '(self.handle'
            for idim in range(self.numDims):
                # order of axes is x3, x2, x1 in functValues
                n = functValues.shape[self.numDims-idim-1]
                args += ', byref(c_int(%d))' % n
            args += ', byref(self.ier))'
            eval(meth + args)
            assert(self.ier.value == 0)
        else:
            # spline or Akima
            meth = 'self.object.czspline_init%d_%s%s' % \
                (self.numDims, self.type, self.underscore)
            args = '(self.handle'
            for idim in range(self.numDims):
                # order of axes is x3, x2, x1 in functValues
                n = functValues.shape[self.numDims-idim-1]
                args += ', byref(c_int(%d))' % n
            for idim in range(self.numDims):
                args += ',  self.bcTypes[%d].ctypes.data_as(POINTER(c_int))' \
                    % idim
            args += ', byref(self.ier))'
            eval(meth + args)
            assert(self.ier.value == 0)

        self.isAllocated = True

        # set axes
        meth = 'self.object.czspline_set_axes%d_%s%s' % \
            (self.numDims, self.type, self.underscore)
        args = '(self.handle'
        for idim in range(self.numDims):
            # order of axes is x3, x2, x1 in functValues
            n = functValues.shape[self.numDims-idim-1]
            args += ', byref(c_int(%d))' % n
        for idim in range(self.numDims):
            args += ', axes[%s].ctypes.data_as(POINTER(%s))' \
                % (idim, self.cType)
        args += ', byref(self.ier))'
        eval(meth + args)
        assert(self.ier.value == 0)
        
        # setup
        meth = 'self.object.czspline_setup%d_%s%s' % \
            (self.numDims, self.type, self.underscore)
        args = '(self.handle'
        for idim in range(self.numDims):
            # order of axes is x3, x2, x1 in functValues
            n = functValues.shape[self.numDims-idim-1]
            args += ', byref(c_int(%d))' % n
        args += ', functValues.ctypes.data_as(POINTER(%s)), byref(self.ier))' \
            % self.cType
        eval(meth + args)
        assert(self.ier.value == 0)

        self.isReady = True

    def getPointInterp(self, pt):
        """
        Interpolation at a single point
        @param pt coordinate point (x1, x2, ...)
        @return interpolated value
        """
        assert(self.isReady)

        fIntrp = eval(self.cType + '()')
        meth = 'self.object.czspline_interp%d_%s%s' % \
            (self.numDims, self.type, self.underscore)
        args = '(self.handle'
        for idim in range(self.numDims):
            args += ', byref(%s(pt[%d]))' % (self.cType, idim)
        args += ', byref(fIntrp), byref(self.ier))'
        eval(meth + args)
        assert(self.ier.value == 0)

        return fIntrp.value

    def getCloudInterp(self, pts):
        """
        Interpolation at arbitrary points
        @param pts list of [x1, x2, ...] points with len(x1)==len(x2)==..
        @return interpolated values
        """
        sz = len(pts[0])
        for idim in range(self.numDims):
            assert(len(pts[idim]) == sz)

        fIntrp = numpy.zeros( (sz,), self.nType)
        meth = 'self.object.czspline_interp%d_cloud_%s%s' % \
            (self.numDims, self.type, self.underscore)
        csz = c_int(sz)
        args = '(self.handle, byref(csz)'
        for idim in range(self.numDims):
            args += ', pts[%d].ctypes.data_as(POINTER(%s))' \
                % (idim, self.cType)
        args += ', fIntrp.ctypes.data_as(POINTER(%s)), byref(self.ier))' \
            % self.cType
        eval(meth + args)
        assert(self.ier.value == 0)

        return fIntrp

    def getArrayInterp(self, pts):
        """
        Interpolation on a structured grid
        @param pts list of [x1, x2, ...] points
        @return interpolated values x1 cross x2 cross ...
        """
        ndims = len(pts)
        assert(ndims == self.numDims)

        dims = '('
        # dimension order is reversed
        for idim in range(ndims-1, -1, -1):
            dims += '%d,' % len(pts[idim])
        dims += ')'
        fIntrp = eval('numpy.zeros(%s, self.nType)' % dims)
        meth = 'self.object.czspline_interp%d_array_%s%s' % \
            (self.numDims, self.type, self.underscore)
        args = '(self.handle'
        for idim in range(self.numDims):
            args += ', byref(c_int(%d))' % len(pts[idim])
        for idim in range(self.numDims):
            args += ', pts[%d].ctypes.data_as(POINTER(%s))' \
                % (idim, self.cType)
        args += ', fIntrp.ctypes.data_as(POINTER(%s)), byref(self.ier))' \
            % self.cType
        eval(meth + args)
        assert(self.ier.value == 0)

        return fIntrp

        
##############################################################

def f3(z): return z**3
def f2(y): return y**3
def f1(x): return numpy.cos(2*numpy.pi*(x-0.2443))

def test2(sharedObject):

    import copy

    def getDataFromString(strng, dtype = numpy.float64):
        """
        Read data from string and return two arrays
        @param strng a string containing 2 space separated columns, 
                     the x and f values
        @param dtype a valid numpy type
        @return x, f numpy arrays
        """
        data = numpy.fromstring(strng, 
                                dtype = dtype, 
                                count = -1, 
                                sep = ' ')
        nrows = len(data) // 2
        data = data.reshape( (nrows, 2) )
        # copy is required, otherwise the returned arrays 
        # are not contiguous
        return copy.copy(data[:,0]), copy.copy(data[:,1])
        

    inputStr = """
0 0.990174
0.0454545 0.99854
0.0909091 0.987857
0.136364 0.982702
0.181818 0.96135
0.227273 0.938445
0.272727 0.91435
0.318182 0.886148
0.363636 0.854768
0.409091 0.818618
0.454545 0.779092
0.5 0.73475
0.545455 0.687039
0.590909 0.63454
0.636364 0.578721
0.681818 0.518193
0.727273 0.454462
0.772727 0.3862
0.818182 0.315012
0.863636 0.239716
0.909091 0.16214
0.954545 0.0814472
1 0
"""

    outputStr = """
0 0.990174
0.030303 0.999836
0.0606061 0.995128
0.0909091 0.987857
0.121212 0.985305
0.151515 0.977205
0.181818 0.96135
0.212121 0.945909
0.242424 0.93079
0.272727 0.91435
0.30303 0.895923
0.333333 0.876112
0.363636 0.854768
0.393939 0.831094
0.424242 0.805872
0.454545 0.779092
0.484848 0.749959
0.515152 0.719274
0.545455 0.687039
0.575758 0.652462
0.606061 0.616356
0.636364 0.578721
0.666667 0.538772
0.69697 0.497378
0.727273 0.454462
0.757576 0.409229
0.787879 0.363096
0.818182 0.315012
0.848485 0.264051
0.878788 0.21645
0.909091 0.16214
0.939394 0.100738
0.969697 0.0702882
1 0.0702882
"""
    xs, fs = getDataFromString(inputStr)
    xi, fi = getDataFromString(outputStr)
    print xs
    # interpolate
    spl = Pspline(sharedObject)
    notAKnotBC = numpy.zeros( (2,), numpy.int )
    spl.setup( axes = [xs,], functValues = fs, boundTypes = [notAKnotBC,] )
    fiNew = spl.getArrayInterp([xi,])
    print fiNew
    

def test(sharedObject):

    MAXERROR = 1.e-8

    nx1, nx2, nx3 = 11, 21, 31
    x1min, x1max = 0., 1.
    x2min, x2max = 0., 2.
    x3min, x3max = 0., 3.
    dx1 = (x1max - x1min) / float(nx1-1)
    dx2 = (x2max - x2min) / float(nx2-1)
    dx3 = (x3max - x3min) / float(nx3-1)
    x1 = x1min + dx1 * numpy.arange(0, nx1)
    x2 = x2min + dx2 * numpy.arange(0, nx2)
    x3 = x3min + dx3 * numpy.arange(0, nx3)
    f = numpy.zeros( (nx3,), numpy.float64 )
    ff = numpy.zeros( (nx3, nx2,), numpy.float64 )
    fff = numpy.zeros( (nx3, nx2, nx1,), numpy.float64 )
    for k in range(nx3):
        z = x3[k]
        f[k] = f3(z)
        for j in range(nx2):
            y = x2[j]
            ff[k, j] = f2(y) * f3(z)
            for i in range(nx1):
                x = x1[i]
                fff[k, j, i] = f1(x) * f2(y) * f3(z)


    bc1 = numpy.array([-1]*2) # periodic
    bc2 = numpy.array([0]*2)  # not-a-knot
    bc3 = numpy.array([0]*2)  # not-a-knot

    #
    # 1d test
    #
    ps1 = Pspline(sharedObject)
    ps1.setup(axes=[x3,], functValues=f, boundTypes=[bc3,])

    z = 1.5
    print ps1.getPointInterp( (z,) )
    assert( abs(ps1.getPointInterp( (z,) ) - f3(z)) < MAXERROR )

    zs = numpy.array([0.5, 1.0, 2.0], numpy.float64)
    print ps1.getCloudInterp([zs,])
    assert( abs(ps1.getArrayInterp([zs,]) - f3(zs)).any() < MAXERROR )

    # should be the same as cloud interpolation in 1D
    print ps1.getArrayInterp([zs,])
    assert( abs(ps1.getArrayInterp([zs,]) - f3(zs)).any() < MAXERROR )

    #
    # 2d test
    #
    ps2 = Pspline(sharedObject)
    ps2.setup(axes=[x2, x3,], functValues=ff, boundTypes=[bc2, bc3,])

    y, z = 1.0, 1.5
    print ps2.getPointInterp( (y, z,) )
    assert( abs(ps2.getPointInterp( (y, z,) ) - f2(y)*f3(z)) < MAXERROR )

    ys = numpy.array([0.3, 0.7, 1.5], numpy.float64)
    zs = numpy.array([0.5, 1.0, 2.0], numpy.float64)
    print ps2.getCloudInterp([ys, zs,])
    assert( abs(ps2.getCloudInterp([ys, zs,]) - f2(ys)*f3(zs)).all() < MAXERROR )

    print ps2.getArrayInterp([ys, zs,])
    print numpy.outer(f3(zs), f2(ys))
    assert( abs(ps2.getArrayInterp([ys, zs,]) - numpy.outer(f3(zs), f2(ys))).all() < MAXERROR )
   
    #
    # 3d test
    #
    ps3 = Pspline(sharedObject)
    ps3.setup(axes=[x1, x2, x3,], functValues=fff, boundTypes=[bc1, bc2, bc3,])

    x, y, z = 0.5, 1.0, 1.5
    print ps3.getPointInterp( (x, y, z,) )
    assert( abs(ps3.getPointInterp( (x, y, z,) ) - f1(x)*f2(y)*f3(z)) < MAXERROR )

    xs = numpy.array([0.1, 0.4, 0.8], numpy.float64)
    ys = numpy.array([0.3, 0.7, 1.5], numpy.float64)
    zs = numpy.array([0.5, 1.0, 2.0], numpy.float64)
    print ps3.getCloudInterp([xs, ys, zs,])
    assert( abs(ps3.getCloudInterp([xs, ys, zs,]) - f1(xs)*f2(ys)*f3(zs)).all() < MAXERROR )

    print ps3.getArrayInterp([xs, ys, zs,])
    ff3 = numpy.multiply.outer( f3(zs), numpy.ones((len(ys), len(xs),), numpy.float64) )
    ff2 = numpy.multiply.outer( numpy.outer( numpy.ones((len(zs),), numpy.float64), f2(ys) ), \
                                    numpy.ones((len(xs),), numpy.float64) )
    ff1 = numpy.multiply.outer( numpy.ones((len(zs), len(ys),), numpy.float64), f1(xs) )
    assert( abs(ps3.getArrayInterp([xs, ys, zs,]) - ff3*ff2*ff1).all() < MAXERROR )
   

if __name__ == '__main__': 
    import sys
    import optparse
    defaultPsplineLib = os.environ['HOME'] + \
        '/software/pspline/lib/libpsplineAll.so'
    parser = optparse.OptionParser()
    parser.add_option('-l', '--lib', action='store', type='string',
                      dest='sharedObject', help='pspline shared library path',
                      default=defaultPsplineLib)
    options, args = parser.parse_args()
    #test(options.sharedObject)
    test2(options.sharedObject)
    
