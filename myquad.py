#QUAD   Numerically evaluate integral, adaptive Simpson quadrature.
#   Q = QUAD(FUN,A,B) tries to approximate the integral of function
#   FUN from A to B to within an error of 1.e-6 using recursive
#   adaptive Simpson quadrature.  The function Y = FUN(X) should
#   accept a vector argument X and return a vector result Y, the
#   integrand evaluated at each element of X.
#
#   Q = QUAD(FUN,A,B,TOL) uses an absolute error tolerance of TOL
#   instead of the default, which is 1.e-6.  Larger values of TOL
#   result in fewer function evaluations and faster computation,
#   but less accurate results.  The QUAD function in MATLAB 5.3 used
#   a less reliable algorithm and a default tolerance of 1.e-3.
#
#   [Q,FCNT] = QUAD(...) returns the number of function evaluations.
#
#   QUAD(FUN,A,B,TOL,TRACE) with non-zero TRACE shows the values
#   of [fcnt a b-a Q] during the recursion.
#
#   QUAD(FUN,A,B,TOL,TRACE,P1,P2,...) provides for additional
#   arguments P1, P2, ... to be passed directly to function FUN,
#   FUN(X,P1,P2,...).  Pass empty matrices for TOL or TRACE to
#   use the default values.
#
#   Use array operators .*, ./ and .^ in the definition of FUN
#   so that it can be evaluated with a vector argument.
#
#   Function QUADL may be more efficient with high accuracies
#   and smooth integrands.
#
#   Example:
#       FUN can be specified three different ways.
#
#       A string expression involving a single variable:
#          Q = quad('1./(x.^3-2*x-5)',0,2);
#
#       An inline object:
#          F = inline('1./(x.^3-2*x-5)');
#          Q = quad(F,0,2);
#
#       A function handle:
#          Q = quad(@myfun,0,2);
#          where myfun.m is an M-file:
#             function y = myfun(x)
#             y = 1./(x.^3-2*x-5);
#
#   See also QUADL, DBLQUAD, INLINE, @.

#   Based on "adaptsim" by Walter Gander.
#   Ref: W. Gander and W. Gautschi, "Adaptive Quadrature Revisited", 1998.
#   http://www.inf.ethz.ch/personal/gander
#   Copyright 1984-2001 The MathWorks, Inc.
#   $Revision: 1.1.1.1 $  $Date: 2002/12/05 20:52:54 $

import numpy as np


def myquad(func, a, b, tol=1e-6, trace=0, *args ):

    if np.size(tol) == 0:
        tol = 1e-6
    if np.size(trace) == 0:
        trace = 0

    # Initialize with three unequal subintervals.
    h = 0.13579*(b-a)
    x = [a, a+h, a+2*h, (a+b)/2, b-2*h, b-h, b]
    argv = list(args)
    y0 = func(x[0], argv[0], argv[1], argv[2], argv[3], argv[4])
    y1 = func(x[1], argv[0], argv[1], argv[2], argv[3], argv[4])
    y2 = func(x[2], argv[0], argv[1], argv[2], argv[3], argv[4])
    y3 = func(x[3], argv[0], argv[1], argv[2], argv[3], argv[4])
    y4 = func(x[4], argv[0], argv[1], argv[2], argv[3], argv[4])
    y5 = func(x[5], argv[0], argv[1], argv[2], argv[3], argv[4])
    y6 = func(x[6], argv[0], argv[1], argv[2], argv[3], argv[4])
    # y1 = feval(f, x(1), varargin{:})
    # y2 = feval(f, x(2), varargin{:})
    # y3 = feval(f, x(3), varargin{:})
    # y4 = feval(f, x(4), varargin{:})
    # y5 = feval(f, x(5), varargin{:})
    # y6 = feval(f, x(6), varargin{:})
    # y7 = feval(f, x(7), varargin{:})
    fcnt = 7

    # Fudge endpoints to avoid infinities.
    y0fin = np.full(len(y0), 2, dtype=int)
    np.isfinite(y0, y0fin)
    if not np.min(np.min(np.min(y0fin))):
        #y1 = feval(f,a+eps*(b-a),varargin{:});
        y0 = func(a+np.spacing(1)*(b-a), argv[0], argv[1], argv[2], argv[3], argv[4])
        fcnt = fcnt+1

    y6fin = np.full(len(y6), 2, dtype=int)
    np.isfinite(y6,y6fin)
    if not np.min(np.min(np.min(y6fin))):
        #y7 = feval(f,b-eps*(b-a),varargin{:});
        y6 = func(b - np.spacing(1)*(b-a), argv[0], argv[1], argv[2], argv[3], argv[4])
        fcnt = fcnt+1

    if len(np.shape(y0))==1:
        size0 = np.size(y0,0)
        size1 = 1
        size2 = 1
    elif len(np.shape(y0))==2:
        size0 = np.size(y0,0)
        size1 = np.size(y0,1)
        size2 = 1
    elif len(np.shape(y0)) == 3:
        size0 = np.size(y0,0)
        size1 = np.size(y0,1)
        size2 = np.size(y0,2)

    Q = np.zeros((size1, size2, 3, size0))
    # Call the recursive core integrator.
    hmin = np.spacing(1)/1024.0*abs(b-a)
    warn = np.zeros(3)
    [Q[:, :, 0, :],fcnt,warn[0]] = quadstep(func, x[0], x[2], np.squeeze(y0), np.squeeze(y1), np.squeeze(y2), tol, trace, fcnt, hmin, argv[0], argv[1], argv[2], argv[3], argv[4])
    [Q[:, :, 1, :],fcnt,warn[1]] = quadstep(func, x[2], x[4], np.squeeze(y2), np.squeeze(y3), np.squeeze(y4), tol, trace, fcnt, hmin, argv[0], argv[1], argv[2], argv[3], argv[4])
    [Q[:, :, 2, :],fcnt,warn[2]] = quadstep(func, x[4], x[6], np.squeeze(y4), np.squeeze(y5), np.squeeze(y6), tol, trace, fcnt, hmin, argv[0], argv[1], argv[2], argv[3], argv[4])
    Q = np.squeeze(np.sum(Q,2))
    warni = np.max(warn)

    if warni == 1:
        print 'warning - Minimum step size reached; singularity possible.'
    elif warni == 2:
        print 'warning - Maximum function count exceeded; singularity likely.'
    elif warni == 3:
        print 'warning - Infinite or Not-a-Number function value encountered.'

    return [Q, fcnt]

def quadstep(func, a, b, fa, fc, fb, tol, trace, fcnt, hmin, *args):
    #QUADSTEP  Recursive core routine for function QUAD.

    maxfcnt = 10000

    # Evaluate integrand twice in interior of subinterval [a,b].
    h = b - a
    c = (a + b)/2.0
    if (abs(h) < hmin) or c == a or c == b:
       # Minimum step size reached; singularity possible.
       Q = h*fc
       warn = 1
       return [Q, fcnt, warn]

    argv = list(args)
    x = [(a + c)/2.0, (c + b)/2.0]
    y0 = func(x[0], argv[0], argv[1], argv[2], argv[3], argv[4])
    y1 = func(x[1], argv[0], argv[1], argv[2], argv[3], argv[4])
    #y1 = feval(f, x(1), varargin{:})
    #y2 = feval(f, x(2), varargin{:})

    fcnt = fcnt + 2
    if fcnt > maxfcnt:
       # Maximum function count exceeded; singularity likely.
       Q = h*fc
       warn = 2
       return [Q, fcnt, warn]
    fd = y0
    fe = y1

    # Three point Simpson's rule.
    Q1 = (h/6.0)*(fa + 4*fc + fb)

    # Five point double Simpson's rule.
    Q2 = (h/12.0)*(fa + 4.0*fd + 2.0*fc + 4.0*fe + fb)

    # One step of Romberg extrapolation.
    Q = Q2 + (Q2 - Q1)/15.0

    Qfin = np.full(len(Q), 2, dtype=int)
    np.isfinite(Q, Qfin)
    if not np.min(np.min(np.min(Qfin))):
        # Infinite or Not-a-Number function value encountered.
        warn = 3
        return [Q, fcnt, warn]

    if trace:
        print '%8.0f %16.10f %18.8e %16.10f' % (fcnt,a,h,Q)

    # Check accuracy of integral over this subinterval.

    if np.max(np.max((abs(Q2 - Q)))) <= tol:
        warn = 0
        return [Q, fcnt, warn]

    # Subdivide into two subintervals.
    else:
        [Qac, fcnt, warnac] = quadstep(func, a, c, fa, fd, fc, tol, trace, fcnt, hmin, argv[0], argv[1], argv[2], argv[3], argv[4])
        [Qcb, fcnt, warncb] = quadstep(func, c, b, fc, fe, fb, tol, trace, fcnt, hmin, argv[0], argv[1], argv[2], argv[3], argv[4])
        Q = Qac + Qcb
        warn = max(warnac, warncb)
        return [Q, fcnt, warn]




