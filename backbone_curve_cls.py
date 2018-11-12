######################
## Version 0.1 #######
## /**********************************************************************
## This is a wrapper, a test environment to the SBAT(Smallest ball for algebraic
## topological analysis)   
## ***********************************************************************/
######################
import sys, time, os, fnmatch, pickle
import numpy as np
import scipy as sp
import scipy.interpolate as spinterp


## ################################################
## ################################################
class backbone_curve_cls:
  """
  Task: wrapper arroun the scipy interpolate for multidimensional splines
        + computes the curvature and torsion of th estimated curve

  """
  ## ------------------------------------------
  def __init__(self,X3d=None):
    """
    Input:  X3d   2d array of 3d points of a sample of the curve 
    """

    self.X3d=X3d

    return
  
  ## ------------------------------------------------
  def splprep(self,X3d):
    """
    Task: it is a wrapper around the scipy.splprep
          it computes the spline polynomials

    Input:   X3d   2d array of 3d points of a sample of the curve

    Output:
             tck   the output of splprep:
                   tck[0]= array of knots
                   tck[1]= list of arrays
                           =the coefficients of the spline polynomials
                           in each point
                   tck[2]= degree of the common polynomials
             u     array of the curve parameter values
    """

    m,n=X3d.shape
    lx=[ X3d[:,_] for _ in range(n)]
    tck,u=spinterp.splprep(lx, k=5)

    return(tck,u)
  ## ------------------------------------------------
  def splade(self,tx,tck):
    """
    Task: wrapper around the scipy splade
          it computes the derivatives of the curves 
    Input:   tx    1d array of parameter values of the points interpolated
                   0 <= tx[i] <=1, and ordered increasingly
                   e.g. np.linspace(0,1,m)
             tck   the output of splprep:
                   tck[0]= array of knots
                   tck[1]= list of arrays
                           =the coefficients of the spline polynomials
                           in each point
                   tck[2]= degree of the common polynomials
    Output: lder   the values of the derivatives
                   list of list of arrays
                   D^k x_i[j]=lder[j][i][k]
                       k degree of derivative
                       i index of a point
                       j component of point vector 
    """

    lder=spinterp.spalde(tx,tck)

    return(lder)

  ## ------------------------------------------------
  def curvature(self,xlder):
    """
    Task:   computes the curvature of a curve in given points

    Input:  xlder  3d array of point wise derivatives
                   D^k x_i[j]=xlder[j,i,k]
                       k degree of derivative
                       i index of a point
                       j component of point vector

    Output: xcurvature   the value of the curvature in each points
                         given by lder 

    """

    n=len(xlder)
    m=len(xlder[0])
    kD=len(xlder[0][0])

    xcurvature=np.zeros(m)
    ## derivatives
    xgamma_1=np.copy(xlder[:,:,1].T)
    xgamma_2=np.copy(xlder[:,:,2].T)
    
    ## scipy interpolate uses parameter range [0,1], thus independent on
    ## the length of the pointset, to apply common speed on different
    ## sequence renormalization is needed
    xgamma_1/=m
    xgamma_2/=m**2

    ## inner product between derivatives
    xgamma_12=np.sum(xgamma_1*xgamma_2,1)

    ## curvature
    xdet=np.sum(xgamma_1**2,1)*np.sum(xgamma_2**2,1)-xgamma_12
    ix=np.where(xdet<0)[0]
    xdet[ix]=0

    xcurvature=np.sqrt(xdet) \
      /np.sum(xgamma_1**2,1)**1.5
      
    return(xcurvature)

  ## ------------------------------------------------
  def torsion(self,xlder):
    """
    Task:   computes the torsion of a curve in given points

    Input:  xlder  3d array of point wise derivatives
                   D^k x_i[j]=xlder[j][i][k]
                       k degree of derivative
                       i index of a point
                       j component of point vector

    Output: xcurvature   the value of the torsion in each points
                         given by xlder 

    """

    n=len(xlder)
    m=len(xlder[0])
    kD=len(xlder[0][0])

    xtorsion=np.zeros(m)
    ## derivatives
    xgamma_1=np.copy(xlder[:,:,1].T)
    xgamma_2=np.copy(xlder[:,:,2].T)
    xgamma_3=np.copy(xlder[:,:,3].T)

    ## scipy interpolate uses parameter range [0,1], thus independent on
    ## the length of the pointset, to apply common speed on different
    ## sequence renormalization is needed
    xgamma_1/=m
    xgamma_2/=m**2
    xgamma_3/=m**3

    ## the determinant at each point of the derivatives
    xnumerator=np.zeros(m)
    xdenominator=np.zeros(m)
    for i in range(m):
      xmat=np.vstack((xgamma_1[i],xgamma_2[i],xgamma_3[i]))
      xnumerator[i]=np.linalg.det(xmat) 
      xnormal=np.cross(xgamma_1[i],xgamma_2[i])
      xdenominator[i]=np.sum(xnormal**2)

    ix0=np.where(xdenominator==0)[0]
    xdenominator[ix0]=1
    xnumerator[ix0]=0
    
    xtorsion=xnumerator/xdenominator
    
    return(xtorsion)
    
## ################################################################
## if __name__ == "__main__":
##   if len(sys.argv)==1:
##     iworkmode=0
##   elif len(sys.argv)>=2:
##     iworkmode=eval(sys.argv[1])
##   ball_main(iworkmode)
