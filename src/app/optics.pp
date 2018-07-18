#!/usr/bin/python3
#####################

#include "language.h"
#include "math.h"
#include "io_util.h"
#include "init.h"
#include "test.h"
#include "spacetimes.h"
#include "runge_kutta.h"
#include "precision.h"

import runge_kutta,fancy,angular,vector,keplerian,transform,schwarzschild
#if "LANG" eq "python"
import sys,os
#endif

verbosity=1

csv=TRUE
csv_file = 'a.csv'

#--------------------------------------------------------------------------------------------------

def optics(r,tol):
  """
  Determine a table of optical aberration angles for an observer in the Schwarzschild spacetime.
  """
  spacetime = SP_SCH
  chart = CH_SCH
  pars = {}
  aa = 1-1/r
  x_obs = [0.0,r,1.0,0.0,0.0]
  # Calculate the velocity vector of the observer, for the standard state of motion:
  v_obs = [1/aa,-sqrt(1-aa),0.0,0.0,0.0]
  #    ... solution of the equations E=Atdot=1, |v|^2=Atdot^2-(1/A)rdot^2=1
  #        We check below that it actually is a solution.
  energy_obs = aa*v_obs[0]
  if abs(energy_obs-1.0)>EPS*10:
    THROW("energy of observer="+str(energy_obs))
  if abs(vector.norm(spacetime,chart,pars,x_obs,v_obs)-1.0)>EPS*10:
    THROW("observer's velocity has bad norm")
  # Calculate the observer's radial vector rho, parallel-transported in from infinity.
  rho = vector.proj(spacetime,chart,pars,x_obs,v_obs,[0.0,1.0,0.0,0.0,0.0])
  # ... take rhat and project out the part orthogonal to the observer's velocity; not yet normalized
  rho = vector.normalize_spacelike(spacetime,chart,pars,x_obs,rho)
  # ... now normalize it to have |rho|^2=-1
  if abs(vector.norm(spacetime,chart,pars,x_obs,rho)+1.0)>EPS*10:
    THROW("norm(rho)!=-1")
  if r>1.0:
    max_le = r/sqrt(aa) # maximum possible L/E for a photon at this r
  else:
    max_le = 10.0
  done = FALSE
  n_angles = 10
  for in_n_out in range(2): # 0=outward, 1=inward
    for i in range(n_angles):
      z = (float(i)/float(n_angles))
      if in_n_out==1:
        z = 1.0-z
      le = z*max_le # photon's ratio L/E of angular momentum to energy
      # Initial position:
      x = x_obs
      if le==0.0:
        v = [1.0,aa,0.0,0.0,0.0]
      else:
        z = r*r/(le*le)-aa
        if z<0.0:
          break # a photon with angular momentum this big can't exist at this r
        dphi_dr = (1.0/r)/sqrt(z)
        if in_n_out==1:
          dphi_dr = -dphi_dr
        dphi_dt = le*aa/(r*r)
        v = [1.0,dphi_dt/dphi_dr,0.0,dphi_dt,0.0] # tangent vector
      norm = vector.norm(spacetime,chart,pars,x,v)
      if abs(norm)>EPS*10:
        THROW("norm="+str(norm))
      v_perp = vector.proj(spacetime,chart,pars,x_obs,v_obs,v)
      # ... part of photon's velocity orthogonal to observer's velocity
      v_perp = vector.normalize_spacelike(spacetime,chart,pars,x_obs,v_perp)
      # ... normalized
      alpha = acos(-vector.inner_product(spacetime,chart,pars,x_obs,rho,v_perp))
      # ... angle at which the observer says the photon is emitted, see docs
      n = 100
      ndebug=0
      if verbosity>=3:
        ndebug=n/10
      ri = r
      lambda_max = 10.0 # fixme, sort of random
      while True:
        opt = {'lambda_max':lambda_max,'ndebug':ndebug,'sigma':1,
                  'future_oriented':FALSE,'tol':tol}
        err,final_x,final_v,final_a,final_lambda,info,sigma  = \
                fancy.trajectory_schwarzschild(spacetime,chart,pars,x,v,opt)
        if err!=0:
          if err==RK_INCOMPLETE:
            done = TRUE
            break
          THROW("error: "+str(err))
        rf = final_x[1]
        if IS_NAN(rf):
          THROW("rf is NaN")
        if rf>1.0e9 and rf>100.0*r:
          break
        if rf>10.0 and final_v[1]>0 and rf>ri: # quick and dirty test, far away and getting farther
          f = (rf-ri)/ri # fractional change in r from the iteration we just completed
          if f<0.1:
            lambda_max = lambda_max*(0.1/f) # try to get at least a 10% change in r with each iteration
        if verbosity>=2:
          PRINT("x=",final_x)
        x = final_x
        v = final_v
        ri = rf
      if done:
        break
      beta = atan2(final_x[3],final_x[2])
      if verbosity>=1:
        PRINT("r=",r,", alpha=",alpha*180.0/math.pi," deg, beta=",beta*180.0/math.pi," deg")
        #PRINT("final_x=",final_x)
        #PRINT("final_v=",final_v)
    if done:
      break

def main():
  optics(2.0,1.0e-6)
  print("done")

main()
