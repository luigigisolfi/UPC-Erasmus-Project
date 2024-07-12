function [tfa,xf]= propTITF(ti,xi,tf,vfield,hmin,hmax,tol,iorb)
%--------------------------------------------------------------------------
% Propagates an inital condition xi at time ti up to time tf in
% the vectorfield vfield. hmin, hmax, tol are the parameters for rk45f. 
% iorb =0 or 1 determines the output. When iorb=0 output tfa ant xf are
% final time tf and its corresponding state, when iorb=1 the output is the
% full propagated trajectory tfa(i), xf(i,:) i=1... integration steps.
%
% NOTE: initial step h is set in the first line of the function.
%--------------------------------------------------------------------------
 h=10;
 if (abs(ti-tf) < 1.e-12), tfa=tf; xf=xi; return; end
 if iorb >0, tfa=[]; xf=[]; end
 if (tf < ti), h=-h; end
 t=ti;  x=xi;
 while ((t<tf && h>0) || (t>tf && h<0)) 
  if iorb > 0, tfa=[tfa; t]; xf=[xf; x]; end
  [t,x,h,err]=rk45f(t,x,h,hmin,hmax,tol,vfield);
 end
 while (abs(t-tf) > 1.e-12)
  h=tf-t;
  [t,x,h,err]=rk45f(t,x,h,hmin,hmax,tol,vfield);
 end
 if iorb <= 0, tfa=t; xf=x; else tfa=[tfa; t]; xf=[xf; x]; end
end
