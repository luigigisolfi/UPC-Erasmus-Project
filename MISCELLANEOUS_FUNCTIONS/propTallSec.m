function [tf, xf] = propTallSec(ti,xi,h,hmin,hmax,tol,vfield,gdgsec,ntall,iorb)
% ------------------------------------------------------------------------------
% Propagates the inital state xi at time ti with stopping condition g(x)=0
% after ntall cuts.
% NOTE: usually ntall > 0, but also one can give ntall=0 and the Newton
%       procedure starts from ti, xi to find tf, xf with g(xf)=0. 
%
% Function g(x) must be provided with its row gradient vector in gdgsec. 
% This is, [g,dg]=gdsec(x).
%
% h,hmin,hmax,tol,vfield: usual rk45f input parameters. In particular the
% sign of h determines the sense of propagation and vfield the vectorfield.
%
% iorb =-1, 0 or 1 determines the output. When iorb=0 output tf ant xf are
% final time and state with g(xf)=0, when iorb>0 the output is the full 
% propagated trajectory tf(i), xf(i,:) i=1... integration steps, when
% iorb<0 is the same as iorb>0 plus information about the procedure is given.
% ------------------------------------------------------------------------------
 t=ti; x=xi;
 if iorb ~= 0, tf=ti; xf=xi; end
 if iorb < 0
   [gva,~]=gdgsec(x);
   fprintf('propTallSec. Initial time and value of g: %f %e\n',t,gva); 
 end
 if ntall > 0
  [t,x,h,~]=rk45f(t,x,h,hmin,hmax,tol,vfield);
  [gv,~]=gdgsec(x);
  if iorb ~= 0, tf=[tf; t]; xf=[xf; x]; end
  if iorb < 0, fprintf('propTallSec, first step. Time and value of g: %f %e\n',t,gv); end
 end
%-----------------------------------------------------------
% Detecting crossing events with the change of sign of g(x)
%-----------------------------------------------------------
 for nt=1:ntall
  gva=gv;
  while gv*gva>0
   [t,x,h,~]=rk45f(t,x,h,hmin,hmax,tol,vfield);
   gva=gv;
   [gv,dgv]=gdgsec(x);
   if iorb ~=0 
    if iorb<0, fprintf('propTallSec, looking for ntall %d,  time and val of g: %f %e\n',nt,t,gv); end
    if gv*gva>0 && nt<=ntall tf=[tf; t]; xf=[xf; x]; end
   end
  end
 end
%------------------------------------------------------
%  Newton-Raphson for the refinement of final g(x)=0.
%------------------------------------------------------
 if ntall == 0, [gv,dgv]=gdgsec(x); end
 iter = 0
 while abs(gv) > 1.e-12
  df=vfield(t,x);
  h=-gv/dot(dgv,df(1:length(dgv)));  % implements h=-g(t)/g'(t)
  [t,x,h,~]=rk45f(t,x,h,hmin,hmax,tol,vfield);
  [gv,dgv]=gdgsec(x);
  iter = iter + 1
  fprintf('iter is\n',iter)
  if iorb < 0, fprintf('propTallSec, final Newton. Time and value of g: %f %e\n',t,gv); end

 end
 if iorb == 0, tf=t; xf=x; else, tf=[tf; t]; xf=[xf; x]; 
 end
end

