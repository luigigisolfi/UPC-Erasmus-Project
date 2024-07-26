function [gv, dgv]=gdgXSec(x)
%-----------------------------------------------------------------------
% Defines section x=0 and its gradient vector at point x
%-----------------------------------------------------------------------
gv=x(1);
dgv=[1 0 0 0 0 0];
end