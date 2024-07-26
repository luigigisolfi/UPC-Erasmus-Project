function [gv, dgv]=gdgYSec(x)
%-----------------------------------------------------------------------
% Defines section y=0 and its gradient vector at point x
%-----------------------------------------------------------------------
gv=x(2);
dgv=[0 1 0 0 0 0];
%dgv =[dgv, zeros(1, 36)];
end