function [gv, dgv]=gdgY2eZ2Sec(x)
%-----------------------------------------------------------------------
% Defines section y^2-z^2=0 and its gradient vector at point x
%-----------------------------------------------------------------------
gv=x(2)*x(2)-x(3)*x(3);
dgv=[0 2*x(2) -2*x(3) 0 0 0];
fprintf('yo')
end