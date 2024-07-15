% =================================================================
% Example for propTallSec
% ==================================================================
clear all; close all; clc;
global mu;
mu=0.0122;

%-----------------------------------------------------
% Finding Periodic Orbit (Howell's Algorithm)
%-----------------------------------------------------
ti=0;
xi=[mu-1+0.15001 0.000000000E+00 .020000000E+00 ...
0.0000000000E+00 +0.4 .0000000000E+00];
xi = [xi, zeros(1,36)]
xi(7) = 1;
xi(14) = 1;
xi(21) = 1;
xi(28) = 1;
xi(35) = 1;
xi(42) = 1;
hmin=1.e-5; hmax=1.e0; tol=1.e-10; h=0.0001;
ntall=2; iorb=-1;
VZ = 1;
VX = 1;

size(xi)
[tf, xf] = propTallSec(ti,xi(1:6),h,hmin,hmax,tol,@rtbp,@gdgYSec,ntall,iorb);

tf_half = tf(end)


fprintf('yaaaa ')

    while abs(VX) >= 1.e-6
     
     [tf, xf] = propTITF(ti,xi,tf_half,@rtbpv,hmin,hmax,tol,iorb); 
     init_42_elements = xi(:)
     last_42_elements = xf(end,:);% Extract the last 36 elements of xi
     vector_state =  last_42_elements(1:6);
     matrix_phi = reshape(last_42_elements(7:42), 6, 6).' % Reshape into a 6x6 matrix
 
     X = vector_state(1);
     Y = vector_state(2);
     Z = vector_state(3);
     VX = vector_state(4);
     VY = vector_state(5);
     VZ = vector_state(6);
     
     %Retrieve x,z accelerations and y velocity at time tf = T/2
     r12= (X-mu+1)^2 + Y^2 + Z^2; % r12: square of distance to P 
     r22= (X-mu)^2 + Y^2 + Z^2;   % r22: square of distance to S
     r13=r12*sqrt(r12);
     r23=r22*sqrt(r22);
    
     Ox = X- ((1-mu)*(X-mu)/r23 + mu*(X-mu+1)/r13);
     Oy = Y - ((1-mu)* Y    /r23 + mu* Y     /r13);
     Oz =      - ((1-mu)* Z    /r23 + mu* Z      /r13);
    
     ax = 2*VY + Ox ;
     ay = -2*VX + Oy ;
     az = Oz;


     %Construct the needed differential
     a11=matrix_phi(4,3);
     a12=matrix_phi(4,5);
     a21=matrix_phi(6,3);
     a22=matrix_phi(6,5);
     a31=matrix_phi(2,3);
     a32=matrix_phi(2,5);
     %the differential of the constraint vector
     DF=[a11, a12, ax;
         a21, a22, az;
         a31, a32, VY];
     %Invert
     D=inv(DF);
    %Compute next (or final/target) initial conditions    
    %compute next step
     Term_1 = [xi(3); xi(5); tf_half];
     tf_half
     Term_2 = -D*[vector_state(4); vector_state(6); vector_state(2)];
     x_star= Term_1 + Term_2;

     xi(3) = x_star(1)
     xi(5) = x_star(2)
     
     tf_half =x_star(3);
     
     fprintf('yooo\n')
end
%------------------------------------------------------------------
% Plot of the trajectory when it has been stored
%-------------------------------------------------------------------
iorb = 1
[tf, xf] = propTITF(ti,xi,2*10^18*tf_half,@rtbpv,hmin,hmax,tol,iorb);
xf(:)
if iorb~=0, plot(xf(:,1),xf(:,2)); end
