function lagrange_points_example
    % These masses represent the Earth-Moon system
    global mu

    function y = collinear_lagrange(xstar)
    r12= (xstar- mu+1)^2; % r12: square of distance to P 
    r22= (xstar-mu)^2;  % r22: square of distance to S
    r13=r12*sqrt(r12);
    r23=r22*sqrt(r22);
    
    y = xstar- ((1-mu)*(xstar-mu)/r23 + mu*(xstar-mu+1)/r13);

    end
    options = optimset('Display','iter');
    L_2 = fzero(@collinear_lagrange,[-1, -1.5],options);
    L_1 = fzero(@collinear_lagrange,[-0.5, -0.97],options);
    L_3 = fzero(@collinear_lagrange,[0.5, 3],options);
    fprintf('L_1=%f, L_2=%f, L_3=%f\n', L_1, L_2, L_3)
    fprintf('L_1_from_Moon = %f, L_2_from_Moon=%f, L_3_from_Earth =%f\n', abs(L_1 - mu +1), abs(L_2 - mu + 1), abs(L_3 - mu))

    function output
        figure()
        title('Earth-Moon')
        hold on
        ylim([-1.2,1.2])
        xlim([-1.2,1.2])
        axis equal
        xlabel('x')
        ylabel('y')
        plot(mu ,0,'bo','MarkerFaceColor','b')
        plot(mu-1,0,'go','MarkerFaceColor','g')
        plot(L_1,0,'rv','MarkerFaceColor','r')
        plot(L_2,0,'r^','MarkerFaceColor','r')
        plot(L_3,0,'rp','MarkerFaceColor','r')
        plot(0.5 + mu- 1,sqrt(3)/2,'rX','MarkerFaceColor','r')
        plot(0.5 + mu- 1,-sqrt(3)/2,'rs','MarkerFaceColor','r')
        %yline(0,'k')
        %coords = linspace(0,pi,100);
        %plot((1-mu).*cos(coords),(1-mu).*sin(coords))
        hold all
        %plot((1-mu).*cos(coords),-(1-mu).*sin(coords))
    end
    output
end