%function [f_value, correction] = parallel_shooting(N, orbit_file)
    global mu
    N = 60 %number of nodes
    mu = 0.0122;
    orbit_file = 'Halo.txt';

    %mu = 0.1 %Masde initial conditions
    %orbit_file = 'orbit_data.txt' Masde orbits
   
    % Retrieve (read) orbit and add noise
    [t, x] = read_orbit(orbit_file);
    
    % Add noise to all orbit
    x_size = size(x);
    x_noised = x + randn(x_size) * 1e-5;
    ti = t(1);
    tf = t(end);
    
    eps = 1e-5;
    hmin = 1.e-4;
    hmax = 1.e0;
    tol = 1.e-10;

    % Create interpolation for both times and positions
    interpol = struct();
    interpol.x = t;
    interpol.y = x_noised;
    interpol.spline = cell(1, 1); 

    % Sample Points
    t_sampled = linspace(ti, tf, N);
    x_sampled = zeros(6, length(t_sampled));
    for dim = 1:6
       interpol.spline{dim} = spline(t, x_noised(:,dim).');
       x_sampled(dim,:) = ppval(interpol.spline{dim}, t_sampled);
    end

    %stop = 'N';
    Q0 = x_sampled;
    for iteration = 1:10
    
    fprintf('iteration %f\n', iteration)
    t_list_ = [];
    F_list = [];
    df = cell(N-1, N);

    if iteration == 10
    figure;
    hold on;
    %scatter(x(:,1), x(:,2), 'HandleVisibility','off', 'Color', 'Black');
    
    end
    fprintf('Q0 to recover %s\n', mat2str(x(1,:)))
    fprintf('Q0 old %s\n', mat2str(round(Q0(:,1), 6)))
    for i = 1:N-1
        
        [t_, phi_Q] = ode78(@HFEM_rtbp, [t_sampled(i), t_sampled(i+1)], Q0(:,i));
        scatter(Q0(1,:),Q0(2,:), 'filled')
        plot(phi_Q(:,1), phi_Q(:,2), 'HandleVisibility','off', 'Color', 'Green')
        t_list_ = [t_list_; t_];
        F = phi_Q(end,:).' - Q0(:,i+1);
        F_list = [F_list; F];
        stm_x_noised_ = numericSTMvfield(t_sampled(i),t_sampled(i+1), Q0(:,i), eps, @HFEM_rtbp,hmin,hmax,tol);
        %xiv=zeros(1,42);
        %xiv(1:6)=Q0(:,i).';
        %for p= 1:6 , xiv(7*p)=1.e0; end  % variational matrix initalized with Id.
        %x_ = propTITF_vfield(t_sampled(i),xiv,t_sampled(i+1),@rtbpv,hmin,hmax,tol);
        %stm_x_noised__ = x_(7:42)
        stm_6x6 = reshape(stm_x_noised_, [6,6]);
        df{i,i} = stm_6x6;
        df{i,i+1} = -eye(6);
        for j = 1:N
            if j ~= i && j ~= i+1
                df{i,j} = zeros(6);
            end
        end
    end
    
    % Initialize M with zeros (or appropriate initial values)
    [row_size, col_size] = size(df{1, 1});  % Assuming all matrices are of the same size
    DF = zeros(size(df, 1) * row_size, size(df, 2) * col_size);

    % Populate DF with data from df
    for i = 1:N-1
        for j = 1:(N)
            % Compute the starting index for each matrix in DF
            start_row = (i - 1) * row_size + 1;
            start_col = (j - 1) * col_size + 1;

            % Assign the matrix from df to the corresponding position in DF
            DF(start_row:start_row + row_size - 1, start_col:start_col + col_size - 1) = df{i, j};
        end
    end

    delta_Q = - DF.' * inv(DF*DF.') * F_list;
    delta_Q = reshape(delta_Q, [6,N]);
    norm(delta_Q.')
    fprintf('Q0 old %f\n', Q0(1,1))
    Q0 = Q0 + delta_Q;
    fprintf('Q0 new %f\n', Q0(1,1))
    scatter(Q0(1,:),Q0(2,:), 'filled')
    end 