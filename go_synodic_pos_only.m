function [rtbp_pos_spacecraft] = go_synodic_pos_only(rel_pos, rel_vel, SEb_pos, inertial_pos)

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nPerforming conversion from synodic (adimensional) system to inertial (physical) system...\n')
rtbp_pos_spacecraft = zeros(3,length(rel_pos));

check_inertial_pos = zeros(3, length(rel_pos));
check_inertial_vel = zeros(3, length(rel_pos));

size(inertial_pos)
size(rel_pos)
for i = 1:length(rel_pos)
    rp_rs = rel_pos(:, i); 
    k = norm(rp_rs);
    vp_vs = rel_vel(:, i); 

    C = construct_C(rp_rs, vp_vs);

    b = SEb_pos(1:3, i);

    
    %Apply the formula to conversion between two reference systems 
    %(see "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points", pp. 137-138)

    rtbp_pos_spacecraft(:,i) =  C\(inertial_pos(i,1:3).'-b)/k; 

    %--------------------------------------------------------------------------------------------%   
    % These two arrays are needed for the check of positions and velocities
    % converted back to rtbp coordinates (if requested)
    check_inertial_pos(:,i) = k*C*rtbp_pos_spacecraft(:,i)+ b ;
   
end

%fprintf('Conversion: Done.\n')
%--------------------------------------------------------------------------------------------%
% Perform a quick check that the obtained inertial positions and velocities
% can be convertedback into the original rtbp positions and velocities and
% they are the same if converted back to rtbp coordinates

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nChecking whether the conversion has been successfull...\n')
array_to_check = check_inertial_pos - inertial_pos(:,1:3).';
check = 0;
for i = 1:3
   for j = 1:length(array_to_check)
      if array_to_check(i,j) > 1e-7
          array_to_check(i,j)
          check = 1 ;
      end
   end
end

if check == 1
   error('Attention. Conversion check has not been successfull. Aborting...\n')
   fprintf('-----------------------------------------------------------\n')
else
   fprintf('All good.\n')
end
%--------------------------------------------------------------------------------------------%
