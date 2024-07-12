function avg_e = get_avg_e(inertial_state_primaries, PRIMARIES)
    
    GM = get_GM_body(PRIMARIES{1}); %of the central body (in the case SUN-EM Barycenter, it is the SUN)
    PRIMARY_str_2 = regexprep(PRIMARIES{2}, [{'\s+'}, {'-'}], '_');

    r_2_list = inertial_state_primaries.(PRIMARY_str_2).position; %position of PRIMARY{2}
    v_2_list = inertial_state_primaries.(PRIMARY_str_2).velocity;%velocity of PRIMARY {2}

e_list = zeros(length(r_2_list),1);
for i = 1:length(r_2_list)
    orbel = rv2cel(r_2_list(1:3,i), v_2_list(1:3,i), GM);
    e_list(i) = orbel(2);
end

avg_e = median(e_list);

end
