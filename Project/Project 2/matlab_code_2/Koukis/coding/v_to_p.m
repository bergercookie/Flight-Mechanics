function data_p = v_to_p(data_v)
% V_TO_P takes the voltage data and converts it to p_values.
% p_values are in [deg/sec]
data_p = data_v/(20*10^-3);
end