function print_eigs_data(altm, mach, eigvectors, eigvalues, fid)
% PRINT_EIGS_DATA prints the output of the eigenvalue analysis of
% the system. 

% decide weather to write to a file or not
if nargin == 4
    fid = 1;
end
% Calculations
wns = abs(imag(eigvalues));
ns = real(eigvalues);
freqs = wns ./ (2*pi);
time2double = log(2) ./ abs(ns); % seee etkin p.164

%% Formatting output
fprintf(fid, '***********************************\n');
fprintf(fid, '** Altitude = %d m**\n', altm);
fprintf(fid, '** Mach = %.2f **\n', mach);

for i = 1:length(eigvalues)
    fprintf(fid, '----------------------------------\n');
    fprintf(fid, '%d. Eigenvalue = %s\n', i, num2str(eigvalues(i)));
    fprintf(fid, '----------------------------------\n');
    fprintf(fid, 'Frequency = %.4f Hz\n', freqs(i));
    fprintf(fid, 'Corresponding Eigenvector:\n');
    fprintf(fid, '\t%f\n', eigvectors(:, i));
    
    if ns(i) > 0
        fprintf(fid, 'Time to double: %f s\n', time2double(i));
    elseif ns(i) < 0
        fprintf(fid, 'Time to half: %f s\n', time2double(i));
    else
        fprintf(fid, 'n = 0; No instability.\n');
    end     
end
fprintf(fid, '***********************************\n\n');
