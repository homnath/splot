% Input parameters time
t0_src = 0.0;   % Origin time (Sec)
f  = 50.0;      % Frequency (Hz)
dt = 1e-4;      % Time interval (Sec)
nt = 5000;      % Total time steps
%--------------------------------------------------------------------------

t0 = 1.0/f + t0_src;  
t = dt*(0:nt-1);


factor = (pi*f*(t - t0)).*(pi*f*(t - t0));

 % Ricker with maximum amplitude 1
stf = (2.*factor - 1.).*exp(-factor);
% Test
fprintf(1,'Ricker source: Maximum amplitude = %f\n',max(abs(stf)));
% Write
inpf = fopen('source_time_function_ricker.txt','w');
fprintf(inpf,'%.6e %.6e\n',[t; stf]);
fclose(inpf);
% Plot
figure
plot(t,stf,'k','LineWidth',1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Ricker source with maximum amplitude 1');

% Gaussian with maximum amplitude 1
stf = exp(-factor);
% Test
fprintf(1,'Gaussian source: Maximum amplitude = %f\n',max(abs(stf)));
% Write
inpf = fopen('source_time_function_gauss_max1.txt','w');
fprintf(inpf,'%.6e %.6e\n',[t; stf]);
fclose(inpf);
% Plot
figure
plot(t,stf,'k');
xlabel('Time (s)');
ylabel('Amplitude');
title('Gaussian source with maximum amplitude 1');

 % Gaussian with integration 1
stf = sqrt(pi)*f*exp(-factor);
% Test
fprintf(1,'Gaussian source: Integration = %f\n',sum(stf)*dt);
% Write
inpf = fopen('source_time_function_gauss_int1.txt','w');
fprintf(inpf,'%.6e %.6e\n',[t; stf]);
fclose(inpf);
% Plot
figure
plot(t,stf,'k');
xlabel('Time (s)');
ylabel('Amplitude');
title('Gaussian source with integration 1');


