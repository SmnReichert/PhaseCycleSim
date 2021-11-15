% % short expample on how to use the simulation for simulating the TQTPPI
% % sequence for agarose pahntoms

B0 = 9.4; % T 

wQbar = 0; % anisotropy (mean value of wQ
wShift = 2*pi* 0; % 2*pi, wenn Omega sonst in Hz gegeben
wShift_RMS = 2*pi* 0; % 2*pi, wenn Omega sonst in Hz gegeben, 

% % 2% agar, 154mM NaCl phantom values calculated from T2 relaxation times
tauC = 1.953610551153281e-08; % single tauC
wQ = 1.497956475676493e+05;   % wQ_RMS

% % 4% agar, 154mM NaCl phantom values calculated from T2 relaxation times
% tauC = 2.489108877198124e-08;
% wQ = 1.777623103329200e+05;

% % 6% agar, 154mM NaCl phantom values calculated from T2 relaxation times
% tauC = 2.833815392207157e-08;
% wQ = 2.077809000319316e+05;

PC = PhaseCycle(B0, tauC, wQ, wQbar, wShift, wShift_RMS ); % create object of PhaseCycle class

PC.tevoStep = 0.200e-3; % set evotultion time step
PC.NumPhaseCycles = 80; % set number of phase cycles

[FIDs_p1, FIDs_p2, tevos, TQs_p1, TQs_p2, SQs_p1, SQs_p2] = PC.TQTPPI_wo180(); % TQTPPI without 180° Puls, 
[tevos, TQTPPI_FID] = PC.get_TQTPPI_FID( tevos, FIDs_p1, FIDs_p2); % create FID for TQTPPI sequence

%%
% plot TQTPPI FID
figure(); hold on
plot(tevos*1e3, real(TQTPPI_FID));
plot(tevos*1e3, imag(TQTPPI_FID));
legend("real(FID)", "imag(FID)");
xlabel("t_{evo} [ms]");
ylabel("signal [a.u.]")

% calculate TQTPPI spectrum and plot spectrum
fft_fid = fftshift(fft(TQTPPI_FID)); 
fft_fid = fft_fid/max(abs(fft_fid)); % normalize on SQ Peak

NumPoints = length(fft_fid); % create x axis (frequencies)
Fs = 1/PC.tevoStep;
dF = Fs/NumPoints;
f_vec = -Fs/2+dF:dF:Fs/2;

figure(); hold on
plot(f_vec*1e-3, real(fft_fid));
plot(f_vec*1e-3, imag(fft_fid));
legend("real(spectrum)", "imag(spectrum)");
xlabel("frequency [kHz]");
ylabel("normalized signal [a.u.]")

