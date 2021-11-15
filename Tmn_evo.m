% % % % % % % % % % % % % %
% % Tmn class 
% %  Simulation of Na Dynamics
% %  Hard Pulses and Relaxation with or without a macroscopic anisotropy
% %  and a potential B0shift/inhomogeneity
% %  
% %  Tmn matrix contains the the fraction of each Tmn basis element

classdef Tmn_evo < handle
properties(Constant)
    Teq = [0 0 0 1 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0]; % thermal equilibrium = T10
end
properties
    
    Tmn = [0 0 0 1 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0];
    % % % % % % % % % % % % % % % % % % % % % 
    % % order m = 1, 2, 3;                   --> changes for relaxation, constant for pulses
    % % rank  n = -3, -2, -1, 0, 1 , 2, 3    --> changes for pulses, constant for relaxation
    % % Tmn  = [T1-3 T1-2 ... T13;
    % %         T2-3 T2-2 ... T23; 
    % %         T3-3 T3-2 ... T33; ]
    % % % % % % % % % % % % % % % % % % % % %
    
    B0 = 9.4; % T
    w0                     % Lamor frequency
    tauC                   % correlation time
    wQ                     % wQ_RMS 
    wQbar = 0;             % mean(wQ(t)) --> anisotropy
    tsDataPoints = 2048;   % Number of points for evolution period
    
    wShift = 0;            % B-field inhomogeneity -> mean(Omega)
    wShift_RMS = 0;        % B-field inhomogeneity -> ~std(Omega)


end
methods(Static)  
    % returns Tmn amplitude for given Tmn Matrix
    function [Tmn_value] = getTmn_matrix(Tmn, m, n)
        Tmn_value = Tmn(m, n+4);
    end

end
methods
    % constructor
    function obj = Tmn_evo(B0, tauC, wQ, wQbar, wShift, wShift_RMS, Tmn_value)
        if nargin == 6 % default is eq matrix = T10
            obj.Tmn = [0 0 0 1 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0];
        else
            obj.Tmn = Tmn_value;
        end
        obj.B0 = B0;
        obj.w0 = getValues.Na_gamma * B0;
        
        obj.tauC = tauC;
        obj.wQ = wQ;
        obj.wQbar = wQbar;
        
        obj.wShift = wShift;
        obj.wShift_RMS = wShift_RMS;
    end
    
    % % copy Tmn_evo oject
    function [Tmn] = copy(obj)
        Tmn = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS, obj.Tmn);
        Tmn.tsDataPoints = obj.tsDataPoints; 
    end
    
    % % getter and setter functions for a single Tmn element
    function [Tmn_value] = getTmn(obj, m, n)
        Tmn_value = obj.Tmn(m, n+4);
    end
    function setTmn(obj, m, n, Tmn_value)
        obj.Tmn(m, n+4) = Tmn_value;
    end
    

    
    % % return Tmn_evo object after pulse application
    function [Tmn_pulse] = pulse(obj, flipAngle, phase)
        Tmn_prePulse = obj.Tmn;
        % % for all m Tmn,new = (P * Tmn^T)^T      ^T = transpose operator 
        Tmn_pulse1 = (getValues.get_pulseOperator(flipAngle, phase, 1) * Tmn_prePulse(1,:).').';
        Tmn_pulse2 = (getValues.get_pulseOperator(flipAngle, phase, 2) * Tmn_prePulse(2,:).').';
        Tmn_pulse3 = (getValues.get_pulseOperator(flipAngle, phase, 3) * Tmn_prePulse(3,:).').';
        Tmn_pulse_matrix  = [Tmn_pulse1 ; Tmn_pulse2; Tmn_pulse3]  ; 
        Tmn_pulse = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS, Tmn_pulse_matrix);
    end
    
    % % relaxation 
    % % without anisotropy, B0 shift, etc
    % % return: Tmn object after relaxation, 
    % %         Tmn Matrix with time evolution
    function [Tmn_relaxation, Tmn_timeSeries] = relaxation(obj, tevo)
        Tmn_preRelax = obj.copy(); 
        Tmn_preRelax.Tmn = Tmn_preRelax.Tmn - Tmn_preRelax.Teq; % recover thermal equilibrium
        
        ts = 0:tevo/obj.tsDataPoints:tevo;
        Tmn_non = zeros([1 length(ts)]); % nonexisting Tmns
        T10_ts = ones([1 length(ts)]);   
                
        fm333  = getValues.get_f(-3, ts, obj.tauC, obj.wQ, obj.w0);
        fm2kks = getValues.get_f(-2, ts, obj.tauC, obj.wQ, obj.w0);
        fm1kks = getValues.get_f(-1, ts, obj.tauC, obj.wQ, obj.w0);
        f0kks  = getValues.get_f( 0, ts, obj.tauC, obj.wQ, obj.w0);
        f1kks  = getValues.get_f( 1, ts, obj.tauC, obj.wQ, obj.w0);
        f2kks  = getValues.get_f( 2, ts, obj.tauC, obj.wQ, obj.w0);        
        f333   = getValues.get_f( 3, ts, obj.tauC, obj.wQ, obj.w0);
          
        % zero quantum relaxation
        f011 = f0kks(1,:); f013 = f0kks(2,:); f033 = f0kks(3,:); f022 = f0kks(4,:);
        T10_pre = Tmn_preRelax.getTmn(1, 0);
        T30_pre = Tmn_preRelax.getTmn(3, 0);
        T20_pre = Tmn_preRelax.getTmn(2, 0);
        T10_relax = f011*T10_pre+f013*T30_pre + T10_ts;  % add T10_ts to recover thermal equilibrium
        T30_relax = f013*T10_pre+f033*T30_pre;
        T20_relax = f022* T20_pre;
        
        % single quantum relaxation
        fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm133 = fm1kks(3,:); fm122 = fm1kks(4,:);
        f111 = f1kks(1,:); f113 = f1kks(2,:); f133 = f1kks(3,:); f122 = f1kks(4,:);
        T1m1_pre = Tmn_preRelax.getTmn(1, -1); T11_pre = Tmn_preRelax.getTmn(1, 1);
        T3m1_pre = Tmn_preRelax.getTmn(3, -1); T31_pre = Tmn_preRelax.getTmn(3, 1);
        T2m1_pre = Tmn_preRelax.getTmn(2, -1); T21_pre = Tmn_preRelax.getTmn(2, 1);
        % n=-1
        T1m1_relax = fm111*T1m1_pre+fm113*T3m1_pre;   
        T3m1_relax = fm113*T1m1_pre+fm133*T3m1_pre;
        T2m1_relax = fm122* T2m1_pre;
        % n=+1
        T11_relax = f111*T11_pre+f113*T31_pre;   
        T31_relax = f113*T11_pre+f133*T31_pre;
        T21_relax = f122* T21_pre;
        
        % multiple quantum relaxation
        % m = 2
        f222 = f2kks(1,:); f233 = f2kks(2,:);
        fm222 = fm2kks(1,:); fm233 = fm2kks(2,:);
        T2m2_pre = Tmn_preRelax.getTmn(2, -2); T22_pre = Tmn_preRelax.getTmn(2, 2);
        T3m2_pre = Tmn_preRelax.getTmn(3, -2); T32_pre = Tmn_preRelax.getTmn(3, 2);
        % n=-2
        T2m2_relax = fm222* T2m2_pre;
        T3m2_relax = fm233* T3m2_pre;
        % n=+2
        T22_relax = f222* T22_pre;
        T32_relax = f233* T32_pre; 
        % m = 3
        T3m3_pre = Tmn_preRelax.getTmn(3, -3); T33_pre = Tmn_preRelax.getTmn(3, 3);
        T3m3_relax = fm333* T3m3_pre;
        T33_relax  = f333* T33_pre;
        
        Tmn_timeSeries = zeros([3 7 length(ts)]);
        Tmn_timeSeries(1,:,:) = [   Tmn_non;    Tmn_non; T1m1_relax; T10_relax; T11_relax;   Tmn_non;   Tmn_non];
        Tmn_timeSeries(2,:,:) = [   Tmn_non; T2m2_relax; T2m1_relax; T20_relax; T21_relax; T22_relax;   Tmn_non];
        Tmn_timeSeries(3,:,:) = [T3m3_relax; T3m2_relax; T3m1_relax; T30_relax; T31_relax; T32_relax; T33_relax];
        
        Tmn_relax = Tmn_timeSeries(:,:,end);        
        Tmn_relaxation = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS, Tmn_relax);
    end  
    
    % % return FID --> only Relaxation in Tmn(1,-1) channel
    % % receiver phase phaseRX     
    function [FID] = getFIDphase(obj, ts, phaseRX)
        Tmn_preRelax = obj.Tmn; % für T10 Relax muss Teq abgezogen und am Ende wieder addiert werden
        fm1kks = getValues.get_f_aniso(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        % f1kks = getValues.get_f(1, ts, obj.tauC, obj.wQ, obj.w0);
        % single quantum relaxation
        fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm112 = -fm1kks(4,:);
        % f111 = f1kks(1); f113 = f1kks(2); 
        T1m1_pre = Tmn_evo.getTmn_matrix(Tmn_preRelax, 1, -1); % T11_pre = Tmn_preRelax(1, 1);
        T3m1_pre = Tmn_evo.getTmn_matrix(Tmn_preRelax, 3, -1); % T31_pre = Tmn_preRelax(1, 1);
        T2m1_pre = Tmn_evo.getTmn_matrix(Tmn_preRelax, 2, -1);
        
%         T1m1_relax = fm111*T1m1_pre+fm113*T3m1_pre;   
        % T11_relax = f111*T11_pre+f113*T31_pre;
        T1m1_relax = (fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre).*exp(-1*(obj.wShift_RMS+1i*obj.wShift)*ts); 
        
        FID = T1m1_relax * exp(1i*deg2rad(phaseRX));
        
%         figure(); plot(ts, real(FID)); 
    end
    
    % % relaxation  with macroscopic anisotropy, without B0 inhomogeneity
    % % return: Tmn object after relaxation, 
    % %         Tmn Matrix with time evolution
    function [Tmn_relaxation, Tmn_timeSeries] = relaxation_aniso(obj, tevo)
        Tmn_preRelax = obj.copy(); 
        Tmn_preRelax.Tmn = Tmn_preRelax.Tmn - Tmn_preRelax.Teq;  % recover thermal equilibrium

        ts = 0:tevo/obj.tsDataPoints:tevo;
        Tmn_non = zeros([1 length(ts)]); % nonexisting Tmns
        T10_ts = ones([1 length(ts)]);   
        
        fm333  = getValues.get_f_aniso(-3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        fm2kks = getValues.get_f_aniso(-2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        fm1kks = getValues.get_f_aniso(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f0kks  = getValues.get_f_aniso( 0, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f1kks  = getValues.get_f_aniso( 1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f2kks  = getValues.get_f_aniso( 2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f333   = getValues.get_f_aniso( 3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
          
        % zero quantum relaxation 
        f011 = f0kks(1,:); f013 = f0kks(2,:); f033 = f0kks(3,:); f022 = f0kks(4,:);
        T10_pre = Tmn_preRelax.getTmn(1, 0);
        T30_pre = Tmn_preRelax.getTmn(3, 0);
        T20_pre = Tmn_preRelax.getTmn(2, 0);
        T10_relax = f011*T10_pre+f013*T30_pre + T10_ts;   % add T10_ts to recover thermal equilibrium
        T30_relax = f013*T10_pre+f033*T30_pre;
        T20_relax = f022* T20_pre;
        
        % single quantum relaxation
        fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm133 = fm1kks(3,:); fm122 = fm1kks(6,:);
        fm112 = -fm1kks(4,:); fm123 = -fm1kks(5,:);
        f111 = f1kks(1,:); f113 = f1kks(2,:); f133 = f1kks(3,:); f122 = f1kks(4,:);
        f112 = f1kks(4,:); f123 = f1kks(5,:);
        T1m1_pre = Tmn_preRelax.getTmn(1, -1); T11_pre = Tmn_preRelax.getTmn(1, 1);
        T3m1_pre = Tmn_preRelax.getTmn(3, -1); T31_pre = Tmn_preRelax.getTmn(3, 1);
        T2m1_pre = Tmn_preRelax.getTmn(2, -1); T21_pre = Tmn_preRelax.getTmn(2, 1);
        % n=-1
        T1m1_relax = fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre;   
        T3m1_relax = fm113*T1m1_pre + fm123*T2m1_pre + fm133*T3m1_pre;
        T2m1_relax = fm112*T1m1_pre + fm122*T2m1_pre + fm123*T3m1_pre;
        % n=+1
        T11_relax = f111*T11_pre + f112*T21_pre + f113*T31_pre;   
        T31_relax = f113*T11_pre + f123*T21_pre + f133*T31_pre;
        T21_relax = f112*T11_pre + f122*T21_pre + f123*T31_pre;
        
        % multiple quantum relaxation
        % m = 2
        f222 = f2kks(1,:); f233 = f2kks(2,:); f223 = f2kks(3,:); 
        fm222 = fm2kks(1,:); fm233 = fm2kks(2,:); fm223 = -f2kks(3,:);
        T2m2_pre = Tmn_preRelax.getTmn(2, -2); T22_pre = Tmn_preRelax.getTmn(2, 2);
        T3m2_pre = Tmn_preRelax.getTmn(3, -2); T32_pre = Tmn_preRelax.getTmn(3, 2);
        % n=-2
        T2m2_relax = fm222*T2m2_pre + fm223*T3m2_pre;
        T3m2_relax = fm233*T3m2_pre + fm223*T2m2_pre;
        % n=+2
        T22_relax = f222*T22_pre + f223*T32_pre;
        T32_relax = f233*T32_pre + f223*T22_pre;
        % m = 3 
        T3m3_pre = Tmn_preRelax.getTmn(3, -3); T33_pre = Tmn_preRelax.getTmn(3, 3);
        T3m3_relax = fm333* T3m3_pre;
        T33_relax  = f333* T33_pre;
        
        Tmn_timeSeries = zeros([3 7 length(ts)]);
        Tmn_timeSeries(1,:,:) = [   Tmn_non;    Tmn_non; T1m1_relax; T10_relax; T11_relax;   Tmn_non;   Tmn_non];
        Tmn_timeSeries(2,:,:) = [   Tmn_non; T2m2_relax; T2m1_relax; T20_relax; T21_relax; T22_relax;   Tmn_non];
        Tmn_timeSeries(3,:,:) = [T3m3_relax; T3m2_relax; T3m1_relax; T30_relax; T31_relax; T32_relax; T33_relax];
        
        Tmn_relax = Tmn_timeSeries(:,:,end);        
        Tmn_relaxation = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS, Tmn_relax);
    end 
 


    % % relaxation  with macroscopic anisotropy and chemical shift (Cauchy/Lorentzian distribution for offset)
    % % return: Tmn object after relaxation, 
    % %         Tmn Matrix with time evolution    
     function [Tmn_relaxation, Tmn_timeSeries] = relaxation_shift(obj, tevo)
        Tmn_preRelax = obj.copy(); 
        Tmn_preRelax.Tmn = Tmn_preRelax.Tmn - Tmn_preRelax.Teq; % recover thermal equilibrium
 
        ts = 0:tevo/obj.tsDataPoints:tevo;
        Tmn_non = zeros([1 length(ts)]); % nonexisting Tmns
        T10_ts = ones([1 length(ts)]);  
  
        fm333  = getValues.get_f_aniso(-3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        fm2kks = getValues.get_f_aniso(-2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        fm1kks = getValues.get_f_aniso(-1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f0kks  = getValues.get_f_aniso( 0, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f1kks  = getValues.get_f_aniso( 1, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f2kks  = getValues.get_f_aniso( 2, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
        f333   = getValues.get_f_aniso( 3, ts, obj.tauC, obj.wQ, obj.wQbar, obj.w0);
          
        % zero quantum relaxation 
        f011 = f0kks(1,:); f013 = f0kks(2,:); f033 = f0kks(3,:); f022 = f0kks(4,:);
        T10_pre = Tmn_preRelax.getTmn(1, 0);
        T30_pre = Tmn_preRelax.getTmn(3, 0);
        T20_pre = Tmn_preRelax.getTmn(2, 0);
        % T10_30_relax = [f011 f013; f013; f033]*[T10_pre; T30_pre];
        % T10_relax = T10_30_relax(1);   T30_relax = T10_30_relax(2);
        T10_relax = f011*T10_pre+f013*T30_pre + T10_ts;   % add T10_ts to recover thermal equilibrium
        T30_relax = f013*T10_pre+f033*T30_pre;
        T20_relax = f022* T20_pre;
        
        % single quantum relaxation
        fm111 = fm1kks(1,:); fm113 = fm1kks(2,:); fm133 = fm1kks(3,:); fm122 = fm1kks(6,:);
        fm112 = -fm1kks(4,:); fm123 = -fm1kks(5,:);
        f111 = f1kks(1,:); f113 = f1kks(2,:); f133 = f1kks(3,:); f122 = f1kks(4,:);
        f112 = f1kks(4,:); f123 = f1kks(5,:);
        T1m1_pre = Tmn_preRelax.getTmn(1, -1); T11_pre = Tmn_preRelax.getTmn(1, 1);
        T3m1_pre = Tmn_preRelax.getTmn(3, -1); T31_pre = Tmn_preRelax.getTmn(3, 1);
        T2m1_pre = Tmn_preRelax.getTmn(2, -1); T21_pre = Tmn_preRelax.getTmn(2, 1);
        % n=-1
        T1m1_relax = (fm111*T1m1_pre + fm112*T2m1_pre + fm113*T3m1_pre).*exp(-1*(obj.wShift_RMS+1i*obj.wShift)*ts);   
        T3m1_relax = (fm113*T1m1_pre + fm123*T2m1_pre + fm133*T3m1_pre).*exp(-1*(obj.wShift_RMS+1i*obj.wShift)*ts);
        T2m1_relax = (fm112*T1m1_pre + fm122*T2m1_pre + fm123*T3m1_pre).*exp(-1*(obj.wShift_RMS+1i*obj.wShift)*ts);
        % n=+1
        T11_relax = (f111*T11_pre + f112*T21_pre + f113*T31_pre).*exp(-1*(obj.wShift_RMS-1i*obj.wShift)*ts);   
        T31_relax = (f113*T11_pre + f123*T21_pre + f133*T31_pre).*exp(-1*(obj.wShift_RMS-1i*obj.wShift)*ts);
        T21_relax = (f112*T11_pre + f122*T21_pre + f123*T31_pre).*exp(-1*(obj.wShift_RMS-1i*obj.wShift)*ts);
        
        % multiple quantum relaxation
        % m = 2
        f222 = f2kks(1,:); f233 = f2kks(2,:); f223 = f2kks(3,:); 
        fm222 = fm2kks(1,:); fm233 = fm2kks(2,:); fm223 = f2kks(3,:);
        T2m2_pre = Tmn_preRelax.getTmn(2, -2); T22_pre = Tmn_preRelax.getTmn(2, 2);
        T3m2_pre = Tmn_preRelax.getTmn(3, -2); T32_pre = Tmn_preRelax.getTmn(3, 2);
        % n=-2
        T2m2_relax = (fm222*T2m2_pre + fm223*T3m2_pre).*exp(-2*(obj.wShift_RMS+1i*obj.wShift)*ts);
        T3m2_relax = (fm233*T3m2_pre + fm223*T2m2_pre).*exp(-2*(obj.wShift_RMS+1i*obj.wShift)*ts);
        % n=+2
        T22_relax = (f222*T22_pre + f223*T32_pre).*exp(-2*(obj.wShift_RMS-1i*obj.wShift)*ts);
        T32_relax = (f233*T32_pre + f223*T22_pre).*exp(-2*(obj.wShift_RMS-1i*obj.wShift)*ts);
        % m = 3 
        T3m3_pre = Tmn_preRelax.getTmn(3, -3); T33_pre = Tmn_preRelax.getTmn(3, 3);
        T3m3_relax = fm333* T3m3_pre.*exp(-3*(obj.wShift_RMS+1i*obj.wShift)*ts);
        T33_relax  = f333* T33_pre.*exp(-3*(obj.wShift_RMS-1i*obj.wShift)*ts);

        Tmn_timeSeries = zeros([3 7 length(ts)]);
        Tmn_timeSeries(1,:,:) = [   Tmn_non;    Tmn_non; T1m1_relax; T10_relax; T11_relax;   Tmn_non;   Tmn_non];
        Tmn_timeSeries(2,:,:) = [   Tmn_non; T2m2_relax; T2m1_relax; T20_relax; T21_relax; T22_relax;   Tmn_non];
        Tmn_timeSeries(3,:,:) = [T3m3_relax; T3m2_relax; T3m1_relax; T30_relax; T31_relax; T32_relax; T33_relax];
        
        Tmn_relax = Tmn_timeSeries(:,:,end);        
        Tmn_relaxation = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS, Tmn_relax);
    end 
end
end
