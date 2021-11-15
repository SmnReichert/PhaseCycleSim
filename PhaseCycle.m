% % % % % % % % % % % % % %
% % Phase cyle class 
% %  Simulation of phase cycle
% %  - soft pulses are not available
% %  - use relaxation_shift to simulate macroscopic anisotropy and B0
% %    inhomogeneities 

classdef PhaseCycle < handle
properties
    B0 = 9.4 % T
    
    w0
    tauC
    wQ
    wQbar = 0
    
    wShift = 0;
    wShift_RMS = 0;
        
    tevoStep = 0.1000 * 1e-3; % s
    tevo0    = 0.1305 * 1e-3; % s 
    tmix     = 0.1500 * 1e-3; % s

    startPhase = 90; % degrees
    alphas 

    NumPhaseCycles 
    PulseAngles = []
    
    TR = 300e-3; % s  repetition time
    
    timeFID = 200e-3; % s
    dataPoints = 2048;
    
    deadtimeFID =   6e-6; % s
    dwelltimeFID = 50e-6; % s
    
    
    flip90  = 90; % can be set to different values to simulate B1+ 
    flip180 = 180;% inhomogeneities
    

end
methods(Static)  
    % % general 3 Pulse sequence phase step
    % % phase of 3rd pulse is always 0 
    function [TmnEnd, TQ, SQ] = general3Pulse( TmnStart, alpha, beta, tevo, tmix, flip1, flip2, flip3)
        [Tmn_pulse90_1] = TmnStart.pulse(flip1, alpha);
        [Tmn_evo, Tmn_timeSeries_evo_1] = Tmn_pulse90_1.relaxation_shift(tevo);
        [Tmn_pulse90_2] = Tmn_evo.pulse(flip2, beta);
        [Tmn_evo_2, Tmn_timeSeries_mix] = Tmn_pulse90_2.relaxation_shift(tmix);
        [Tmn_pulse90_3] = Tmn_evo_2.pulse(flip3, 0);

        TmnEnd = Tmn_pulse90_3;
        
        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_evo_2.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, 90); 
        TQm = Tmn_evo_2.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, 90);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_evo_2.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, 90) + Tmn_evo_2.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, 90); 
        SQm = Tmn_evo_2.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, 90) + Tmn_evo_2.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, 90); 
        SQ = SQp + SQm;
    end  
    
    % % general 3 Pulse sequence phase step with 180° refocussing pulse
    % % phase of 3rd pulse is always 0 
    function [TmnEnd, TQ, SQ] = general3Pulse_w180(TmnStart, alpha, beta, tevo, tmix, flip1, flip2, flip3, flipRefocus)
        [Tmn_pulse90_1] = TmnStart.pulse(flip1, alpha);
        [Tmn_evo_1, Tmn_timeSeries_evo_1] = Tmn_pulse90_1.relaxation_shift(tevo/2);
        [Tmn_pulse180] = Tmn_evo_1.pulse(flipRefocus, alpha+90);
        [Tmn_evo_2, Tmn_timeSeries_evo_2] = Tmn_pulse180.relaxation_shift(tevo/2);
        [Tmn_pulse90_2] = Tmn_evo_2.pulse(flip2, beta);
        [Tmn_evo_3, Tmn_timeSeries_mix] = Tmn_pulse90_2.relaxation_shift(tmix);
        [Tmn_pulse90_3] = Tmn_evo_3.pulse(flip3, 0);
        TmnEnd = Tmn_pulse90_3;
        
        % calculate TQ signal (T3-3, T33 after mixing period)
        TQp = Tmn_evo_2.getTmn(3, +3)*getValues.get_WignerD(3, -1, +3, 90); 
        TQm = Tmn_evo_2.getTmn(3, -3)*getValues.get_WignerD(3, -1, -3, 90);
        TQ = TQp+TQm;
        
        % calculate SQ signal (T1-1, T11, T3-1, T31 after mixing period)
        SQp = Tmn_evo_2.getTmn(1, +1)*getValues.get_WignerD(1, -1, +1, 90) + Tmn_evo_2.getTmn(3, +1)*getValues.get_WignerD(3, -1, +1, 90); 
        SQm = Tmn_evo_2.getTmn(1, -1)*getValues.get_WignerD(1, -1, -1, 90) + Tmn_evo_2.getTmn(3, -1)*getValues.get_WignerD(3, -1, -1, 90); 
        SQ = SQp + SQm;
    end

    
    % % correct phase of 1st dimenstion FIDs
    function [FID] = phase_cor_FID(fid)
        spectrum = fftshift(fft(fid(1,:)));
        [~,pos] = max(abs(spectrum));
        banana = @(x) -real(spectrum(pos)*exp(1i*x));
        [x(1),~] = fminsearch(banana,[0]);
        if max(real(fid(1,:)* exp(1i*x(1)))) < 0
            FID = fid * exp(-1i*x(1));
        else
            FID = fid * exp(1i*x(1));
        end
    end
    
    % % TQTPPI FID
    % % both FIDs are added for DQ suppression
    % % Peak hight of 1st dimension spectra used for 2nd dimension FID
    function [tevos, TQTPPI_FID] = get_TQTPPI_FID(tevos, FIDs_p1, FIDs_p2)
        sizeP1 = size(FIDs_p1);
        FIDs = FIDs_p1 + FIDs_p2;
        
        FIDs = PhaseCycle.phase_cor_FID(FIDs);
        
        FIDs_ft =  zeros(sizeP1);
        for j = 1:sizeP1(1)
            FID = FIDs(j,:);
            spec = fftshift(fft(FID.'));
            FIDs_ft(j,:) =   spec;
        end
        maxFIDs = FIDs_ft(:,fix(sizeP1(2)/2)+1);
        TQTPPI_FID = maxFIDs;        
    end 
    
    
end
methods
    % constructor
    function obj = PhaseCycle(B0, tauC, wQ, wQbar, wShift, wShift_RMS )
        obj.wQbar = wQbar;
        obj.B0 = B0;
        obj.tauC = tauC;
        obj.wQ = wQ;
        obj.alphas =  [0:45:359] + obj.startPhase; % degrees   
        
        obj.wShift = wShift;
        obj.wShift_RMS = wShift_RMS;
        obj.NumPhaseCycles = 90;
    end   
   
    % % complete TQTPPI sequence with 180° refocussing pulse and 
    % % DQ filter: Run sequence twice with phase shift of 180° and add FIDs 
    % % returns: FIDs_p1  : FID with beta = alpha + 90°
    % %          FIDs_p2  : FID with beta = alpha - 90°
    % %          tevos    : vector of evolution times
    % %          TQs_p1,2 : TQ signal (T33+T3-3 -->obs pulse--> T3-1)  
    function [ FIDs_p1, FIDs_p2, tevos, TQs_p1, TQ_p2, SQs_p1, SQs_p2] = TQTPPI_w180(obj)
        TmnStart = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS);
        lengthAlpha = length(obj.alphas);
        times = 0:obj.timeFID/obj.dataPoints:obj.timeFID;
        Tmn_p1 = TmnStart;
        Tmn_p2 = TmnStart;
        Tmns_p1 = [];
        Tmns_p2 = [];
        FIDs_p1 = [];
        FIDs_p2 = [];
        tevos = zeros([1 obj.NumPhaseCycles*lengthAlpha]);
        for j = 1:obj.NumPhaseCycles
            for alphaIndex = 1:lengthAlpha
                tevo = obj.tevo0 +  ( (j-1)*lengthAlpha + alphaIndex - 1 ) *obj.tevoStep; 
                tevos((j-1)*lengthAlpha + alphaIndex ) = tevo;
                alpha = obj.alphas(alphaIndex);
                beta_p1 = alpha + 90;
                beta_p2 = alpha - 90;
                
                [Tmn_p1, TQ_p1, SQ_p1] = PhaseCycle.general3Pulse_w180(Tmn_p1, alpha, beta_p1, tevo, obj.tmix, obj.flip90, obj.flip90, obj.flip90, obj.flip180);
                [Tmn_p2, TQ_p2, SQ_p2] = PhaseCycle.general3Pulse_w180(Tmn_p2, alpha, beta_p2, tevo, obj.tmix, obj.flip90, obj.flip90, obj.flip90, obj.flip180);
                
                Tmns_p1 = [Tmns_p1; Tmn_p1];
                Tmns_p2 = [Tmns_p2; Tmn_p2]; 
                
                FIDs_p1 = [FIDs_p1; Tmn_p1.getFIDphase(times, 0)];
                FIDs_p2 = [FIDs_p2; Tmn_p2.getFIDphase(times, 0)]; 
                
                TQs_p1((j-1)*lengthAlpha + alphaIndex ) = TQ_p1;
                TQs_p2((j-1)*lengthAlpha + alphaIndex ) = TQ_p2; 
                
                SQs_p1((j-1)*lengthAlpha + alphaIndex ) = SQ_p1;
                SQs_p2((j-1)*lengthAlpha + alphaIndex ) = SQ_p2; 
                
                [Tmn_p1, Tmn_timeSeries_TR] =  Tmn_p1.relaxation_shift(obj.TR);
                [Tmn_p2, Tmn_timeSeries_TR] =  Tmn_p2.relaxation_shift(obj.TR);                      
            end            
        end
    end 

    % % complete TQTPPI sequence without 180° refocussing pulse and 
    % % DQ filter: Run sequence twice with phase shift of 180° and add FIDs 
    % % returns: FIDs_p1  : FID with beta = alpha + 90°
    % %          FIDs_p2  : FID with beta = alpha - 90°
    % %          tevos    : vector of evolution times
    % %          TQs_p1,2 : TQ signal (T33+T3-3 -->obs pulse--> T3-1)  
    function [ FIDs_p1, FIDs_p2, tevos, TQs_p1, TQ_p2, SQs_p1, SQs_p2] = TQTPPI_wo180(obj)
        TmnStart = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS);
        lengthAlpha = length(obj.alphas);
        times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);
        Tmn_p1 = TmnStart;
        Tmn_p2 = TmnStart;
        Tmns_p1 = [];
        Tmns_p2 = [];
        FIDs_p1 = [];
        FIDs_p2 = [];
        tevos = zeros([1 obj.NumPhaseCycles*lengthAlpha]);

        for j = 1:obj.NumPhaseCycles
            for alphaIndex = 1:lengthAlpha
                tevo = obj.tevo0 +  ( (j-1)*lengthAlpha + alphaIndex - 1 ) *obj.tevoStep; 
                tevos((j-1)*lengthAlpha + alphaIndex ) = tevo;
                alpha = obj.alphas(alphaIndex);
                beta_p1 = alpha + 90;
                beta_p2 = alpha - 90;

                [Tmn_p1, TQ_p1,SQ_p1] = PhaseCycle.general3Pulse(Tmn_p1, alpha, beta_p1, tevo, obj.tmix, obj.flip90, obj.flip90, obj.flip90);
                [Tmn_p2, TQ_p2,SQ_p2] = PhaseCycle.general3Pulse(Tmn_p2, alpha, beta_p2, tevo, obj.tmix, obj.flip90, obj.flip90, obj.flip90);
                
                Tmns_p1 = [Tmns_p1; Tmn_p1];
                Tmns_p2 = [Tmns_p2; Tmn_p2]; 

                FIDs_p1 = [FIDs_p1; Tmn_p1.getFIDphase(times, 0)];
                FIDs_p2 = [FIDs_p2; Tmn_p2.getFIDphase(times, 0)]; 
                
                TQs_p1((j-1)*lengthAlpha + alphaIndex ) = TQ_p1;
                TQs_p2((j-1)*lengthAlpha + alphaIndex ) = TQ_p2;        
                
                SQs_p1((j-1)*lengthAlpha + alphaIndex ) = SQ_p1;
                SQs_p2((j-1)*lengthAlpha + alphaIndex ) = SQ_p2; 
                
                [Tmn_p1, Tmn_timeSeries_TR] =  Tmn_p1.relaxation_shift(obj.TR);
                [Tmn_p2, Tmn_timeSeries_TR] =  Tmn_p2.relaxation_shift(obj.TR);            
            end   
        end
    end 

    
    % % complete IRTQTPPI sequence without DQ filter 
    % % returns: FIDs_p1  : FID with beta = 0°
    % %          FIDs_p2  : FID with beta = 0°
    % %          tevos    : vector of evolution times
    % %          TQs_p1,2 : TQ signal (T33+T3-3 -->obs pulse--> T3-1)  
    function [ FIDs_p1, FIDs_p2, tevos, TQs_p1, TQ_p2, SQs_p1, SQs_p2] = IRTQTPPI(obj)
        TmnStart = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS);
        lengthAlpha = length(obj.alphas);
        times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);
        Tmn_p1 = TmnStart;
        Tmn_p2 = TmnStart;
        Tmns_p1 = [];
        Tmns_p2 = [];
        FIDs_p1 = [];
        FIDs_p2 = [];
        tevos = zeros([1 obj.NumPhaseCycles*lengthAlpha]);

        for j = 1:obj.NumPhaseCycles
            for alphaIndex = 1:lengthAlpha
                tevo = obj.tevo0 +  ( (j-1)*lengthAlpha + alphaIndex - 1 ) *obj.tevoStep; 
                tevos((j-1)*lengthAlpha + alphaIndex ) = tevo;
                
                alpha_p1 = 0;
                alpha_p2 = 0;
                beta_p1 = obj.alphas(alphaIndex);
                beta_p2 = obj.alphas(alphaIndex);

                [Tmn_p1, TQ_p1, SQ_p1] = PhaseCycle.general3Pulse(Tmn_p1, alpha_p1, beta_p1, tevo, obj.tmix, obj.flip180, obj.flip90, obj.flip90);
                [Tmn_p2, TQ_p2, SQ_p2] = PhaseCycle.general3Pulse(Tmn_p2, alpha_p2, beta_p2, tevo, obj.tmix, obj.flip180, obj.flip90, obj.flip90);
                
                Tmns_p1 = [Tmns_p1; Tmn_p1];
                Tmns_p2 = [Tmns_p2; Tmn_p2]; 

                FIDs_p1 = [FIDs_p1; Tmn_p1.getFIDphase(times, 0)];
                FIDs_p2 = [FIDs_p2; Tmn_p2.getFIDphase(times, 0)]; 
                
                TQs_p1((j-1)*lengthAlpha + alphaIndex ) = TQ_p1;
                TQs_p2((j-1)*lengthAlpha + alphaIndex ) = TQ_p2;   
                
                SQs_p1((j-1)*lengthAlpha + alphaIndex ) = SQ_p1;
                SQs_p2((j-1)*lengthAlpha + alphaIndex ) = SQ_p2;   
                
                [Tmn_p1, Tmn_timeSeries_TR] =  Tmn_p1.relaxation_shift(obj.TR);
                [Tmn_p2, Tmn_timeSeries_TR] =  Tmn_p2.relaxation_shift(obj.TR); 
                
            end   
        end
    end  
    
    % % complete IRTQTPPI sequence with "0 180" DQ filter
    % % DQ filter: Run sequence twice with excitation pulse phase
    % % alternating between 0 and 180°
    % % returns: FIDs_p1  : FID with alpha = 0° 
    % %          FIDs_p2  : FID with alpha = 180°
    % %          tevos    : vector of evolution times
    % %          TQs_p1,2 : TQ signal (T33+T3-3 -->obs pulse--> T3-1)  
    function [ FIDs_p1, FIDs_p2, tevos, TQs_p1, TQ_p2, SQs_p1, SQs_p2] = IRTQTPPI_0180DQsupr(obj)
        TmnStart = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS); 
        lengthAlpha = length(obj.alphas);
        times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);
        Tmn_p1 = TmnStart;
        Tmn_p2 = TmnStart;
        Tmns_p1 = [];
        Tmns_p2 = [];
        FIDs_p1 = [];
        FIDs_p2 = [];
        tevos = zeros([1 obj.NumPhaseCycles*lengthAlpha]);

        for j = 1:obj.NumPhaseCycles
            for alphaIndex = 1:lengthAlpha
                tevo = obj.tevo0 +  ( (j-1)*lengthAlpha + alphaIndex - 1 ) *obj.tevoStep; 
                tevos((j-1)*lengthAlpha + alphaIndex ) = tevo;
                
                alpha_p1 = 0;
                alpha_p2 = 180;
                beta_p1 = obj.alphas(alphaIndex);
                beta_p2 = obj.alphas(alphaIndex);

                [Tmn_p1, TQ_p1, SQ_p1] = PhaseCycle.general3Pulse(Tmn_p1, alpha_p1, beta_p1, tevo, obj.tmix, obj.flip180, obj.flip90, obj.flip90);
                [Tmn_p2, TQ_p2, SQ_p2] = PhaseCycle.general3Pulse(Tmn_p2, alpha_p2, beta_p2, tevo, obj.tmix, obj.flip180, obj.flip90, obj.flip90);
                
                Tmns_p1 = [Tmns_p1; Tmn_p1];
                Tmns_p2 = [Tmns_p2; Tmn_p2]; 

                FIDs_p1 = [FIDs_p1; Tmn_p1.getFIDphase(times, 0)];
                FIDs_p2 = [FIDs_p2; Tmn_p2.getFIDphase(times, 0)]; 
                
                TQs_p1((j-1)*lengthAlpha + alphaIndex ) = TQ_p1;
                TQs_p2((j-1)*lengthAlpha + alphaIndex ) = TQ_p2;   
                
                SQs_p1((j-1)*lengthAlpha + alphaIndex ) = SQ_p1;
                SQs_p2((j-1)*lengthAlpha + alphaIndex ) = SQ_p2;   
                
                [Tmn_p1, Tmn_timeSeries_TR] =  Tmn_p1.relaxation_shift(obj.TR);  
                [Tmn_p2, Tmn_timeSeries_TR] =  Tmn_p2.relaxation_shift(obj.TR);                   
            end   
        end
    end    
    
    % % complete IRTQTPPI sequence with "PhaseCycle" DQ filter
    % % DQ filter: Run sequence twice with phase shift of 270° and add FIDs
    % % returns: FIDs_p1  : FID with beta + 135° 
    % %          FIDs_p2  : FID with beta - 135°
    % %          tevos    : vector of evolution times
    % %          TQs_p1,2 : TQ signal (T33+T3-3 -->obs pulse--> T3-1)  
    function [ FIDs_p1, FIDs_p2, tevos, TQs_p1, TQ_p2, SQs_p1, SQs_p2] = IRTQTPPI_PCDQsupr(obj)
        TmnStart = Tmn_evo(obj.B0, obj.tauC, obj.wQ, obj.wQbar, obj.wShift, obj.wShift_RMS);
        lengthAlpha = length(obj.alphas);
        times = obj.deadtimeFID + obj.dwelltimeFID.*(0:obj.dataPoints-1);
        Tmn_p1 = TmnStart;
        Tmn_p2 = TmnStart;
        Tmns_p1 = [];
        Tmns_p2 = [];
        FIDs_p1 = [];
        FIDs_p2 = [];
        tevos = zeros([1 obj.NumPhaseCycles*lengthAlpha]);

        for j = 1:obj.NumPhaseCycles
            for alphaIndex = 1:lengthAlpha
                tevo = obj.tevo0 +  ( (j-1)*lengthAlpha + alphaIndex - 1 ) *obj.tevoStep; 
                tevos((j-1)*lengthAlpha + alphaIndex ) = tevo;
                
                alpha_p1 = 0;
                alpha_p2 = 0;
                beta_p1 = obj.alphas(alphaIndex) + 135;
                beta_p2 = obj.alphas(alphaIndex) - 135;

                [Tmn_p1, TQ_p1, SQ_p1] = PhaseCycle.general3Pulse(Tmn_p1, alpha_p1, beta_p1, tevo, obj.tmix, obj.flip180, obj.flip90, obj.flip90);
                [Tmn_p2, TQ_p2, SQ_p2] = PhaseCycle.general3Pulse(Tmn_p2, alpha_p2, beta_p2, tevo, obj.tmix, obj.flip180, obj.flip90, obj.flip90);

                Tmns_p1 = [Tmns_p1; Tmn_p1];
                Tmns_p2 = [Tmns_p2; Tmn_p2]; 

                FIDs_p1 = [FIDs_p1; Tmn_p1.getFIDphase(times, 0)];
                FIDs_p2 = [FIDs_p2; Tmn_p2.getFIDphase(times, 0)]; 
                
                TQs_p1((j-1)*lengthAlpha + alphaIndex ) = TQ_p1;
                TQs_p2((j-1)*lengthAlpha + alphaIndex ) = TQ_p2;   
                
                SQs_p1((j-1)*lengthAlpha + alphaIndex ) = SQ_p1;
                SQs_p2((j-1)*lengthAlpha + alphaIndex ) = SQ_p2;   
                
                [Tmn_p1, Tmn_timeSeries_TR] =  Tmn_p1.relaxation_shift(obj.TR); 
                [Tmn_p2, Tmn_timeSeries_TR] =  Tmn_p2.relaxation_shift(obj.TR);                
            end   
        end
    end   
  
end
end