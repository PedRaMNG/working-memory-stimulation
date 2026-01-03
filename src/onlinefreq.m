function [dFreq, smoOoth] = onlinefreq(sign, winwid, steps, sampletime, pOl, maxAmp)
    % (signal, window_width, each_step_shifted_window, sampling_time,
    % smoother_pole, maximum_amplitude);
    
    sign    = sign(:);
    n       = numel(sign);
    dFreq   = zeros(n, 1);
    window  = ones(winwid,1);
    h2      = sampletime*winwid;
    nimwid  = floor(winwid/2);
    
    for t = 1:steps:(n-winwid)
        h1    = sign(t:(t+winwid-1)).*window;
        Fre   = sum(h1 == maxAmp)/h2;
        dFreq(t:(t + steps - 1)) = repmat( Fre , 1, steps);
        
        %% FFT proccess
            % fr       = fft(h1);
            % fr       = fr(1:nimwin);
            % amp      = abs(fr);
            % freqs    = linspace(0, pi, nimwin) / (2*pi*sampletime);
            % dFreq(t:(t + steps - 1)) = repmat(freqs(amp == max(amp)), 1, steps);
        
        % % Plots -----------------------------------------------------------
        % subplot(311);
        %     plot(h1);
        %     xlabel('Time');
        %     ylabel('Amp');
        % subplot(212);
        %     semilogx(freqs, amp);
        %     grid on;
        %     xlabel('Frequency (Hz)');
        %     ylabel('Amplitude');
        % 
        % pause(0.001);
        % % -----------------------------------------------------------------
    end
    g        = tf(1, [pOl, 1]);
    timeline = linspace(0, numel(dFreq)*sampletime, numel(dFreq));
    smoOoth  = lsim(g, dFreq, timeline);
end

