% spectrum_Baz.m
%
% Author: Luis Miguel Bazdresch Sierra
%         Departamento de Electronica, Sistemas e Informática
%         Universidad ITESO
%Modified by: Omar Longoria  2012, for  Matlab version
%
% Finds the spectrum of a signal, using the Fast Fourier Transform. Originally
% based on 'plotspec.m' from the book Telecommunications Breakdown, by Johnson
% and Sothares, this funcion allows control over many of the fft and plotting
% parameters, and adds phase spectrum plotting.
%
% input arguments:
%
% x  : signal
% Ts : sampling period (real > 0)
%      0.5 : (default)
%      y : Ts = y
% N  : fft size (integer > 0)
%      0 : (default) N = 2^nextpow2(length(x))
%      y : N = y
% win : window input data (0 - 1)
%      0 : (default) no windowing
%      1 : hanning window
%      2 : kaiser window
% dp : control plotting (0 - 4)
%      0 : (default) no plotting
%      1 : plot mag spec and phase spec in two figures
%      2 : plot mag spec only
%      3 : plot phase spec only
%      4 : plot mag spec and phase spec in one figure
% lang : language of plot labels
%      0 : (default) english
%      1 : spanish
%     -1 : do not write plot labels
% clip : clip amplitudes before finding phase (0 - 100)
%      0 : (default) no clipping
%      y : clip to zero if amplitude is less than y% of max amplitude
% pctrl : which command to plot with
%      0 : (default) use stem
%      1 : use plot
%      2 : use semilogy
%
% output:
%
% ssf  : vector of frequencies
% afxs : magnitude spectrum
% pfxs : phase spectrum
%
% Spectrums may be plotted, for example, with plot(ssf,afxs)

function [ssf,afxs, pfxs] = spectrum_Baz(x, Ts, N, win, dp, lang, clip, pctrl)
	if( nargin == 0 )
		usage( 'At least one argument is needed' );
	end
	if( N == 0 )
		N = 2^nextpow2(length(x));
	end

	ssf=(-N/2:N/2-1)/(N*Ts);              % frequency vector
	if( win == 1 )
		x = hanning(length(x))'.*x;       % hanning window
    elseif (win == 2)
        x = kaiser(length(x),8)'.*x;
    end
	fx=fftshift(fft(x,N));                % FFT
	afxs=abs(fx)/length(x);               % magnitude spectrum (scaled)
	if( clip ~= 0 )                       % clip
		c = max(afxs)/clip;
		for i=1:length(afxs)
			if afxs(i) < c
				fx(i) = 0;
			end
		end
	end
	pfxs=unwrap(angle(fx));               % phase spectrum
	% plots
	if( dp == 1 || dp == 2 || dp == 4 )
		figure();
		if( dp == 4 )
			subplot(2,1,1);
		end
		if( pctrl == 0)
			stem(ssf,afxs);
		elseif( pctrl == 1 )
			plot(ssf,afxs);
		elseif( pctrl == 2 )
			semilogy(ssf,afxs);
        end
		if( lang == 0 )
			title('Magnitude spectrum');
			xlabel('Frequency (Hz)');
			ylabel('Magnitude X(f)');
		elseif( lang == 1 )
			title('Espectro de Magnitud');
			xlabel('Frecuencia (Hz)');
			ylabel('Magnitude X(f)');
		end
	end
	if( dp == 1 || dp == 3 || dp == 4)
		if( dp == 4 )
			subplot(2,1,2);
		else
			figure();
		end
		if( pctrl == 0)
			stem(ssf,pfxs);
		elseif( pctrl == 1 )
			plot(ssf,pfxs);
		elseif( pctrl == 2 )
			semilogy(ssf,pfxs);
		end
		if( lang == 0 )
			title('Phase spectrum');
			xlabel('Frequency (Hz)');
			ylabel('Phase (radians)');
		elseif( lang == 1 )
			title('Espectro de fase');
			xlabel('Frecuencia (Hz)');
			ylabel('Fase (radianes)');
		end
	end
end
