%% Copyright (c) 2012 Miguel Bazdresch
%%
%% Permission is hereby granted, free of charge, to any person obtaining a
%% copy of this software and associated documentation files (the "Software"),
%% to deal in the Software without restriction, including without limitation
%% the rights to use, copy, modify, merge, publish, distribute, sublicense,
%% and/or sell copies of the Software, and to permit persons to whom the
%% Software is furnished to do so, subject to the following conditions:
%%
%% The above copyright notice and this permission notice shall be included in
%% all copies or substantial portions of the Software.
%%
%% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
%% DEALINGS IN THE SOFTWARE.

%% usage: [P symsadj) = eyed(s, M [, Ts [, neye [, first [, syms [, doplot]]]]])
%%
%% Plot a signal's eye diagram.
%%
%% Input arguments:
%%
%% s      : signal
%% M      : number of samples in each symbol interval
%% Ts     : sampling interval (default = 0.5)
%% neye   : how many symbols to include in diagram (default = 3)
%% first  : symbol number where diagram starts (default = 50)
%% syms   : how many symbols to include in diagram (default = 100)
%% doplot : if set to zero, do not plot (default = 1)
%%
%% Output:
%%
%% P : a matrix, where each column corresponds to one line in the eye diagram.
%%     May be useful to 'export' the eye diagram for plotting with an external
%%     utility.
%% symsadj : how many symbols were actually included in the diagram
%%
%% Example: plot the eye diagram of an NRZ signal with additive gaussian noise:
%%
%% bits = round(rand(1,160));
%% mp = 10;
%% pnrz = ones(1,mp);
%% s = zeros(1,(numel(bits)-1)*mp+1);
%% s(1:mp:end) = bits;
%% xnrz = conv(s,pnrz);
%% xnrzn = xnrz + sqrt(0.01)*randn(1,numel(xnrz));
%% eyed(xnrzn,mp);

%% Author: Miguel Bazdresch

function [P symsadj] = eyed(s, M, Ts, neye, first, syms, doplot)

	if( nargin < 2)
		usage( 'At least two arguments are needed' );
	end
	if( nargin < 7 )
		doplot = 1;
	end
	if( nargin < 6 )
		syms = 100;
	end
	if( nargin < 5 )
		first = 50;
	end
	if(nargin < 4)
		neye = 3;
	end
	if(nargin < 3)
		Ts = 0.5;
	end

	%% time vector
	t = (0:neye*M-1)*Ts;

	%% if necessary, adjust syms down so that the number of symbols plotted
	%% is a multiple of neye
	symsadj = floor(syms/neye)*neye;

	x = s(first*M+1:(first+symsadj)*M); % slice x from first to first+syms
	c = symsadj/neye;                   % number of lines in diagram
	P = reshape(x,neye*M,c);            % P has one diagram line per column
	if doplot
		figure();
		plot(t,P,'k');   % without 'k', each line is a different color
		xlabel('time');
		ylabel('Amplitude');
		axis('tight');   % sometimes plot leaves white space around diagram
	end

end
