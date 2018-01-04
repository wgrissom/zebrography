function d = KL(I,J,nIter,smoothWinSize,display)

if nargin < 5
	display = 0;
end

if nargin < 4
	smoothWinSize = 10;
end

if nargin < 3
	nIter = 20;
end

if nargin < 2
	error('must input two signals');
end
xxx = 50;

win = hamming(smoothWinSize); win = win(:);
% 
Is =I;% filter2(win,I(:),'valid');
Js =J;% filter2(win,J(:),'valid');

Jsx = Js([2:end,end])/2 - Js([1,1:end-1])/2;

Js = crop(Js,xxx);
Jsx = crop(Jsx,xxx);
Jsx2 = Jsx(:)'*Jsx(:);

d = 0;

for n=1:nIter

	Isc = crop(subShift(Is,d),xxx);
	delta_d = real(Jsx(:)'*(Js(:) - Isc(:)) /Jsx2);
	d = d - delta_d;
%     hold on; figure(display), plot([Js(:),Isc(:)]); pause(0.0005);
end

% plot(Is,'linewidth',1);
% hold on; plot(JS,'linewidth');
% hold on; plot(Isc,'linewidth',1);

function x = subShift(x,d)
	n = 1:length(x);
	nd = n - d;
	x = interp1(n,x,nd,'linear');
	x(find(isnan(x))) = 0;

function x = crop(x,N)
%function crops N samples from either sides of the signal
x = x(N+1:end-N);


