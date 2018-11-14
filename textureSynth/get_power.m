function pow  = get_phase(img, param)
% function [oimg ppyr2] = randomize_phase_real(img, param)
% 
% param.nsc = [5]  number of scale
% param.ppyr2 = seed matrix

defaults.nsc = 5; % number of scales;
defaults.addphase = 0; % whether add or replace phase.
param = propval(param, defaults);
if isfield(param, 'ppyr2')
    ppyr2 = param.ppyr2;
end

if max(img(:)) > 1
    img = double(img)/255;
end
img = power(img, 2.2);

oimg = zeros(size(img));
for n=1:size(img,3)
    %%
    Nsc = param.nsc; % Number of scales
    Nor = 4; % Number of orientations
    [pyr0,pind0] = buildSCFpyr(img(:,:,n),Nsc,Nor-1);
    pow = abs(pyr0);
	return;
end


