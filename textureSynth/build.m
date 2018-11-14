% Example 1: Synthesis of a "text" texture image, using
% Portilla-Simoncelli texture analysis/synthesis code, based on
% alternate projections onto statistical constraints in a complex
% overcomplete wavelet representation.
%
% See Readme.txt, and headers of textureAnalysis.m and
% textureSynthesis.m for more details.
%
% Javier Portilla (javier@decsai.ugr.es).  March, 2001

close all

im0 = imread('original2_sm.bmp');	% im0 is a double float matrix!
im0 = double(im0);
im0 = power(im0/255, 2);

im0 = im0 * 255;

im0r = im0(:,:,1);
%im0g = im0(:,:,2);
%im0b = im0(:,:,3);
%im1 = pgmRead('test2.pgm');	% im0 is a double float matrix!

Nsc = 4; % Number of scales
Nor = 4; % Number of orientations
Na = 9;  % Spatial neighborhood is Na x Na coefficients
	 % It must be an odd number!
paramsr = textureAnalysis(im0r, Nsc, Nor, Na);
%paramsg = textureAnalysis(im0g, Nsc, Nor, Na);
%paramsb = textureAnalysis(im0b, Nsc, Nor, Na);



Niter = 25;	% Number of iterations of synthesis loop
Nsx = 256;	% Size of synthetic image is Nsy x Nsx
Nsy = 256;	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)

resr = textureSynthesis(paramsr, [Nsy Nsx], Niter);
%resg = textureSynthesis(paramsg, [Nsy Nsx], Niter);
%resb = textureSynthesis(paramsb, [Nsy Nsx], Niter);

res = zeros([size(resr,1), size(resr,2), 3]);
for i=1:1:size(resr,1)
    for j=1:1:size(resr,2)
        res(i,j,1) = resr(i,j);
        res(i,j,2) = resr(i,j);
        res(i,j,3) = resr(i,j);
    end
end
res = power(double(res)/255, 1/2);
res = res * 255;

close all
figure(1)
showIm(im0, 'auto', 1, 'Original texture');
figure(2)
showIm(res, 'auto', 1, 'Synthesized texture');

% Can you read the NEW text? ;-)
