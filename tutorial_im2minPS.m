clear;

% This is tutorial of im2minPS
% G Okazawa (2015)

%% load image

% The input image should satisfy following requirements:
% 1) should be 128 x 128 pixels grayscale
% 2) should have a unit of luminance (cd/m2)
% 3) the mean and SD of luminance should be 15 and 6cd/m2
%
% These are needed to meet the standard in Okazawa et al. (2015).
S = load('demo_img.mat');
Limg = S.Limg;

figure;
imagesc(Limg);
axis square;
set(gca, 'visible', 'off');
h = colorbar;
ylabel(h, 'Luminance (cd/m2)');
colormap gray;

%% run im2minPS
[minPS, minPSraw, minPSstruct] = im2minPS(Limg);

% minPS has a vector of 29 parameters.
% for the definition of each parameter, see minPS_list.pdf
% The values are normalized with respect to the dataset in Okazawa et al.
% (2015)

% minPSraw contains unnormalized data, but note that, during the
% calculation of minPS, there are several steps where the normalization of
% values is necessary to obtain meaningful parameters. Those normalization
% is performed using the dataset in Okazawa et al, so the parameters cannot
% be completely immuned to the dependency on the dataset. Moreover, the
% dimensionality of linear cross position and energy cross position was
% reduced using PCA, which is also dependent on the dataset. Thus, it is
% impossible to generate a "raw" value of minPS parameters.

% minPSstruct stores the same data (minPSraw) in a struct format, which
% is directly related to the output of textureAnalysis in Portilla &
% Simoncelli's code.





