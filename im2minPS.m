function [minPS, minPSraw, minPSstruct] = im2minPS(im)
%% Convert images to minPS statistics
%
% [minPS minPSraw minPSstruct] = im2minPS(im)
%
% Input:
%       im - should be a grayscale image (128 x 128 double).
%            Or, this can be a cell including an image in each component (N x
%            1). (N = number of images)
%            The unit is supposed to be luminance (cd/m2).
%            The average and SD of luminances are supposed to be 15 and 6 cd/m2, respectively.
%
% Output:
%       minPS    - a matrix [image num x 29] or a vector [1 x 29]
%                  that contains minPS statistics (29 parameters)
%                  The values are normalized such that the image database used in Okazawa et al.(2015).
%                  have 0 mean and 1 SD.
%       minPSraw - Unnormalized version of minPS statistics
%                   (note: still the values are normalized in each
%                   calculation step to eliminate artifacts due to large
%                   difference in values across filter scales. See comments 
%                   in the code for the details).
%       minPSstruct - a struct containing minPS values
%
%
% G Okazawa (2015)
    addpath matlabPyrTools;
    addpath matlabPyrTools/MEX;
    addpath textureSynth;
    addpath textureSynth/MEX;


    Nsc = 4;
    Nor = 4;
    Na = 7;


    if iscell(im)
        for n=1:length(im)
            if ~isequal(size(im{n}), [128 128])
                error('Image size should be 128 x 128');
            end
            PS(n) = textureAnalysis_acr(im{n}, Nsc, Nor, Na); %#ok<AGROW>
        end
    else
        if ~isequal(size(im), [128 128])
            error('Image size should be 128 x 128');
        end
        PS = textureAnalysis_acr(im, Nsc, Nor, Na);
    end

    d = load('parameters_Okazawa2015.mat');
    param = d.param;

    PS = reduce_pixelStats(PS);
    PS = reduce_magMeans(PS, param);
    PS = reduce_autoCorrReal(PS, param);
    PS = reduce_parentRealCorr(PS, param);
    PS = reduce_cousinMagCorr(PS, param);
    PS = reduce_autoCorrMag(PS, param);
    PS = reduce_parentMagCorr(PS, param);

    PS = rmfield(PS, 'pixelLPStats');
    PS = rmfield(PS, 'varianceHPR');
    PS = rmfield(PS, 'cousinRealCorr');

    minPSstruct = PS;

    minPSraw = vectorize_minPS(minPSstruct);

    minPS = normalize_values(minPSraw, param);

end


%%  reduce_pixelStats 
% From marginal statistics, extract skewness only
function PS = reduce_pixelStats(PS)
    for n=1:length(PS)
        PS(n).pixelStats = PS(n).pixelStats(3); % skew only
    end
end

%%  reduce_magMeans 
% extract spectral statistics.
% reduced to 2 scale and 2 orientation
function [PS, param] = reduce_magMeans(PS, param)
    magM = [PS.magMeans]';
    magM = magM(:, 2:end-1); % remove high and low pass filter

    Nimg = size(magM,1);

    % normalization of magnitude of each subband
    % this is necessary because, without normalization, averaging across scales and
    % orientations generate results dominated by low frequency scales
    % (because they have large amplitudes).
    if ~exist('param', 'var')
        [magM, param.magMeans_mn, param.magMeans_sd] = zscore(magM);
    else
        magM = (magM - (ones(Nimg,1) * param.magMeans_mn)) ./ (ones(Nimg,1) * param.magMeans_sd); % normalize
    end

    magM = reshape(magM, Nimg, 4, 4); % Nimg x ori x scale


    newMagM = zeros(Nimg, 4);
    newMagM(:,1) = magM(:,1,1) + magM(:,1,2) + (magM(:,2,1) + magM(:,2,2)) * .5 + (magM(:,4,1) + magM(:,4,2)) * .5; % fine vertical
    newMagM(:,2) = magM(:,3,1) + magM(:,3,2) + (magM(:,2,1) + magM(:,2,2)) * .5 + (magM(:,4,1) + magM(:,4,2)) * .5; % fine horizontal
    newMagM(:,3) = magM(:,1,3) + magM(:,1,4) + (magM(:,2,3) + magM(:,2,4)) * .5 + (magM(:,4,3) + magM(:,4,4)) * .5; % coarse vertical
    newMagM(:,4) = magM(:,3,3) + magM(:,3,4) + (magM(:,2,3) + magM(:,2,4)) * .5 + (magM(:,4,3) + magM(:,4,4)) * .5; % coarse horizontal

    for n=1:Nimg
        PS(n).magMeans = newMagM(n,:);
    end
end

%%  reduce_autoCorrReal 
% extract linear cross position parameters
function [PS, param] = reduce_autoCorrReal(PS, param)

    param_mode = exist('param', 'var');

    Nimg = length(PS);

    % concatenate
    acrmat = zeros(Nimg, 7*7, 4,4); % Nimg x pixel(49) x 4 scale x 4 ori
    for n=1:Nimg
        acr = PS(n).autoCorrReal;
        acrmat(n,:,:,:) = reshape(acr, 7*7, 4, 4);
    end

    % normalize and average
    if ~param_mode
        MN = zeros(Nimg, 4); % Nimg x 4 scale
        SD = zeros(Nimg, 4);
        for n=1:Nimg
            acr = PS(n).autoCorrReal;
            for m=1:4
                a = acr(:,:,m,:);
                MN(n,m) = mean(a(:));
                SD(n,m) = std(a(:));
            end
        end
        MN = mean(MN,1);
        SD = mean(SD,1);
        param.autoCorrReal_mn = MN;
        param.autoCorrReal_sd = SD;
    else
        MN = param.autoCorrReal_mn;
        SD = param.autoCorrReal_sd;
    end
    for n=1:4
        acrmat(:,:,n,:) = (acrmat(:,:,n,:) - MN(n))/SD(n);
            % normalize each scale because different scales have very
            % different mean.
    end

    acrmat = permute(acrmat, [1 3 4 2]); % Nimg x 4 scale x 4 ori x 49
    acrmat = reshape(acrmat, Nimg * 4 * 4, 49);

    % the mean should be 0 before running pca
    if ~param_mode
        param.autoCorrReal_acrmat_mn = mean(acrmat,1);
    end
    acrmat = acrmat - ones(size(acrmat,1),1) * param.autoCorrReal_acrmat_mn;

    % PCA
    if ~param_mode
        [acr_coef, acr_score, acr_latent] = pca(acrmat);
        ndim_ac = 4;
        acr_score = acr_score(:, 1:ndim_ac); % (Nimg*4*4) x ndim_ac
        param.autoCorrReal_coef = acr_coef;
        param.autoCorrReal_latent = acr_latent;
        param.autoCorrReal_ndim = ndim_ac;
    else
        ndim_ac = param.autoCorrReal_ndim;
        acr_score = acrmat * param.autoCorrReal_coef;
        acr_score = acr_score(:, 1:ndim_ac); % (Nimg*4*4) x ndim_ac
    end

    % averaging
    acr_score = reshape(acr_score, Nimg, 4, 4, ndim_ac); % Nimg x 4scale x 4ori x N pc
    acr_score = mean(mean(acr_score, 2),3); % Nimg x 1 x 1 x N pc
    acr_score = reshape(acr_score, [Nimg ndim_ac]);

    % back to param
    for n=1:Nimg
        PS(n).autoCorrReal = acr_score(n,:);
    end
end


%%  reduce_parentRealCorr 
% extract linear cross scale 
function [PS, param] = reduce_parentRealCorr(PS, param)

    Nimg = length(PS);

    % concatenate
    parentCorr = zeros(Nimg, 4, 4, 3); % Nimg x 4ori x 4ori x 3scale
    MN = zeros(Nimg, 3);
    SD = zeros(Nimg, 3);
    for n=1:Nimg
        a = PS(n).parentRealCorr(1:4,1:4,1:3); % discard 'imaginary correlation'
        parentCorr(n,:,:,:) = a;
        a = reshape(a, 16,3);
        MN(n,:) = mean(a,1);
        SD(n,:) = std(a,1,1);
    end
    % normalize
    if ~exist('param', 'var')
        MN = mean(MN,1);
        SD = mean(SD,1);
        param.parentRealCorr_mn = MN;
        param.parentRealCorr_sd = SD;
    else
        MN = param.parentRealCorr_mn;
        SD = param.parentRealCorr_sd;
    end

    for n=1:3
        parentCorr(:,:,:,n) = (parentCorr(:,:,:,n) - MN(n))/SD(n);
            % normalize by the mean of each scale
    end

    % average scale
    parentCorr2 = zeros(Nimg, 4, 4, 2); % 4ori x 4ori x 2scale
    parentCorr2(:,:,:,1) = parentCorr(:,:,:,1) + parentCorr(:,:,:,2) * .5;
    parentCorr2(:,:,:,2) = parentCorr(:,:,:,3) + parentCorr(:,:,:,2) * .5;


    % average ori
    parentCorr3 = zeros(Nimg, 2, 2); % 2 ori x 2 scale

    for n=1:2 % extract correlation between filter responses with the same orientation
        parentCorr3(:, 1, n) = parentCorr2(:, 1, 1, n) + parentCorr2(:, 2, 2, n) * .5 + parentCorr2(:, 4, 4, n) * .5;
        parentCorr3(:, 2, n) = parentCorr2(:, 3, 3, n) + parentCorr2(:, 2, 2, n) * .5 + parentCorr2(:, 4, 4, n) * .5;
    end


    % back
    for n=1:Nimg
        PS(n).parentRealCorr = squeeze(parentCorr3(n,:,:));
    end
end

%%  reduce_cousinMagCorr 
% extract energy cross orientation
function [PS, param] = reduce_cousinMagCorr(PS, param)

    Nimg = length(PS);

    % concatenate
    cousinCorr = zeros(Nimg, 6, 4); % Nimg x 6combination x 4scale
    for n=1:Nimg
        a = PS(n).cousinMagCorr(:,:,1:4);
        a = reshape(a, 16,4);
        a = a([2;3;4;7;8;12;], :);
            % [2;3;4;7;8;12;] are non-diagonal combinations
            % 2.. A x B, 3.. A x C, 4.. A x D
            % 7.. B x C, 8.. B x D, 12.. C x D
            % (A.. vertical, B.. right diagonal, 
            %    C.. horizontal, D.. left diagonal)
        cousinCorr(n,:,:) = a;
    end
    % normalize
    if ~exist('param', 'var')
        MN = zeros(Nimg, 4);
        SD = zeros(Nimg, 4);
        for n=1:Nimg
            a = squeeze(cousinCorr(n,:,:));
            MN(n,:) = mean(a,1);
            SD(n,:) = std(a,1,1);
        end
        MN = mean(MN,1);
        SD = mean(SD,1);
        param.cousinMagCorr_mn = MN;
        param.cousinMagCorr_sd = SD;
    else
        MN = param.cousinMagCorr_mn;
        SD = param.cousinMagCorr_sd;
    end

    for n=1:3
        cousinCorr(:,:,n) = (cousinCorr(:,:,n) - MN(n))/SD(n); % normalize each scale
    end

    % average scale
    cousinCorr2 = zeros(Nimg, 6, 2); % Nimg x 6combination x 2scale
    cousinCorr2(:,:,1) = cousinCorr(:,:,1) + cousinCorr(:,:,2); % fine
    cousinCorr2(:,:,2) = cousinCorr(:,:,3) + cousinCorr(:,:,4); % coarse


    % average ori
    cousinCorr3 = zeros(Nimg, 3, 2); % Nimg x 3combination x 2scale

    for n=1:2
        cousinCorr3(:, 1, n) = cousinCorr2(:, 1, n) + cousinCorr2(:, 3, n); % A-B, A-D vertical vs. oblique
        cousinCorr3(:, 2, n) = cousinCorr2(:, 2, n) + cousinCorr2(:, 5, n); % A-C, B-D vertical vs. horizontal and R oblique vs. L oblique
        cousinCorr3(:, 3, n) = cousinCorr2(:, 4, n) + cousinCorr2(:, 6, n); % B-C, C-D horizontal vs. oblique
            % (A.. vertical, B.. right diagonal, 
            %    C.. horizontal, D.. left diagonal)
    end


    % back
    for n=1:Nimg
        PS(n).cousinMagCorr = squeeze(cousinCorr3(n,:,:));
    end
end


%%  reduce_autoCorrMag 
% extract energy cross position
function [PS, param] = reduce_autoCorrMag(PS, param)

    param_mode = exist('param', 'var');

    Nimg = length(PS);

    % concatenate
    acmat = zeros(Nimg, 7*7, 4,4);
    for n=1:Nimg
        ac = PS(n).autoCorrMag;
        acmat(n,:,:,:) = reshape(ac, 7*7, 4, 4);
    end

    % normalize and average
    if ~param_mode
        MN = zeros(Nimg, 4);
        SD = zeros(Nimg, 4);
        for n=1:Nimg
            ac = PS(n).autoCorrMag;
            for m=1:4
                a = ac(:,:,m,:);
                MN(n,m) = mean(a(:));
                SD(n,m) = std(a(:));
            end
        end
        MN = mean(MN,1);
        SD = mean(SD,1);
        param.autoCorrMag_mn = MN;
        param.autoCorrMag_sd = SD;
    else
        MN = param.autoCorrMag_mn;
        SD = param.autoCorrMag_sd;
    end
    for n=1:4
        acmat(:,:,n,:) = (acmat(:,:,n,:) - MN(n))/SD(n); % normalize each scale
    end


    acmat = permute(acmat, [1 3 4 2]); % Nimg x 4 scale x 4 ori x 49
    acmat = reshape(acmat, Nimg * 4 * 4, 49);

    % the mean should be 0 before running pca
    if ~param_mode
        param.autoCorrMag_acmat_mn = mean(acmat,1);
    end
    acmat = acmat - ones(size(acmat,1),1) * param.autoCorrMag_acmat_mn;

    % PCA
    if ~param_mode
        [ac_coef, ac_score, ac_latent] = pca(acmat);
        ndim_ac = 3;
        ac_score = ac_score(:, 1:ndim_ac); % Nimg*4 x ndim_ac
        param.autoCorrMag_coef = ac_coef;
        param.autoCorrMag_latent = ac_latent;
        param.autoCorrMag_ndim = ndim_ac;
    else
        ndim_ac = param.autoCorrMag_ndim;
        ac_score = acmat * param.autoCorrMag_coef;
        ac_score = ac_score(:, 1:ndim_ac);
    end

    % averaging
    ac_score = reshape(ac_score, Nimg, 4, 4, ndim_ac); % Nimg x 4scale x 4ori x N pc
    ac_score = reshape(mean(ac_score, 2), [Nimg 4 ndim_ac]); % averaged across scale: Nimg x 1 x 4ori x N pc
    ac_score2 = zeros(Nimg, 2, ndim_ac); % Nimg x 2ori x 3 pc
    
        % averaging across ori
    for n=1:ndim_ac
        ac_score2(:,1,n) = ac_score(:,1,n) + ac_score(:,2,n) * .5 + ac_score(:,4,n) * .5;
        ac_score2(:,2,n) = ac_score(:,3,n) + ac_score(:,2,n) * .5 + ac_score(:,4,n) * .5;
    end

    % back to param

    for n=1:Nimg
        ac = squeeze(ac_score2(n, :,:));
        PS(n).autoCorrMag = ac;
    end
end

%%  reduce_parentMagCorr 
% extract energy cross scale
function [PS, param] = reduce_parentMagCorr(PS, param)

    Nimg = length(PS);

    % concatenate
    parentCorr = zeros(Nimg, 4, 4, 3); % Nimg x 4ori x 4ori x 3scale
    for n=1:Nimg
        a = PS(n).parentMagCorr(1:4,1:4,1:3);
        parentCorr(n,:,:,:) = a;
    end
    
    % normalize
    if ~exist('param', 'var')
        MN = zeros(Nimg, 3);
        SD = zeros(Nimg, 3);
        for n=1:Nimg
            a = parentCorr(n,:,:,:);
            a = reshape(a, 16,3);
            MN(n,:) = mean(a,1);
            SD(n,:) = std(a,1,1);
        end
        MN = mean(MN,1);
        SD = mean(SD,1);
        param.parentMagCorr_mn = MN;
        param.parentMagCorr_sd = SD;
    else
        MN = param.parentMagCorr_mn;
        SD = param.parentMagCorr_sd;
    end

    for n=1:3
        parentCorr(:,:,:,n) = (parentCorr(:,:,:,n) - MN(n))/SD(n); % normalize each scale
    end

    % average scale
    parentCorr2 = zeros(Nimg, 4, 4, 2); % Nimg x 4ori x 4ori x 2scale
    parentCorr2(:,:,:,1) = parentCorr(:,:,:,1) + parentCorr(:,:,:,2) * .5;
    parentCorr2(:,:,:,2) = parentCorr(:,:,:,3) + parentCorr(:,:,:,2) * .5;


    % extract cross correlation of the same orientation filter across
    % scales
    parentCorr3 = zeros(Nimg, 2, 2); % Nimg x 2ori x 2scale

    for n=1:2
        parentCorr3(:,1,n) = parentCorr2(:, 1, 1, n) + parentCorr2(:, 2, 2, n) * .5 + parentCorr2(:, 4, 4, n) * .5; % vertical
        parentCorr3(:,2,n) = parentCorr2(:, 3, 3, n) + parentCorr2(:, 2, 2, n) * .5 + parentCorr2(:, 4, 4, n) * .5; % horizontal
    end


    % back
    for n=1:Nimg
        PS(n).parentMagCorr = squeeze(parentCorr3(n,:,:));
    end
end

%%  vectorize_minPS 
function vector = vectorize_minPS(PS)
        
    vector = zeros(length(PS), 29);

    for n=1:length(PS)
        vector(n,:) = [PS(n).magMeans(:)', ...
                        PS(n).pixelStats(:)', ...
                        PS(n).autoCorrReal(:)', ...
                        PS(n).parentRealCorr(:)', ...
                        PS(n).autoCorrMag(:)', ...
                        PS(n).cousinMagCorr(:)', ...
                        PS(n).parentMagCorr(:)'];
    end
end

%%  normalize_values   
function [minPS, param] = normalize_values(minPSraw, param)
    if ~exist('param', 'var')
        param.Nmn = mean(minPSraw,1);
        param.Nsd = std(minPSraw, 1,1);
    end
    minPS = minPSraw - ones(size(minPSraw,1),1) * param.Nmn;
    minPS = minPS ./ (ones(size(minPS,1),1) * param.Nsd);
end







