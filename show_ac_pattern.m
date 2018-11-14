clear;


S = load('parameters');
param = S.param;

figurea4;

for n=1:param.autoCorrReal_ndim
    subplot(4,2, n);
    imagesc(reshape(param.autoCorrReal_coef(:,n), 7,7));
    axis square;
    set(gca, 'XTick', [], 'YTick', []);
    title(sprintf('Linear cross position: PC %d', n));
end

for n=1:param.autoCorrMag_ndim
    subplot(4,2, n + 4);
    imagesc(reshape(param.autoCorrMag_coef(:,n), 7,7));
    axis square;
    set(gca, 'XTick', [], 'YTick', []);
    title(sprintf('Energy cross position: PC %d', n));
end

