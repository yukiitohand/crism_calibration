function [RT14h2,dc] = crmcal_dark_column_subtract(RT14h,DMdata)
% [RT14h2,dc] = crmcal_dark_column_subtract(RT14h,DMdata)
%  subtract median of dark columns. The pixel of the dark columns are
%  specified by detector mask (CRISM DM data)
%   Input Parameters
%    RT14h: image data[counts/ms]
%    DMdata: CRISMdata obj, CDR DM data
%   Output Parameters
%    RT14h2: image dark column subtracted
%    dc: median of dark column

if isempty(DMdata.img), DMdata.readimg(); end

[L,S,B] = size(RT14h);
% darkMask = double(DMdata.img == 2);
% darkMask(darkMask==0) = nan;
% imgDark = RT14h;
% % for l=1:L
% %     imgDark(l,:,:) = imgDark(l,:,:) .* darkMask;
% % end
% 
% imgDark = imgDark .* darkMask;
% % imgDark = reshape(permute(imgDark,[2,3,1]),[S*B,L]);
% % dc = median(imgDark,1,'omitnan');
% % dc = squeeze(dc)';
% % RT14h2 = RT14h - repmat(dc,[1,S,B]);

darkMask = (DMdata.img == 2);
% Dark columns are assumed to be same for all bands.
darkMask_1d = darkMask(1,:,round(B/2));
imgDark = RT14h(:,darkMask_1d,:);

dc = median(imgDark,[2,3],'omitnan');
RT14h2 = RT14h - dc;

end