function [ RT14i,SL ] = crmcal_straylight_correction( RT14h,DMdata,rownum_table,varargin )
% [ RT14i,SL ] = crmcal_straylight_correction( RT14h,DMdata,rownum_table,varargin )
%   stray light (scattered light) correction for VNIR correction
%  INPUTS
%   RT14h: [LxSxB], input image
%   DMdata: CDR DM data (detector mask)
%   rownum_table: detector ROWNUM table
%   binnig_id: binning identifier
%  OUTPUTS
%   RT14i: processed image, same size as RT14g

binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = crism_get_binx(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});   
        end
    end
end

switch binx
    case 1
        c2d = 541:620;
    case 2
        c2d = 271:310;
    case 5
        c2d = 109:124;
    case 10
        c2d = 55:62;
    otherwise
        error('Undefined binx %d',binx);
end

if isempty(DMdata.img), DMdata.readimg(); end
[L,S,B] = size(RT14h);
SL = nan(L,S,B);

%%
%-------------------------------------------------------------------------%
% stripe removal of the band lam0
%-------------------------------------------------------------------------%
% lambda0 == (detector row 187)
[~,lam0] = ismember([186 187 188],rownum_table);
% [~,lam0] = ismember([187],rownum_table);
if any(lam0==0), error('detector 187 is not used'); end

% s_w = zeros(1,S);
% assuming scattered mask is same for all the wavelength bands
scatMask = (DMdata.img(:,:,lam0(1)) == 4);
% s_w(squeeze(scatMask)) = 3;
% scene mask
sceneMask = DMdata.img(:,:,lam0(1))==1;
% s_w(squeeze(sceneMask)) = 3;

% imlam0_mean = nanmean(RT14h(:,:,lam0),1);
% [imlam0_mean_smooth] = movemean_robust1d_batchsort(imlam0_mean,'Window_Size',s_w);
% coeff = imlam0_mean ./ imlam0_mean_smooth;
% RT14h(:,:,lam0) = RT14h(:,:,lam0) ./ coeff;


%%
%-------------------------------------------------------------------------%
% the first wavelength zone (lambda < 563nm)
%-------------------------------------------------------------------------%
% For detector rows 185-215
[~,blt563] = ismember((185:215)',rownum_table);
blt563eq0 = (blt563==0);
% remove any wavelength doesn't match
if any(blt563eq0), blt563 = blt563(blt563eq0); end


% scatMask = (DMdata.img(:,:,lam0(1)) == 4);

RT14h_ab_mean = mean(RT14h(:,scatMask,:),2,'omitnan'); % L x 1 x B
SW = mean(RT14h(:,:,lam0),3,'omitnan') ./ mean(RT14h_ab_mean(:,:,lam0),3,'omitnan'); % L x S x 1

if L==1
    SL(:,:,blt563) = SW.*RT14h_ab_mean(:,:,blt563);
else
    % SL(:,:,blt563) = SW.*RT14h_ab_mean(:,:,blt563);
    SL(:,:,blt563) = crmcal_movemean_robust(SW.*RT14h_ab_mean(:,:,blt563),1);
end
% SL(:,:,:) = crmcal_movemean_robust(SW.*RT14h_ab_mean(:,:,:),1);

% RT14i(:,:,blt563) = RT14h(:,:,blt563) - crmcal_movemean_robust(SW.*RT14h_ab_mean(:,:,blt563),1);

%%
%-------------------------------------------------------------------------%
% the second wavelength zone (lambda > 563nm)
%-------------------------------------------------------------------------%
% For detector rows 216-291
[~,bgt563] = ismember((216:291)',rownum_table);
bgt563eq0 = (bgt563==0);
% remove any wavelength doesn't match
if any(bgt563eq0), bgt563 = bgt563(bgt563eq0); end

xS = 1:S; % (size of xS) = [1,S,1]
% assuming scene pixels are same for all the wavelength bands
xk_valid = double(sceneMask);
xk_valid(xk_valid==0) = nan;
xk = xS .* xk_valid; % (size of xk) = [1,S,1]
% xk4 = reshape(xk,[1,1,1,S]); % (size of xk4) = [1,1,1,S]
% 
% CB_coeff = 1 ./ (1+((xS-xk4)./10).^2); % (size of CB_coeff) = [1,S,1,S]
% Bx = 1 ./ nansum( 1 ./ (1+((xk-xk4).*binx./10).^2) ,4); % (size of Bx) = [1,S,1]
% CB = nan(L,S,B);
% for l=1:L
%     CBlk = reshape(permute(RT14h(l,:,:),[1,3,2]),[1,1,B,S]);
%     CB(l,:,:) = Bx .* nansum( CBlk .* CB_coeff , 4);
% end

CB_coeff = 1 ./ (1+((xS-xk').*binx./10).^2);
Bx = 1./sum(CB_coeff,1,'omitnan');
CB_coeff_nrmed = CB_coeff.*Bx;
CB = nan(L,S,B);
for l=1:L
%     CBlk = reshape(permute(RT14h(l,:,:),[1,3,2]),[1,1,B,S]);
%     CB(l,:,:) = Bx .* nansum( CBlk .* CB_coeff , 4);
%     CBlk = squeeze(RT14h(l,:,:)); % S x B
%     for b=1:B
%         CBlb = sum(CB_coeff_nrmed .* CBlk(:,b),1,'omitnan'); % (SxS) .* (Sx1)
%         CB(l,:,b) = CBlb;
%     end
    CBlk = permute(RT14h(l,:,:),[2,1,3]); % S x 1 x B
    CB(l,:,:) = sum(CB_coeff_nrmed .* CBlk,1,'omitnan');
end

% SW_cd_mean = nanmean(SW(:,c2d,:),2); % L x 1 x B
RT14h_cd_mean = mean(RT14h(:,c2d,:),2,'omitnan'); % L x 1 x B
CS = SW .* (CB./mean(CB(:,:,lam0),3,'omitnan')) .* (mean(RT14h_cd_mean(:,:,lam0),3,'omitnan')./RT14h_cd_mean);

if L==1
    SL(:,:,bgt563) = CS(:,:,bgt563).*RT14h_ab_mean(:,:,bgt563);
else
    % SL(:,:,bgt563) = CS(:,:,bgt563).*RT14h_ab_mean(:,:,bgt563);
    SL(:,:,bgt563) = crmcal_movemean_robust(CS(:,:,bgt563).*RT14h_ab_mean(:,:,bgt563),1);
end
% SL = crmcal_movemean_robust(CS(:,:,:).*RT14h_ab_mean(:,:,:),1);

%%
RT14i = RT14h - SL;


end


