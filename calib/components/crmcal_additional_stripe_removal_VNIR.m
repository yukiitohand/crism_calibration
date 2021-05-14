function [RDmm,coeff] = crmcal_additional_stripe_removal_VNIR(RDm,DMdata,rownum_table)
% [RDmm,coeff] = crmcal_additional_stripe_removal_VNIR(RDm,DMdata,rownum_table)
%  Perform additional destriping of the image. The destriping is performed
%  by multiplying a factor estimated from the image itself
%  INPUTS
%    RDm: radiance image [L x S x B]
%  OUTPUTS
%    RDmm: stripe removed radiance image [L x S x B]
%    coeff: multiplication factor for destrping [1 x S x B]
%
[L,S,B] = size(RDm);
[~,lam0] = ismember([187],rownum_table);
s_w = ones(1,S);
% assuming scattered mask is same for all the wavelength bands
% scatMask = (DMdata.img(:,:,lam0(1)) == 4);
% s_w(squeeze(scatMask)) = 3;
% scene mask
sceneMask = DMdata.img(:,:,lam0(1))==1;
s_w(squeeze(sceneMask)) = 9;

imgmean = nanmean(RDm,1);
imgmean_smooth = nan(size(imgmean));
for b=1:B
    [imgmean_smooth(:,:,b)] = crmcal_movemean_robust1d_batchsort(imgmean(:,:,b)','Window_Size',s_w);
end
coeff = imgmean ./ imgmean_smooth;
RDmm = RDm./ coeff;

end
