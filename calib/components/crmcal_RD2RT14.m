function [RT14j] = crmcal_RD2RT14(RDm,RSPl,NUdata)
% [RT14j] = crmcal_RD2RT14(RDm,RSPl,NUdata)
%  reverse the scene radiance at the instrument aperture flat field applied
%  to RT14j.
%   Input Parameters
%     RDn: scene image radiance 
%     RSPl : binned spectral radiometric responsitivity
%     NUdata: CRISMdata obj, CDR NU data
%   Output Parameters
%     RT14j  : scene [counts/ms]
%     FF   : flat fielding component
%   *Detail*
%    see Section 2.15 of "CRISM_DPSIS_Appendix_L_v5_2016-09-01.pdf"


if isempty(NUdata.img), NUdata.readimg(); end;

[L,S,Bands] = size(RDm);
[L1,S1,Bands1] = size(RSPl);

if (S~=S1) || (Bands~=Bands1)
    error('Size of RTj and RSPl does not match');
end

FL = NUdata.img(1,:,:);
EP = NUdata.img(2,:,:);

if any(EP~=0)
    error('EP~=0 is not supported');
end

RT14j = RDm .* repmat(FL.*RSPl,[L,1,1]);

end