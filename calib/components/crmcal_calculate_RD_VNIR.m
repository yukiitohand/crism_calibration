function [RDm,FF] = crmcal_calculate_RD_VNIR(RTj,RSPl,SPdata,SSdata,NUdata,rowNumTableRSPj)
% [RDm,FF] = crmcal_calculate_RD_VNIR(RTj,RSPl,SPdata,SSdata,NUdata,rowNumTableRSPj)
%  calculate scene radiance at the instrument aperture and divide by flat
%  field for VNIR.
%   Input Parameters
%     RTj: scene image [counts/ms]
%     RSPl : binned spectral radiometric responsitivity
%     SPdata: CRISMdata obj, CDR SP data
%     SSdata: CRISMdata obj, CDR SS data
%     NUdata: CRISMdata obj, CDR NU data
%     rowNumTableRSPj: rownumtable for RSP
%   Output Parameters
%     RDm  : calculated scene radiance
%     FF   : flat fielding component
%   *Detail*
%    see Section 2.15 of "CRISM_DPSIS_Appendix_L_v5_2016-09-01.pdf"


if isempty(NUdata.img), NUdata.readimg(); end
if isempty(SPdata.img), SPdata.readimg(); end
if isempty(SSdata.img), SSdata.readimg(); end

% propNU7 = NUdata.prop;
% propNU7.version = 3;
% propNU7.version = 7;
% basenameNU7 = crism_get_basenameCDR4_fromProp(propNU7);
% NUdata7 = CRISMdata(basenameNU7,'');
% NUdata7.readimg();

[L,S,Bands] = size(RTj);
[L1,S1,Bands1] = size(RSPl);

if (S~=S1) || (Bands~=Bands1)
    error('Size of RTj and RSPl does not match');
end

% For detector rows 185-215, SSdata is only applied
[~,blt563] = ismember((185:215)',rowNumTableRSPj);
blt563eq0 = (blt563==0);
% remove any wavelength doesn't match
if any(blt563eq0), blt563 = blt563(blt563eq0); end

FL = NUdata.img(1,:,:);
EP = NUdata.img(2,:,:);
nu1 = FL;
% nu1mean = nanmean(nu1,2); % level adjustement
FF = ones(size(RTj));
% for l=1:L
%      FF(l,:,blt563) = FL(:,:,blt563) + EP(:,:,blt563) .* log(RTj(l,:,blt563) ./ RSPl(:,:,blt563)) ;
% end

FF(:,:,blt563) = FL(:,:,blt563) + EP(:,:,blt563) .* log(RTj(:,:,blt563) ./ RSPl(:,:,blt563)) ;

% the flat fielding introduces additional stripes...
% for l=1:L
%      FF= FL + EP.* log(RTj(l,:,:) ./ RSPl) ;
% end

% it is tricky for % For detector rows 216-291
[~,bgt563] = ismember((216:291)',rowNumTableRSPj);
bgt563eq0 = (bgt563==0);
% remove any wavelength doesn't match
if any(bgt563eq0), bgt563 = bgt563(bgt563eq0); end
nu1 = FL;
nu1mean = mean(nu1,2,'omitnan'); % level adjustement
K = SPdata.img./SSdata.img(1,:,:); 
coeff = nu1./K;
coeff_mean = mean(coeff,2,'omitnan');
coeff_nrmed = coeff ./ coeff_mean; % stripe removal
% for l=1:L
%     FF(l,:,bgt563) = nu1mean(:,:,bgt563) .* coeff_nrmed(:,:,bgt563);
% end 
FF(:,:,bgt563) = nu1mean(:,:,bgt563) .* coeff_nrmed(:,:,bgt563);



%FF = repmat(FL,[L,1,1]) + repmat(EP,[L,1,1]) .* log( RTj ./ repmat(RSPl,[L,1,1]) );

% RDm = RTj ./ (FF .* repmat(RSPl,[L,1,1]));

RDm = RTj ./ (FF .* RSPl);

end
