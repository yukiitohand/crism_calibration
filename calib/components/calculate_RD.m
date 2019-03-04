function [RDm,FF] = calculate_RD(RTj,RSPl,NUdata)
% [RD] = calculate_RD(RTj,RSPl,NUdata)
%  calculate scene radiance at the instrument aperture and divide by flat
%  field
%   Input Parameters
%     RTj: scene image [counts/ms]
%     RSPl : binned spectral radiometric responsitivity
%     NUdata: CRISMdata obj, CDR NU data
%   Output Parameters
%     RDm  : calculated scene radiance
%     FF   : flat fielding component
%   *Detail*
%    see Section 2.15 of "CRISM_DPSIS_Appendix_L_v5_2016-09-01.pdf"


if isempty(NUdata.img), NUdata.readimg(); end;

[L,S,Bands] = size(RTj);
[L1,S1,Bands1] = size(RSPl);

if (S~=S1) || (Bands~=Bands1)
    error('Size of RTj and RSPl does not match');
end

FL = NUdata.img(1,:,:);
EP = NUdata.img(2,:,:);

FF = zeros(size(RTj));
for l=1:L
    FF(l,:,:) = FL + EP .* log(RTj(l,:,:) ./ RSPl) ;
end

%FF = repmat(FL,[L,1,1]) + repmat(EP,[L,1,1]) .* log( RTj ./ repmat(RSPl,[L,1,1]) );

RDm = RTj ./ (FF .* repmat(RSPl,[L,1,1]));

end
