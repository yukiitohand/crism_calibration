function [RDm,FF] = crmcal_calculate_RD(RTj,RSPl,NUdata,varargin)
% [RDm,FF] = crmcal_calculate_RD(RTj,RSPl,NUdata,varargin)
%  calculate scene radiance at the instrument aperture and divide by flat
%  field
%   Input Parameters
%     RTj: scene image [counts/ms]
%     RSPl : binned spectral radiometric responsitivity
%     NUdata: CRISMdata obj, CDR NU data
%   Output Parameters
%     RDm  : calculated scene radiance
%     FF   : flat fielding component
%   OPTIONAL PARAMETERS
%     'FLAT_FIELD' : 
%      boolean, whether or not to perform flat field correction using 
%      NUdata or not.
%      (default) true
%   
%   *Detail*
%    see Section 2.15 of "CRISM_DPSIS_Appendix_L_v5_2016-09-01.pdf"

flat_field = true;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'FLAT_FIELD'
                flat_field = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

if isempty(NUdata.img), NUdata.readimg(); end

[L,S,Bands] = size(RTj);
[L1,S1,Bands1] = size(RSPl);

if (S~=S1) || (Bands~=Bands1)
    error('Size of RTj and RSPl does not match');
end

if flat_field
    FL = NUdata.img(1,:,:);
    EP = NUdata.img(2,:,:);

    % FF = zeros(size(RTj));
    % for l=1:L
    %     FF(l,:,:) = FL + EP .* log(RTj(l,:,:) ./ RSPl) ;
    % end
    
    if all(abs(EP)<1e-10)
        FF = FL;
    else
        FF = FL + EP .* log(RTj./RSPl) ;
    end
        
    % FF = repmat(FL,[L,1,1]) + repmat(EP,[L,1,1]) .* log( RTj ./ repmat(RSPl,[L,1,1]) );

    % RDm = RTj ./ (FF .* repmat(RSPl,[L,1,1]));
    RDm = RTj ./ (FF.*RSPl);
else
    RDm = RTj ./ RSPl;
    FF = [];
end

end
