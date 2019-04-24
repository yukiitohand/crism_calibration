function [ DN14g ] = nonlinearity_correction( DN14e,LCdata,hkt,varargin )
% [ DN14g ] = nonlinearity_correction( DN14b,LCdata )
%  The forth step of the calibration: nonliearity correcion
%  Input Parameters
%    DN14e  : 14bit DN image (L,S,B) quadrant ghost corrected
%    LCdata : CRISMdata obj, CDR LC data
%    hkt    : housekeeping table data (from TRR3)
%  Output parameters
%    DN14g    : processed DN14e data
%  Optional Parameters
%    'BINNING' : ID of binning mode {0,1,2,3}
%                (default) 0
%    'BINX'    : binning size (PIXEL_AVERAGING_WIDTH in LBL)
%                (default) 1 
%     *Note: The relationship between BINNING and BINX
%        BINNING   BINX
%              0      1
%              1      2
%              2      5
%              3     10
%    If you only specify 'BINNING', 'BINX' is obtained automatically. If
%    you specify 'BINX', whatever 'BINNING' parameter will not be
%    considered.
%  Detail:
%   read section 2.6 of "CRISM_DPSIS_Appendix_L_v5_2016-09-01.pdf" and
%   description of "CDR6_1_0000000000_LC_L_1.LBL"
%   Currently only the simplest case:
%   Eff = epsilon * log(DN) + phi
%    (Eff = a * log(DN) + b in the LBL file)
%   is implemented.
%   
if isempty(LCdata.tab), LCdata.readTAB(); end

binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = get_binning(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end

[L,S,Bands] = size(DN14e);

rate = [hkt.data.RATE];
rate = rate(:);
if L==1
    rate = unique(rate);
    if length(rate)>1
        error('Something wrong');
    end
end

a = rateQuadrantTABformatter(rate,LCdata.tab,'A','BINX',binx);
b = rateQuadrantTABformatter(rate,LCdata.tab,'B','BINX',binx);
c = rateQuadrantTABformatter(rate,LCdata.tab,'C','BINX',binx);
d = rateQuadrantTABformatter(rate,LCdata.tab,'D','BINX',binx);
e = rateQuadrantTABformatter(rate,LCdata.tab,'E','BINX',binx);

DN14e_1 = DN14e;
DN14e_1(DN14e<0) = nan;

% currently only
if all(e(:)==0)
    % small value correction is implemented.
    DN14e_1_lt100  = DN14e_1 < 100;
    DN14e_1(DN14e_1_lt100) = 5*log(exp(4)+exp(DN14e_1(DN14e_1_lt100)/5));
    F = repmat(a,[1,1,Bands]) .* log(DN14e_1) + repmat(b,[1,1,Bands]);
else
    error('The complex case is not implemented yet, only F=a*log(DN)+b is implemented');
end

DN14g = DN14e ./ F;

end
