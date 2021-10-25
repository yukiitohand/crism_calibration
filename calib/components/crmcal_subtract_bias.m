function [ DN14a,BI_m ] = crmcal_subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,varargin )
% [ DN14a,BI_m ] = crmcal_subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,varargin )
%  The second step of the calibration: subtract Bias
%  Input parameters:
%    DN14   : 14bit DN image (L,S,B)
%    BIdata : CRISMdata obj, CDR BI data
%    BSdata : CRISMdata obj, CDR BS data
%    DBdata : CRISMdata obj, CDR DB data
%    EBdata : CRISMdata obj, CDR EB data
%    hkt    : housekeeping table data (from TRR3)
%    rownum_table : ROWNUM_TABLE
%  Output parameters
%    DN14a    : processed DN data
%    BI_m   : calculated bias
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
%    Please specify only one of 'BINNING', and 'BINX'
%
%  * Detail *
%  DN14a = DN14 - BI - a*Hv(row_lambda+1-(502/480)*(480-integ_t))
%          - beta_I*(T_aI-T_bI) - beta_J*(T_aJ-T_bJ)
%
%   where
%       DN: 
%           DN value of the image
%       BI: 
%           CDR BI data
%       Hv(): 
%           heaviside step function
%       a0I: = 31
%           coefficient for the step function (CDR BS data)
%           currently fixed to 31 because it is same for everything
%       row_lambda: 
%           detector row number (counting from zero), 
%       integ_t: 
%           MRO:EXPOSURE_PARAMETER (integer, 0-480) (found in LABEL file)
%           EXPOSURE in HKT
%           t[ms]=1000*(502*floor((502/480)*(480-integ))/(502*[frame rate])
%           t: integration time [mili seconds]
%           frame rate: [Hz]
%           maybe 502 (>480) considers detector readout time?
%       beta_I: 
%           coefficient (CDR DB)
%       T_aI:
%           MRO_DETECTOR_TEMPERATURE in lbl ?
%           IR_DETECTOR_TEMP1, IR_DETECTOR_TEMP2 in HKT ?
%           IR detector temperature for the image
%       T_bI:
%           MRO_DETECTOR_TEMPERATURE in lbl
%           IR detector temperature for the CDR BI
%       beta_J:
%           coefficient (CDR EB)
%       T_aJ:
%           MRO_FPE_TEMPERATURE in lbl ?
%           IR_FPU_BOARD_TEMP in HKT ?
%           IR focal plane board temperature for the image (DN)
%       T_bJ:
%           IR focal plane board for the CDR BI
%       
% if isempty(EDRdata.img),EDRdata.readimg(); end;
% if isempty(EDRdata.hkt),EDRdata.readHKT(); end;
if isempty(BIdata.img),BIdata.readimg(); end
if isempty(BSdata.tab),BSdata.readTAB(); end
if isempty(DBdata.tab),DBdata.readTAB(); end
if isempty(EBdata.tab),EBdata.readTAB(); end

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
                error('Unrecognized option: %s', varargin{i});   
        end
    end
end


[L,S,B] = size(DN14);


term2 = repmat(BIdata.img,[L,1,1]);

rate = cat(1,hkt.data.RATE);
a0I = crism_rateQuadrantTABformatter(rate,BSdata.tab,'A0','BINX',binx);
% a0I = 31;
% fprintf('a0I=31 is currently manually set becase it is same for all the cases now');

% row_lambdaList = reshape(1:B,[1,1,B]);
row_lambdaList = reshape(rownum_table,[1 1 B])+1;
% +1 is already performed.
integ_t = cat(1,hkt.data.EXPOSURE);
term_integ_t = (502/480)*(480-integ_t);

term3 = repmat(a0I,[1,1,B]) ...
    .*heaviside(repmat(row_lambdaList,[L,S,1])-repmat(term_integ_t,[1,S,B]));

% what is the right way to do? mean or 
beta_I = crism_rateQuadrantTABformatter(rate,DBdata.tab,'A','BINX',binx);
T_aI = cat(1,hkt.data.IR_DETECTOR_TEMP1);
% T_aI2 = [hkt.data.IR_DETECTOR_TEMP2];
% T_aI = mean([T_aI1;T_aI2],1);
% T_aI = T_aI(:);
T_bI = BIdata.lbl.MRO_DETECTOR_TEMPERATURE;

term4 = repmat(beta_I,[1,1,B]) .* repmat(T_aI-T_bI,[1,S,B]);

beta_J = crism_rateQuadrantTABformatter(rate,EBdata.tab,'A','BINX',binx);
T_aJ = cat(1,hkt.data.IR_FPU_BOARD_TEMP);
T_bJ = BIdata.lbl.MRO_FPE_TEMPERATURE;

term5 = repmat(beta_J,[1,1,B]) .* repmat(T_aJ-T_bJ,[1,S,B]);

BI_m = term2 - term3 - term4 + term5; % + or -?

DN14a = DN14 - BI_m;

end

