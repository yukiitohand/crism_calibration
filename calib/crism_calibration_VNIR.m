function [] = crism_calibration_VNIR(obs_id,varargin)
% [] = crism_calibration(crism_obs,varargin)
%   Perform calibration, currently only sensorID='L' is supported. By
%   default, the main scene will be processed.
%
%   Input arguments
%    obs_id: observation ID
%   Optional Parameters
%    'TRRIFdata','TRRRAdata','EDRdata','DFdata1','DFdata2','BKdata1','BKdata2'
%       If you perform calibration on any EPF image, you need to specify
%       those optional parameters manually. If you process the main image
%       (in case of FRT, obs_counter '07', you do not need to do anything)
%    'SAVE_PDIR'
%       parent directory to be stored
%       (default) [localCRISM_PDSrootDir]/../YUK/'
%    'SAVE_DIR_YYYY_DOY'
%       Boolean, if the YYYY_DOY is inserted in the saving directory.
%       (default) true;
%    'VERSION'
%       version of the data
%       (default) '6'
%    'PRODUCT_TYPE'
%       product type of the data
%       (default) 'TRR'
%    'MODE'
%       mode of the calibration {'yuki','original'}
%       'yuki': pipeline_calibration_IR_yuki will be used.
%       'original': pipeline_calibration_IR_original will be used.
%       (default) 'yuki'
%    'SAVE_MEMORY'
%       saving memory or not. true or false
%       (default) true
%     'APBPRMVL'
%       option for a priori bad pixel removal {'HighOrd', 'None'}
%       'HighOrd': the image where a priori bad pixel removal is performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       'None'   : the image where a priori bad pixel removal is NOT performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       (default) 'HighOrd'
%   Output
%    None

global crism_env_vars
dir_yuk = crism_env_vars.dir_YUK;


%% setup initial values
save_pdir = dir_yuk;
save_dir_yyyy_doy = true;
force = false;
skip_ifexist = false;
save_mem = true;
mode_calib = 'yuki';
% apbprmvl = 'HighOrd';

crism_obs = CRISMObservation(obs_id,'SENSOR_ID','S');
if ~isempty(crism_obs.info.basenameIF)
    crism_obs.load_data(crism_obs.info.basenameIF,crism_obs.info.dir_trdr,'if');
    TRRIFdata = crism_obs.data.if;
    TRRIF_is_empty = false;
else
    TRRIFdata = '';
    TRRIF_is_empty = true;
end
if ~isempty(crism_obs.info.basenameRA)
    crism_obs.load_data(crism_obs.info.basenameRA,crism_obs.info.dir_trdr,'ra');
    TRRRAdata = crism_obs.data.ra;
else
    TRRRAdata = '';
end

if ~isempty(crism_obs.info.basenameSC)
    crism_obs.load_data(crism_obs.info.basenameSC,crism_obs.info.dir_edr,'sc');
    EDRdata = crism_obs.data.sc;
else
    EDRdata = '';
end

if ~isempty(crism_obs.info.basenameDDR)
    crism_obs.load_ddr(crism_obs.info.basenameDDR,crism_obs.info.dir_ddr,'ddr');
    DDRdata = crism_obs.data.ddr;
else
    DDRdata = '';
end

if TRRIF_is_empty
    TRRIFdata = TRRRAdata;
end

obs_id_scene = TRRIFdata.get_obsid();
TRRIFdata.load_basenamesCDR();
TRRIFdata.readCDR('BI');
switch EDRdata.lbl.OBSERVATION_TYPE
    case {'FRT','HRL','HRS'}
        crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
        DFdata1 = crism_obs.data.df1;
        crism_obs.load_data(crism_obs.info.basenameDF{2},crism_obs.info.dir_edr,'df2');
        DFdata2 = crism_obs.data.df2;      
        for j=1:length(TRRIFdata.cdr.BI)
            bidata = TRRIFdata.cdr.BI(j);
            if strcmpi(bidata.get_obsid, obs_id_scene)
                if strcmpi(bidata.get_obs_number,'06')
                    BIdata1 = bidata;
                elseif strcmpi(bidata.get_obs_number,'08')
                    BIdata2 = bidata;
                end
            end
        end
    case {'FRS','ATO'}
        if ischar(crism_obs.info.basenameDF)
            crism_obs.load_data(crism_obs.info.basenameDF,crism_obs.info.dir_edr,'df1');
        elseif iscell(crism_obs.info.basenameDF)
            crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
        end
        DFdata1 = crism_obs.data.df1;
        DFdata2 = crism_obs.data.df1;
        for j=1:length(TRRIFdata.cdr.BI)
            bidata = TRRIFdata.cdr.BI(j);
            if strcmpi(bidata.get_obsid, obs_id_scene)
                if strcmpi(bidata.get_obs_number,'00')
                    BIdata1 = bidata; BIdata2 = bidata;
                end
            end
        end
    case {'MSP','HSP'}
        crism_obs.load_data(crism_obs.info.basenameDF{1},crism_obs.info.dir_edr,'df1');
        DFdata1 = crism_obs.data.df1;
        crism_obs.load_data(crism_obs.info.basenameDF{2},crism_obs.info.dir_edr,'df2');
        DFdata2 = crism_obs.data.df2;      
        for j=1:length(TRRIFdata.cdr.BI)
            bidata = TRRIFdata.cdr.BI(j);
            if strcmpi(bidata.get_obsid, obs_id_scene)
                if strcmpi(bidata.get_obs_number,'00')
                    BIdata1 = bidata;
                elseif strcmpi(bidata.get_obs_number,'02')
                    BIdata2 = bidata;
                end
            end
        end
    otherwise
        error('Please define for other cases')
end
vr = '';
product_type = 'TRR';

%% variable input arguments
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'TRRIFDATA'
                TRRIFdata = varargin{i+1};
            case 'TRRRADATA'
                TRRRAdata = varargin{i+1};
            case 'EDRDATA'
                EDRdata = varargin{i+1};
            case 'DFDATA1'
                DFdata1 = varargin{i+1};
            case 'DFDATA2'
                DFdata2 = varargin{i+1};
            case 'BIDATA1'
                BIdata1 = varargin{i+1};
            case 'BIDATA2'
                BIdata2 = varargin{i+1};
            case 'SAVE_PDIR'
                save_pdir = varargin{i+1};
            case 'SAVE_DIR_YYYY_DOY'
                save_dir_yyyy_doy = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case 'VERSION'
                vr = varargin{i+1};
            case 'PRODUCT_TYPE'
                product_type = varargin{i+1};
            case 'SKIP_IFEXIST'
                skip_ifexist = varargin{i+1};
            case 'SAVE_MEMORY'
                save_mem = varargin{i+1};
            case 'MODE'
                mode_calib = varargin{i+1};
%             case 'APBPRMVL'
%                 apbprmvl = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if force && skip_ifexist
    error('You are forcing or skipping? Not sure what you want');
end

if isempty(vr) % vr is hard coded if not set.
    switch lower(mode_calib)
        case 'yuki'
            vr = 'Y';
        case 'yuki2'
            vr = 'B';
        case 'yuki3'
            vr = 'C';
        case 'original'
            vr = 'O';
        otherwise
            error('mode_calib %s is not defined.',mode_calib);
    end
end

%% setting saving directory
if save_dir_yyyy_doy
    save_dir = joinPath(save_pdir,crism_obs.info.yyyy_doy,crism_obs.info.dirname);
else
    save_dir = joinPath(save_pdir,crism_obs.info.dirname);
end

propIF = crism_getProp_basenameOBSERVATION(TRRIFdata.basename);
propIF.product_type = product_type;
propIF.version = vr;
if TRRIF_is_empty
    propIF.activity_id = 'IF';
end
bnameIF = crism_get_basenameOBS_fromProp(propIF);

propRA = crism_getProp_basenameOBSERVATION(TRRRAdata.basename);
propRA.product_type = product_type;
propRA.version = vr;
bnameRA = crism_get_basenameOBS_fromProp(propRA);

propDF1_IF = crism_getProp_basenameOBSERVATION(DFdata1.basename);
propDF1_IF.activity_id = 'IF';
propDF1_IF.product_type = product_type;
propDF1_IF.version = vr;
bnameDF1_IF = crism_get_basenameOBS_fromProp(propDF1_IF);


if ~any(strcmpi(EDRdata.lbl.OBSERVATION_TYPE,{'FRS','ATO'}))
    propDF2_IF = crism_getProp_basenameOBSERVATION(DFdata2.basename);
    propDF2_IF.activity_id = 'IF';
    propDF2_IF.product_type = product_type;
    propDF2_IF.version = vr;
    bnameDF2_IF = crism_get_basenameOBS_fromProp(propDF2_IF);

else
   bnameDF2_IF = '';
end

fpath_TRRYIF_img = joinPath(save_dir,[bnameIF,'.IMG']);
fpath_TRRYIF_lbl = joinPath(save_dir,[bnameIF,'.LBL']);
fpath_TRRYRA_img = joinPath(save_dir,[bnameRA,'.IMG']);
fpath_TRRYRA_lbl = joinPath(save_dir,[bnameRA,'.LBL']);
fpath_TRRYIFDF1 = joinPath(save_dir,[bnameDF1_IF,'.mat']);
if ~isempty(bnameDF2_IF)
    fpath_TRRYIFDF2 = joinPath(save_dir,[bnameDF2_IF,'.mat']);
end
% fpath_TRRYIF_mat = joinPath(save_dir,[bnameIF,'.mat']);

% outputs_fpath = {fpath_TRRYIF_img,fpath_TRRYIF_lbl,fpath_TRRYRA_img,fpath_TRRYRA_lbl,fpath_TRRYIF_mat};
outputs_fpath = {fpath_TRRYIF_img,fpath_TRRYIF_lbl,fpath_TRRYRA_img,fpath_TRRYRA_lbl,fpath_TRRYIFDF1,fpath_TRRYIFDF2};

% examine if all the output files exist.
exist_flg = all(cellfun(@(x) exist(x,'file'),outputs_fpath));

if exist_flg
    if skip_ifexist
        return;
    elseif ~force
        flg = 1;
        while flg
            prompt = sprintf('There exists processed images. Do you want to continue to process and overwrite?(y/n)');
            ow = input(prompt,'s');
            if any(strcmpi(ow,{'y','n'}))
                flg=0;
            else
                fprintf('Input %s is not valid.\n',ow);
            end
        end
        if strcmpi(ow,'n')
            fprintf('Process aborted...\n');
            return;
        elseif strcmpi(ow,'y')
            fprintf('processing continues and will overwrite...\n');
        end
    end
end


if ~exist(save_dir,'dir'), mkdir(save_dir); end

%% main processing
switch lower(mode_calib)
%     case 'yuki'
%         [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_yuki(...
%             TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
%             'APBPRMVL','HighOrd','SAVE_MEMORY',save_mem,'SATURATION_RMVL',0,...
%             'BK_SATURATION_RMVL',1,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
%     case 'yukiz'
%             case 'yuki'
%         [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_yuki(...
%             TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
%             'APBPRMVL','none','SAVE_MEMORY',save_mem,'SATURATION_RMVL',0,...
%             'BK_SATURATION_RMVL',1,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
%     case 'yuki2'
%         [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_yuki(...
%             TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
%             'APBPRMVL','HighOrd','SAVE_MEMORY',save_mem,'SATURATION_RMVL',1,...
%             'BK_SATURATION_RMVL',1,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',0);
%     case 'yuki3'
%         %custom SPdata
%         TRRIFdata.load_basenamesCDR();
%         TRRIFdata.readCDR('SP');
%         for i=1:length(TRRIFdata.cdr.SP)
%             spdata = TRRIFdata.cdr.SP(i);
%             spdata_prop = crism_getProp_basenameCDR4(spdata.basename);
%             switch upper(spdata_prop.sensor_id)
%                 case 'L'
%                     SPdata = spdata;
%                 case 'S'
%                     SPdataVNIR = spdata;
%                 otherwise
%                     error('sensor_id %s is wrong',sensor_id);
%             end
%         end
%         [SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o]...
%           = minipipeline_calibration_IR_SP_wCDRSP_yuki(...
%             SPdata,2,'SAVE_MEMORY',save_mem,'APBPRMVL','HighOrd',...
%             'MEAN_DN14',1,'SATURATION_RMVL',2,...
%             'BK_SATURATION_RMVL',2,'BK_BPRMVL',0,...
%             'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',1);
%         [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_yuki(...
%             TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
%             'APBPRMVL','HighOrd','SAVE_MEMORY',save_mem,'SATURATION_RMVL',2,...
%             'BK_SATURATION_RMVL',2,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',1,...
%             'BK_BPRMVL',0,...
%             'SPdata_o',SPdata_o);
    case 'original'
        [RDn,RDn_woc] = pipeline_calibration_VNIR_original(TRRIFdata,EDRdata,...
            DFdata1,DFdata2,BIdata1,BIdata2,1,'APBPRMVL','HighOrd',...
            'SAVE_MEMORY',save_mem);
        %RDn_woc = RDn;
    otherwise
        error('undefined mode %s',mode_calib);
end

% raidance to if
d_km = crism_obs.data.ddr.lbl.SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );
SFdata = TRRIFdata.readCDR('SF');
[IoF]     = crmcal_rd2if(RDn,SFdata,d_au);
[IoF_woc] = crmcal_rd2if(RDn_woc,SFdata,d_au);
% save(fpath_TRRYIF_mat,'IoF','IoF_woc');

%% saving files
[hdrif_cat] = crism_const_cathdr(TRRIFdata,false);
fprintf('Saving %s ...\n',fpath_TRRYIF_lbl);
% envihdrwritex(hdrif_cat,joinPath(save_dir,[bnameIF '.hdr']),'OPT_CMOUT','false');
copyfile(TRRIFdata.lblpath,fpath_TRRYIF_lbl);
fprintf('Done\n');
fprintf('Saving %s ...\n',fpath_TRRYIF_img);
envidatawrite(single(IoF_woc),fpath_TRRYIF_img,hdrif_cat);
fprintf('Done\n');

[hdrra_cat] = crism_const_cathdr(TRRRAdata,false);
fprintf('Saving %s ...\n',fpath_TRRYRA_lbl);
% envihdrwritex(hdrra_cat,joinPath(save_dir,[bnameRA '.hdr']),'OPT_CMOUT','false');
copyfile(TRRRAdata.lblpath,fpath_TRRYRA_lbl);
fprintf('Done\n');
fprintf('Saving %s ...\n',fpath_TRRYRA_img);
envidatawrite(single(RDn_woc),fpath_TRRYRA_img,hdrra_cat);
fprintf('Done\n');

%%
% conversion for 
switch lower(mode_calib)
    case {'yuki','yuki2','yuki3'}
        [IoF_bk1_o] = crmcal_rd2if(RDn_bk1_o,SFdata,d_au);
        [IoF_bk2_o] = crmcal_rd2if(RDn_bk2_o,SFdata,d_au);
        fprintf('Saving %s ...\n',fpath_TRRYIFDF1);
        IoF_bk1_o = single(IoF_bk1_o);
        save(fpath_TRRYIFDF1,'IoF_bk1_o');
        fprintf('Done\n');
        if ~isempty(bnameDF2_IF)
            fprintf('Saving %s ...\n',fpath_TRRYIFDF2);
            IoF_bk2_o = single(IoF_bk2_o);
            save(fpath_TRRYIFDF2,'IoF_bk2_o');
            fprintf('Done\n');
        end
    otherwise
end
