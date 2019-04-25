function [] = crism_calibration_IR_v2(obs_id,varargin)
% [] = crism_calibration_IR_v2(obs_id,varargin)
%   Perform calibration for CRISM IR image. %he main scene will be processed.
%   FFC is also supported
%
%   Input arguments
%    obs_id: observation ID
%   Optional Parameters
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
vr = '';
product_type = 'TRR';
ffc_counter = 1;
%% variable input arguments
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'FFC_IF_COUNTER'
                ffc_counter = varargin{i+1};
                if ischar(ffc_counter)
                    ffc_counter = str2num(ffc_counter);
                end
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

%%
crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L');
% crism_obsS = CRISMObservation(obs_id,'SENSOR_ID','S');
switch upper(crism_obs.info.obs_classType)
    case {'FRT','HRL','HRS','FRS','ATO','MSP','HSP'}
        if ~isempty(crism_obs.info.basenameIF)
            TRRIF_is_empty = false;
            TRRIFdata = CRISMdata(crism_obs.info.basenameIF,'');
            % TRRIFdataS = CRISMdata(crism_obsS.info.basenameIF,'');
            TRRRAdata = CRISMdata(crism_obs.info.basenameRA,'');
        elseif ~isempty(crism_obs.info.basenameRA)
            TRRIF_is_empty = true;
            TRRIFdata = CRISMdata(crism_obs.info.basenameRA,'');
            % TRRIFdataS = CRISMdata(crism_obsS.info.basenameRA,'');
            TRRRAdata = CRISMdata(crism_obs.info.basenameRA,'');
        else
            error('Check data');
        end
        DDRdata = CRISMdata(crism_obs.info.basenameDDR,'');
    case {'FFC'}
        switch ffc_counter
            case 1
                if ~isempty(crism_obs.info.basenameIF)
                    TRRIF_is_empty = false;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameIF{1},'');
                    % TRRIFdataS = CRISMdata(crism_obsS.info.basenameIF{1},'');
                    TRRRAdata = CRISMdata(crism_obs.info.basenameRA{1},'');
                elseif ~isempty(crism_obs.info.basenameRA)
                    TRRIF_is_empty = true;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameRA{1},'');
                    % TRRIFdataS = CRISMdata(crism_obsS.info.basenameRA{1},'');
                    TRRRAdata = CRISMdata(crism_obs.info.basenameRA{1},'');
                else
                    error('Check data');
                end
                DDRdata = CRISMdata(crism_obs.info.basenameDDR{1},'');
            case 3
                if ~isempty(crism_obs.info.basenameIF)
                    TRRIF_is_empty = false;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameIF{2},'');
                    % TRRIFdataS = CRISMdata(crism_obsS.info.basenameIF{2},'');
                    TRRRAdata = CRISMdata(crism_obs.info.basenameRA{2},'');
                elseif ~isempty(crism_obs.info.basenameRA)
                    TRRIF_is_empty = true;
                    TRRIFdata = CRISMdata(crism_obs.info.basenameRA{2},'');
                    % TRRIFdataS = CRISMdata(crism_obsS.info.basenameRA{2},'');
                    TRRRAdata = CRISMdata(crism_obs.info.basenameRA{2},'');
                else
                    error('Check data');
                end
                DDRdata = CRISMdata(crism_obs.info.basenameDDR{2},'');
            otherwise
                error('Check data');
        end
end
[DFdata1,DFdata2] = get_DFdata4SC(TRRIFdata,crism_obs);
%% setting saving directory
if save_dir_yyyy_doy
    save_dir = joinPath(save_pdir,crism_obs.info.yyyy_doy,crism_obs.info.dirname);
else
    save_dir = joinPath(save_pdir,crism_obs.info.dirname);
end

propIF = getProp_basenameOBSERVATION(TRRIFdata.basename);
propIF.product_type = product_type;
propIF.version = vr;
if TRRIF_is_empty
    propIF.activity_id = 'IF';
end
bnameIF = get_basenameOBS_fromProp(propIF);

propRA = getProp_basenameOBSERVATION(TRRRAdata.basename);
propRA.product_type = product_type;
propRA.version = vr;
bnameRA = get_basenameOBS_fromProp(propRA);

propDF1_IF = DFdata1.prop;
propDF1_IF.activity_id = 'IF';
propDF1_IF.product_type = product_type;
propDF1_IF.version = vr;
bnameDF1_IF = get_basenameOBS_fromProp(propDF1_IF);


if ~any(strcmpi(TRRIFdata.lbl.OBSERVATION_TYPE,{'FRS','ATO'}))
    propDF2_IF = DFdata2.prop;
    propDF2_IF.activity_id = 'IF';
    propDF2_IF.product_type = product_type;
    propDF2_IF.version = vr;
    bnameDF2_IF = get_basenameOBS_fromProp(propDF2_IF);

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
    case 'yuki'
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_IF_wTRRIF_yuki(...
        TRRIFdata,TRRIFdataS,crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',2,'MODE_SP','SOC',...
        'APBPRMVL','HighOrd','SATURATION_RMVL',0,...
        'BK_SATURATION_RMVL',1,'BK_BPRMVL',false,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
        % [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_yuki(...
        %     TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
        %     'APBPRMVL','HighOrd','SAVE_MEMORY',save_mem,'SATURATION_RMVL',0,...
        %     'BK_SATURATION_RMVL',1,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
    case 'yukiz'
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_IF_wTRRIF_yuki(...
        TRRIFdata,TRRIFdataS,crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',2,'MODE_SP','SOC',...
        'APBPRMVL','none','SATURATION_RMVL',0,...
        'BK_SATURATION_RMVL',1,'BK_BPRMVL',false,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
        % [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_yuki(...
        %    TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
        %     'APBPRMVL','none','SAVE_MEMORY',save_mem,'SATURATION_RMVL',0,...
        %     'BK_SATURATION_RMVL',1,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
    case 'yuki2'
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_IF_wTRRIF_yuki(...
        TRRIFdata,[],crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',2,'MODE_SP','SOC',...
        'APBPRMVL','HighOrd','SATURATION_RMVL',2,...
        'BK_SATURATION_RMVL',2,'BK_BPRMVL',false,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',0);
        % [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_yuki(...
        %     TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
        %     'APBPRMVL','HighOrd','SAVE_MEMORY',save_mem,'SATURATION_RMVL',2,...
        %     'BK_SATURATION_RMVL',2,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',0);
    case 'yuki3'
        %custom SPdata
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_IF_wTRRIF_yuki(...
        TRRIFdata,TRRIFdataS,crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',2,'MODE_SP','MAN',...
        'APBPRMVL','HighOrd','SATURATION_RMVL',2,...
        'SP_SATURATION_RMVL',2,'SP_APBPRMVL','HighOrd','SP_MEAN_ROBUST',1,'SP_MEAN_DN14',1,...
        'BK_SATURATION_RMVL',2,'BK_BPRMVL',false,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',1);
    case 'original'
        error('original mode is not supported yet');
        %[RDn,RDn_woc] = pipeline_calibration_IR_original(TRRIFdata,EDRdata,...
        %    DFdata1,DFdata2,BKdata1,BKdata2,1,'APBPRMVL','HighOrd',...
        %    'SAVE_MEMORY',save_mem);
        %RDn_woc = RDn;
    otherwise
        error('undefined mode %s',mode_calib);
end

% raidance to if
d_km = DDRdata.lbl.SOLAR_DISTANCE{1};
[ d_au ] = km2au( d_km );
SFdata = TRRIFdata.readCDR('SF');
[IoF] = rd2if(RDn,SFdata,d_au);
[IoF_woc] = rd2if(RDn_woc,SFdata,d_au);
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
        [IoF_bk1_o] = rd2if(RDn_bk1_o,SFdata,d_au);
        [IoF_bk2_o] = rd2if(RDn_bk2_o,SFdata,d_au);
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

end