function [IoF_woc,RDn_woc,IoF_bk1_o,IoF_bk2_o] = crism_calibration_IR_v2(obs_id,varargin)
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
%    'SAVE_FILE'
%       Boolean, save files or not
%       (default) true;
%    'VERSION'
%       version of the data
%       (default) '6'
%    'PRODUCT_TYPE'
%       product type of the data
%       (default) 'TRR'
%    'MODE'
%       mode of the calibration {'yuki','original'}
%       'yuki': crmcal_pipeline_IR_IF_wTRRIF_yuki will be used.
%       'original': crmcal_pipeline_IR_original will be used.
%       (default) 'yuki'
%    'SAVE_MEMORY'
%       saving memory or not. true or false
%       (default) true
%    'DWLD','DOWNLOAD'
%       if download the data or not, 2: download, 1:access an only show the
%       path, 0: nothing
%       (default) 0
%   'Force_dwld'
%       binary, whether or not to force performing pds_downloader.
%       (default) false
%   'DWLD_OVERWRITE'
%       if overwrite the file if exists
%       (default) 0
%   'DWLD_INDEX_CACHE_UPDATE'
%       boolean, whether or not to update index.html 
%       (default) false
%   'VERBOSE'
%       boolean, whether or not to show operations.
%       (default) true
%   OUTPUTS
%    IoF_woc     : [L x S x B] I/F image cube, single
%    RDn_woc     : [L x S x B] Radiance image cube, single
%    IoF_bk1_o   : [L_D x S x B] I/F image cube of DF1 measurement, single
%    IoF_bk2_o   : [L_D x S x B] I/F image cube of DF2 measurement, single
%                  empty if not exist.

global crism_env_vars
if ~isfield(crism_env_vars,'dir_TRRX')
    dir_trrx = -1;
else
    dir_trrx = crism_env_vars.dir_TRRX;
end


%% setup initial values
save_pdir = dir_trrx;
save_dir_yyyy_doy = true;
save_file = true;
force = false;
skip_ifexist = false;
save_mem = true;
mode_calib = 'yuki';
% apbprmvl = 'HighOrd';
vr = '';
product_type = 'TRR';
ffc_counter = 1;
OBS_COUNTER_SCENE_custom = 0;
OBS_COUNTER_DF_custom = 0;
verbose = true;
dwld = 0;
force_dwld = 0;
dwld_index_cache_update = 0;
dwld_overwrite = 0;
force_dwld = 0;
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
            case 'SAVE_FILE'
                save_file = varargin{i+1};
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
            case 'OBS_COUNTER_SCENE'
                obs_counter_tmp = varargin{i+1};
                OBS_COUNTER_SCENE_custom = 1;
            case 'OBS_COUNTER_DF'
                obs_counter_df_tmp = varargin{i+1};
                OBS_COUNTER_DF_custom = 1;
                
            case 'FORCE_DWLD'
                force_dwld = varargin{i+1};
            case {'DOWNLOAD','DWLD'}
                dwld = varargin{i+1};
            case 'DWLD_INDEX_CACHE_UPDATE'
                dwld_index_cache_update = varargin{i+1};
            case 'DWLD_OVERWRITE'
                dwld_overwrite = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

if save_file
    if save_pdir==-1
        error(['No default directory crism_env_vars.dir_YUK is defined. ' ...
            'Please update crismToolbox.json and rerun crism_init, ' ...
            'or manually specify "SAVE_PDIR".']);
    elseif ~exist(save_pdir,'dir')
        [status] = mkdir(save_pdir);
        if status
            if verbose, fprintf('"%s" is created.\n',save_pdir); end
            chmod777(save_pdir,verbose);
        else
            error('Failed to create %s',save_pdir);
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
        case 'yuki4'
            vr = 'D';
        case 'original'
            vr = 'O';
        otherwise
            error('mode_calib %s is not defined.',mode_calib);
    end
end

%%

[ yyyy_doy,obs_classType ] = crism_searchOBSID2YYYY_DOY_v2(obs_id);

switch obs_classType
    case {'FRT','HRL','HRS'}
        obs_counter = '07';
        obs_counter_epf = '[0-689A-Za-z]{2}';
        obs_counter_epfdf = '0[0E]{1}';
        obs_counter_df = '0[68]{1}';
    case {'FRS','ATO'}
        obs_counter = '01';
        obs_counter_df = '0[03]{1}';
        obs_counter_epf = '';
        obs_counter_epfdf = '';
        obs_counter_un = '02';
    case 'FFC'
        obs_counter = '0[13]{1}';
        obs_counter_df = '0[024]{1}';
        % this could be switched.
        obs_counter_epf = '';
        
    case 'CAL'
        obs_counter = '[0-9a-fA-F]{2}';
        obs_counter_df = '[0-9a-fA-F]{1}';
    case 'ICL'
        obs_counter = '[0-9a-fA-F]{2}';
        obs_counter_df = '[0-9a-fA-F]{1}';
    case {'MSP','HSP'}
        obs_counter = '01';
        obs_counter_df = '0[02]{1}';
        obs_counter_epf = '';
        obs_counter_epfdf = '';
    otherwise
        error('OBS_TYPE %s is not supported yet.',obs_classType);
end

if OBS_COUNTER_SCENE_custom
    obs_counter = obs_counter_tmp;
end
if OBS_COUNTER_DF_custom
    obs_counter_df = obs_counter_df_tmp;
end

%%
crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L',...
    'obs_counter_scene',obs_counter,'obs_counter_df',obs_counter_df);
% crism_obsS = CRISMObservation(obs_id,'SENSOR_ID','S');
switch upper(crism_obs.info.obs_classType)
    case {'FRT','HRL','HRS','FRS','ATO','MSP','HSP'}
        TRRIFdata = get_CRISMdata(crism_obs.info.basenameIF ,'');
        TRRRAdata = get_CRISMdata(crism_obs.info.basenameRA ,'');
        DDRdata   = get_CRISMdata(crism_obs.info.basenameDDR,'','VERSION',1);
    case {'FFC'}
        [TRRIFdata] = get_scene_CRISMdata_FFC(crism_obs.info.basenameIF ,'',ffc_counter);
        [TRRRAdata] = get_scene_CRISMdata_FFC(crism_obs.info.basenameRA ,'',ffc_counter);
        [DDRdata]   = get_scene_CRISMdata_FFC(crism_obs.info.basenameDDR,'',ffc_counter);
end
TRRIF_is_empty = isempty(TRRIFdata);
if TRRIF_is_empty
    TRRIFdata = TRRRAdata;
end
[DFdata1,DFdata2] = crism_get_DFdata4SC(TRRIFdata,crism_obs);
%% setting saving directory
if save_file
    if save_dir_yyyy_doy
        dirpath_yyyy_doy = joinPath(save_pdir,crism_obs.info.yyyy_doy);
        if ~exist(dirpath_yyyy_doy,'dir')
            status = mkdir(dirpath_yyyy_doy);
            if status
                if verbose, fprintf('"%s" is created.\n',dirpath_yyyy_doy); end
                chmod777(dirpath_yyyy_doy,verbose);
            else
                error('Failed to create %s',dirpath_yyyy_doy);
            end
        end
        save_dir = joinPath(dirpath_yyyy_doy,crism_obs.info.dirname);
        if ~exist(save_dir,'dir')
            status = mkdir(save_dir);
            if status
                if verbose, fprintf('"%s" is created.\n',save_dir); end
                chmod777(save_dir,verbose);
            else
                error('Failed to create %s',save_dir);
            end
        end
    else
        save_dir = joinPath(save_pdir,crism_obs.info.dirname);
        if ~exist(save_dir,'dir')
            status = mkdir(save_dir);
            if status
                if verbose, fprintf('"%s" is created.\n',save_dir); end
                chmod777(save_dir,verbose);
            else
                error('Failed to create %s',save_dir);
            end
        end
    end
end
propIF = TRRIFdata.prop;
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

propDF1_IF = DFdata1.prop;
propDF1_IF.activity_id = 'IF';
propDF1_IF.product_type = product_type;
propDF1_IF.version = vr;
bnameDF1_IF = crism_get_basenameOBS_fromProp(propDF1_IF);


if ~any(strcmpi(TRRIFdata.lbl.OBSERVATION_TYPE,{'FRS','ATO'}))
    propDF2_IF = DFdata2.prop;
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
else
    fpath_TRRYIFDF2 = './';
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

%% main processing
switch lower(mode_calib)
    case 'yuki'
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_wTRRIF_yuki(...
        TRRIFdata,TRRIFdataS,crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',dwld,'FORCE_DWLD',force_dwld, ...
        'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose, ...
        'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,...
        'MODE_SP','SOC', 'APBPRMVL','HighOrd','SATURATION_RMVL',0,...
        'BK_SATURATION_RMVL',1,'BK_BPRMVL',false,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
        % [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_yuki(...
        %     TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
        %     'APBPRMVL','HighOrd','SAVE_MEMORY',save_mem,'SATURATION_RMVL',0,...
        %     'BK_SATURATION_RMVL',1,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
    case 'yukiz'
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_wTRRIF_yuki(...
        TRRIFdata,TRRIFdataS,crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',dwld,'FORCE_DWLD',force_dwld,...
        'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose, ...
        'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,...
        'MODE_SP','SOC', 'APBPRMVL','none','SATURATION_RMVL',0,...
        'BK_SATURATION_RMVL',1,'BK_BPRMVL',false,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
        % [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_yuki(...
        %    TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
        %     'APBPRMVL','none','SAVE_MEMORY',save_mem,'SATURATION_RMVL',0,...
        %     'BK_SATURATION_RMVL',1,'BK_MEAN_ROBUST',0,'BK_MEAN_DN14',0);
    case 'yuki2'
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_wTRRIF_yuki(...
        TRRIFdata,[],crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',dwld,'FORCE_DWLD',force_dwld,...
        'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose, ...
        'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,...
        'MODE_SP','SOC', 'APBPRMVL','HighOrd','SATURATION_RMVL',2,...
        'BK_SATURATION_RMVL',2,'BK_BPRMVL',false,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',0);
        % [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_yuki(...
        %     TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,2,...
        %     'APBPRMVL','HighOrd','SAVE_MEMORY',save_mem,'SATURATION_RMVL',2,...
        %     'BK_SATURATION_RMVL',2,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',0);
    case 'yuki3'
        %custom SPdata
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_wTRRIF_yuki(...
        TRRIFdata,[],crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',dwld,'FORCE_DWLD',force_dwld, ...
        'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose, ...
        'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,...
        'MODE_SP','MAN', 'APBPRMVL','HighOrd','SATURATION_RMVL',2,...
        'SP_SATURATION_RMVL',2,'SP_APBPRMVL','HighOrd','SP_MEAN_ROBUST',1,'SP_MEAN_DN14',1,...
        'BK_SATURATION_RMVL',2,'BK_BPRMVL',false,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',1);
    case 'yuki4'
        %custom SPdata
        bkoption = 2;
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_wTRRIF_yuki(...
        TRRIFdata,[],crism_obs,bkoption,...
        'SAVE_MEMORY',save_mem,'dwld',dwld,'FORCE_DWLD',force_dwld, ...
        'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose, ...
        'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,...
        'MODE_SP','MAN', 'APBPRMVL','HighOrd','SATURATION_RMVL',2,'FLAT_FIELD',false,...
        'SP_SATURATION_RMVL',2,'SP_APBPRMVL','HighOrd','SP_MEAN_ROBUST',1,'SP_MEAN_DN14',1,...
        'BK_SATURATION_RMVL',2,'BK_BPRMVL',false,'BK_MEAN_ROBUST',1,'BK_MEAN_DN14',1);
    case 'original'
        error('original mode is not supported yet');
        %[RDn,RDn_woc] = crmcal_pipeline_IR_original(TRRIFdata,EDRdata,...
        %    DFdata1,DFdata2,BKdata1,BKdata2,1,'APBPRMVL','HighOrd',...
        %    'SAVE_MEMORY',save_mem);
        %RDn_woc = RDn;
    otherwise
        error('undefined mode %s',mode_calib);
end

% raidance to if
d_km = DDRdata.lbl.SOLAR_DISTANCE.value;
[ d_au ] = km2au( d_km );
SFdata = TRRIFdata.readCDR('SF');
[IoF]     = crmcal_rd2if(RDn,SFdata,d_au);
[IoF_woc] = crmcal_rd2if(RDn_woc,SFdata,d_au);
% save(fpath_TRRYIF_mat,'IoF','IoF_woc');

%% saving files
IoF_woc = single(IoF_woc);
RDn_woc = single(RDn_woc);
if save_file
    [hdrif_cat] = crism_const_cathdr(TRRIFdata,false);
    fprintf('Saving %s ...\n',fpath_TRRYIF_lbl);
    % envihdrwritex(hdrif_cat,joinPath(save_dir,[bnameIF '.hdr']),'OPT_CMOUT','false');
    copyfile(TRRIFdata.lblpath,fpath_TRRYIF_lbl);
    chmod777(fpath_TRRYIF_lbl,verbose);
    fprintf('Done\n');
    fprintf('Saving %s ...\n',fpath_TRRYIF_img);
    envidatawrite(IoF_woc,fpath_TRRYIF_img,hdrif_cat);
    chmod777(fpath_TRRYIF_img,verbose);
    fprintf('Done\n');

    [hdrra_cat] = crism_const_cathdr(TRRRAdata,false);
    fprintf('Saving %s ...\n',fpath_TRRYRA_lbl);
    % envihdrwritex(hdrra_cat,joinPath(save_dir,[bnameRA '.hdr']),'OPT_CMOUT','false');
    copyfile(TRRRAdata.lblpath,fpath_TRRYRA_lbl);
    chmod777(fpath_TRRYRA_lbl,verbose);
    fprintf('Done\n');
    fprintf('Saving %s ...\n',fpath_TRRYRA_img);
    envidatawrite(RDn_woc,fpath_TRRYRA_img,hdrra_cat);
    chmod777(fpath_TRRYRA_img,verbose);
    fprintf('Done\n');
end

%%
% conversion for 
switch lower(mode_calib)
    case {'yuki','yuki2','yuki3','yuki4'}
        [IoF_bk1_o] = crmcal_rd2if(RDn_bk1_o,SFdata,d_au);
        if ~isempty(RDn_bk2_o)
            [IoF_bk2_o] = crmcal_rd2if(RDn_bk2_o,SFdata,d_au);
        else
            IoF_bk2_o = [];
        end
        
        IoF_bk1_o = single(IoF_bk1_o);
        if save_file
            fprintf('Saving %s ...\n',fpath_TRRYIFDF1);
            save(fpath_TRRYIFDF1,'IoF_bk1_o');
            chmod777(fpath_TRRYIFDF1,verbose);
            fprintf('Done\n');
        end
        
        if ~isempty(bnameDF2_IF)
            IoF_bk2_o = single(IoF_bk2_o);
            if save_file
                fprintf('Saving %s ...\n',fpath_TRRYIFDF2);
                save(fpath_TRRYIFDF2,'IoF_bk2_o');
                chmod777(fpath_TRRYIFDF2,verbose);
                fprintf('Done\n');
            end        
        end
    otherwise
end

end