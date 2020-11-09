
% Create a protocol
ProtocolName='infant_templates';
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

study   = bst_get('Study');
ProtocolSubjects = bst_get('ProtocolSubjects');   
ProtocolInfo = bst_get('ProtocolInfo');

fid_dict = load('/home/christian/brainstorm_db/fiducials.mat');
files = dir("/home/christian/synchedin/infants_atlas_modeling/infant_template_paper/Templates/ANTS*3T");

montages = {'10-5', '10-10', '10-20', 'HGSN128'};

delete '/home/christian/.brainstorm/defaults/anatomy/*'

for noSubject=1:size(files, 1)
    
    subjectName = files(noSubject).name;
    
    nbSubjects = length(ProtocolSubjects.Subject);
    iSubject = noSubject; %nbSubjects + 1;
    sSubject = db_template('Subject');
    sSubject.Name     = subjectName;
    sSubject.FileName = bst_fullfile(subjectName, 'brainstormsubject.mat');
    % Get subjects defaults from protocol
    sSubject.UseDefaultAnat    = ProtocolInfo.UseDefaultAnat;
    sSubject.UseDefaultChannel = ProtocolInfo.UseDefaultChannel;    
    
    sSubject = db_add_subject(sSubject, iSubject);
   
    FsDir = [files(noSubject).folder '/' files(noSubject).name];    
    [BstMriFile, sMri] = import_mri(noSubject, [FsDir '/mri/T1.mgz']);
    
    fiducials_world = getfield(fid_dict, strrep(strrep(subjectName,'ANTS','s'), '-', '_'));
    fiducials_mri = cs_convert(sMri, 'world', 'mri', fiducials_world ./ 1000) .* 1000;
    
	sFid = struct('NAS', fiducials_mri(1, :), ...
                  'LPA', fiducials_mri(2, :), ...
                  'RPA', fiducials_mri(3, :), ...
                  'AC',  fiducials_mri(4, :), ...
                  'PC',  fiducials_mri(5, :), ...
                  'IH',  fiducials_mri(6, :));  

    import_anatomy_fs(iSubject, FsDir, 15000, 0, sFid);

    iStudy = db_add_condition(subjectName, 'test');
    sStudy = bst_get('Study', iStudy);    
    
    % Create folder for channel file
    for no_montage = 1:length(montages)
        montage = [montages{no_montage} '-montage.fif'];
        RefChannelMat = db_template('channelmat');

        % Read channel file
        FifMat = import_channel([], [FsDir '/' montage], 'FIF', 0, 0, 0);

        % Create destination channel file
        ChannelMat = RefChannelMat;
        ChannelMat.SCS.NAS = FifMat.HeadPoints.Loc(:,2);
        ChannelMat.SCS.LPA = FifMat.HeadPoints.Loc(:,1);
        ChannelMat.SCS.RPA = FifMat.HeadPoints.Loc(:,3);
        for i = 15:size(FifMat.HeadPoints.Loc, 2)
            ChannelMat.Channel(i-14).Loc = FifMat.HeadPoints.Loc(:,i);
        end
        ChannelMat.HeadPoints = FifMat.HeadPoints;

        % Save channel file
        ChannelFile = bst_fullfile(bst_fileparts(file_fullpath(sStudy.FileName)), montages{no_montage});
        bst_save(ChannelFile, ChannelMat, 'v7');
    end
    
    db_reload_studies(iStudy);    

    bst_memory('UnloadAll', 'Forced');
    db_reload_subjects(iSubject);
    panel_protocols('UpdateTree');
    db_save();    
    
    export_default_anat(iSubject, subjectName, 1);    
end
