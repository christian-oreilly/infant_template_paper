% Script for generating Brainstorm anatomy templates for O'Reilly infant atlases.
% Reference: https://github.com/christian-oreilly/infant_template_paper
% 
% Author: Francois Tadel, 2021

atlas_dir = 'C:\Work\RawData\Test\Oreilly\2021';
output_dir = 'C:\Work\RawData\Test\Oreilly\bst_templates';

% Atlases to import
AtlasNames = {...
    'ANTS1-0Months3T', 'Oreilly_1m_2021'; ...
    'ANTS2-0Months3T', 'Oreilly_2m_2021'; ...
    'ANTS2-0Weeks3T', 'Oreilly_0.5m_2021'; ...
    'ANTS2-0Years3T', 'Oreilly_24m_2021'; ...
    'ANTS3-0Months3T', 'Oreilly_3m_2021'; ...
    'ANTS4-5Months3T', 'Oreilly_4m_2021'; ...
    'ANTS6-0Months3T', 'Oreilly_6m_2021'; ...
    'ANTS7-5Months3T', 'Oreilly_7.5m_2021'; ...
    'ANTS9-0Months3T', 'Oreilly_9m_2021'; ...
    'ANTS10-5Months3T', 'Oreilly_10.5m_2021'; ...
    'ANTS12-0Months3T', 'Oreilly_12m_2021'; ...
    'ANTS15-0Months3T', 'Oreilly_15m_2021'; ...
    'ANTS18-0Months3T', 'Oreilly_18m_2021'};


% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end
% Create output folder if doesn't exist
if ~file_exist(output_dir)
    mkdir(output_dir);
end

% Import atlas in database
for iAtlas = 1:size(AtlasNames,1)
    % Delete existing subject
    SubjectName = AtlasNames{iAtlas,1};
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if ~isempty(sSubject)
        db_delete_subjects(iSubject);
    end
    
    % FreeSurfer folder
    FsDir = fullfile(atlas_dir, SubjectName);
    
    % Read T1 MRI, to get the world transformation
    T1File = fullfile(FsDir, 'mri', 'T1.mgz');
    sMri = in_mri(T1File, 'MGH', 0);
    
    % Change folder to access private function
    curDir = pwd;
    cd(fullfile(fileparts(which('import_channel')), 'private'));
    
    % Read fiducials
    FidFile = fullfile(FsDir, 'bem', [SubjectName '-fiducials.fif']);
    [ fid, tree ] = fiff_open(FidFile);
    dig = fif_read_headpoints(fid, tree);
    fclose(fid);
    
    % Restore folder
    cd(curDir);
    
    % Convert fiducials to MRI coordinates
    fidWorld = double([dig.r]');
    fidMri = cs_convert(sMri, 'world', 'mri', fidWorld);
    iNas = find([dig.ident] == 2);  % FIFF.FIFFV_POINT_NASION
    iLpa = find([dig.ident] == 1);  % FIFF.FIFFV_POINT_LPA
    iRpa = find([dig.ident] == 3);  % FIFF.FIFFV_POINT_RPA

    % Process: Import FreeSurfer folder
    bst_process('CallProcess', 'process_import_anatomy', [], [], ...
        'subjectname', SubjectName, ...
        'mrifile',     {FsDir, 'FreeSurfer+Thick'}, ...
        'nvertices',   15000, ...
        'nas', fidMri(iNas,:) * 1000, ...
        'lpa', fidMri(iLpa,:) * 1000, ...
        'rpa', fidMri(iRpa,:) * 1000);

    % Get subject
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    % Delete unwanted files
    AnatDir = bst_fileparts(file_fullpath(sSubject.FileName));
    MidFiles = file_find(AnatDir, 'tess_cortex_mid_*.mat', 1, 0);
    file_delete(MidFiles, 1, 3);
    db_reload_subjects(iSubject);
    
    % Import BEM surfaces
    BemFiles = {fullfile(FsDir, 'bem', 'inner_skull.surf'), ...
                fullfile(FsDir, 'bem', 'outer_skull.surf'), ...
                fullfile(FsDir, 'bem', 'outer_skin.surf')};
    import_surfaces(iSubject, BemFiles, 'FS', 1);
    
    % List all electrodes files
    TsvFiles = file_find(FsDir, '*_electrodes.tsv', 2, 0);
    % Import each channel file
    for iTsv = 1:length(TsvFiles)
        % Split filename
        [fPath, MontageName, fExt] = bst_fileparts(TsvFiles{iTsv});
        MontageName = strrep(MontageName, '_electrodes', '');
        % Create folder
        iStudy = db_add_condition(SubjectName, MontageName);
        % Import channel file
        ChannelFile = import_channel(iStudy, TsvFiles{iTsv}, 'FREESURFER-TSV', 0, 0, 1, 0, 0);
        % Update file
        if ~isempty(ChannelFile)
            % Get full path
            ChannelFile = file_fullpath(ChannelFile);
            % Update file comment
            ChannelMat = load(ChannelFile);
            ChannelMat.Comment = [MontageName, sprintf(' (%d)', length(ChannelMat.Channel))];
            bst_save(ChannelFile, ChannelMat, 'v7');
            % Rename file based on input file name
            [fPath, fBase, fExt] = bst_fileparts(ChannelFile);
            movefile(ChannelFile, fullfile(fPath, ['channel_', file_standardize(MontageName), '.mat']));
            % Reload folder
            db_reload_studies(iStudy);
        end
    end
    
    % Rename subject
    db_rename_subject(SubjectName, AtlasNames{iAtlas,2});
    % Export as anatomy template
    ZipFile = export_default_anat(iSubject, AtlasNames{iAtlas,2}, 1);
    % Move zip file to output folder
    movefile(ZipFile, output_dir);
end

