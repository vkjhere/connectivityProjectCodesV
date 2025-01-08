% This GUI helps compare PPC vs freq between 2 protocols.
function runComparePPCbetweenSt_Bl()
    % Create the GUI
    fig = figure('Name', 'PPC comparision between stimulus and baseline', 'NumberTitle', 'off', ...
                 'Position', [50, 100, 800, 600]);

    % Panel for Controlpanel
    controlPanel = uipanel('Parent', fig, 'Position', [0, 0.92, 1, 0.15], ...
                           'Title', 'Controls');

    % Dropdown for Subject Group
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [10, 20, 120, 20], 'String', 'Subject Group:');
    groupDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [130, 20, 120, 20], ...
                               'String', {'Meditators', 'Controls'});

    % Dropdown for Protocol
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [260, 20, 120, 20], 'String', 'Protocol :');
    protocolDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [380, 20, 120, 20], ...
                                  'String', { 'G1', 'G2', 'M2' });

    % Dropdown for Reference Electrode Group
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [550, 20, 190, 20], 'String', 'Reference Electrode Group:');
    referenceDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [750, 20, 120, 20], ...
                                   'String', {'Left Occipital','Right Occipital','Back Occipital','Left Frontal', 'Right Frontal'});

    % Plot Button
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [950, 20, 80, 20], 'String', 'Plot', ...
              'Callback', @(~, ~) plotConnectivity());
      % Clear Button
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [1100, 20, 80, 20], 'String', 'Clear', ...
              'Callback', @(~, ~) cla_Callback());


    % Axes for Plots
    axesArray = getPlotHandles(1,8,[0.05 0.1 0.9 0.65], 0.02, 0.01, 0);
    for n = 1:8
        axesArray(n) = subplot(1, 8, n, 'Parent', fig);
        pos = get(axesArray(n), 'Position');
        pos(2) = pos(2) - 0.05; % Adjust y-position downwards
        set(axesArray(n), 'Position', pos);
        title(axesArray(n), sprintf('Electrode Group %d', n));
    
    end

    % Plot Connectivity Function
    function plotConnectivity()
        % Get selected group, protocol, and analysis choice
        subjectsChoice = groupDropdown.Value;
        protocolName = protocolDropdown.Value;
  
  
        refElectrodes = referenceDropdown.Value;
        % Map selections to actual strings
        subjectGroups = {'med', 'con'};
        protocols = {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2' };
        
        refElectrodeGroups = {[14 44 47], [19 49 52], [16 17 18 48], [34 36 3],[62 63 30]};
        subjectsChoice = subjectGroups{subjectsChoice};
        protocolName = protocols{protocolName};
        
        
        refElectrodes = refElectrodeGroups{refElectrodes};

      
        %refElectrodes = [16 17 18 48]; % back occipital
        groupType = 'rel';
        cutoffList = [2 25];
        pairedSubjectNameList = getPairedSubjectsBK1;            
        subjectNameLists{1} = pairedSubjectNameList(:,1);
        subjectNameLists{2} = pairedSubjectNameList(:,2);
        pairedDataFlag      = 1;
       
analysisChoice1 = 'bl';
analysisChoice2 = 'st';

if ~exist('refElectrodes','var');         refElectrodes = [];           end
if ~exist('groupType','var');             groupType='rel';              end
if ~exist('connMethod','var');            connMethod = 'ppc';           end
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end

if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [20 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end


if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100]; % freqLims
    axisRangeList{2} = [0 1]; % yLims for 
    axisRangeList{3} = [0 1]; % cLimsTopo
end
if ~exist('cutoffList','var')
    cutoffList = [2 30]; % elcs & trials
end

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 0;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

saveFolderName = 'savedData';
capType = 'actiCap64_UOL';

 %function [dataToPlot,freqVals,connDataElectrodeGroup2,electrodeGroupList,groupNameList, binnedCenters] = fetchConnectivityData(groupChoice, protocolName1, analysisChoice)
[~, ~,connDataElectrodeGroupSt,~,~, ~] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice1,cutoffList,pairedDataFlag,saveFolderName,capType);
[~,freqVals,connDataElectrodeGroupBl,electrodeGroupList,groupNameList, ~] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice2,cutoffList,pairedDataFlag,saveFolderName,capType);

%%% remove line noise %%%%
lineNoiseRange = [46 54];
lineNoiseFreqVals = 46:54;
freqValsNoLineNoise = setdiff(freqVals, lineNoiseFreqVals); % No line noise now.
freqValsPosNoLineNoise = find(ismember(freqVals, freqValsNoLineNoise));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];

titleStr{1} = 'Stimulus ';
titleStr{2} = 'Baseline ';
%titleStr{3} = 'Med - Con';
freqLims = axisRangeList{1};
yLims = axisRangeList{2};
cLimsTopo = axisRangeList{3};
numElectrodeGroups = length(electrodeGroupList);
hConn = axesArray;
%hConn = getPlotHandles(1,numElectrodeGroups,[0.05 0.1 0.9 0.75], 0.02, 0.01, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for m=1:numElectrodeGroups
    if strcmp(subjectsChoice, 'med')
    dataToPlot{1} = squeeze(connDataElectrodeGroupSt{1}(:,m,:)); % med, protocol 1
    dataToPlot{1} = dataToPlot{1}(:,freqValsPosNoLineNoise); % to avoid line noise points
    dataToPlot{2} = squeeze(connDataElectrodeGroupBl{1}(:,m,:)); % med , prot 2
    dataToPlot{2} = dataToPlot{2}(:,freqValsPosNoLineNoise);
    else
    dataToPlot{1} = squeeze(connDataElectrodeGroupSt{2}(:,m,:)); % con, protocol 1
    dataToPlot{1} = dataToPlot{1}(:,freqValsPosNoLineNoise);
    dataToPlot{2} = squeeze(connDataElectrodeGroupBl{2}(:,m,:)); % con , prot 2
    dataToPlot{2} = dataToPlot{2}(:,freqValsPosNoLineNoise);
    end
    
    if displayDataFlag
        % Connectivity as a function of frequency__________________
   displaySettings.colorNames(1,:) = [0 1 0]; % green for protocol 1
   displaySettings.colorNames(2,:) =  [1 0 0]; % red for protocol 2
        displayAndcompareData(hConn(m), dataToPlot,freqVals,displaySettings,yLims,1,useMedianFlag,~pairedDataFlag);
        sgtitle(["PPC of " + subjectsChoice + " during " + protocolName + " -- baseline and stimulus wrt elecs " + num2str(refElectrodes)]);

        title(groupNameList{m});
        xlim(hConn(m),freqLims);
        for j=1:2
            text(10,yLims(2)-0.1*j,[titleStr{j} '(' num2str(size(dataToPlot{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hConn(m));
            
        end
    end
end
end 

function [connData,freqVals,connDataElectrodeGroup,electrodeGroupList,groupNameList,binnedCenters] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice,cutoffList,pairedDataFlag,saveFolderName,capType)

[electrodeGroupList,groupNameList,binnedCenters] = getElectrodeGroupsConn(groupType,refElectrodes,capType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
badSubjectList = cell(1,2);
connDataTMP = cell(1,2);
connDataElectrodeGroupTMP = cell(1,2);
numRefElectrodes = length(refElectrodes);
numElectrodeGroups = length(electrodeGroupList);

for i=1:2 % med and con
    for j=1:length(subjectNameLists{i}) % no. of subjects in 'i'
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']));
        numGoodTrials = tmpData.numGoodTrials;

        if numGoodTrials<cutoffList(2)
            badSubjectList{i}(j) = 1;
        else
            if strcmp(analysisChoice,'bl')
                connDataTMP2 = tmpData.connPre(refElectrodes,:,:);
                freqVals = tmpData.freqPre;
            elseif strcmp(analysisChoice,'st')
                connDataTMP2 = tmpData.connPost(refElectrodes,:,:);
                freqVals = tmpData.freqPost;
            else
                connDataTMP2 = (tmpData.connPre(refElectrodes,:,:) + tmpData.connPost(refElectrodes,:,:))/2;
                freqVals = tmpData.freqPost;
            end

            numGoodElectodes = trace(~isnan(squeeze(connDataTMP2(:,refElectrodes,1))));

            if numGoodElectodes >= cutoffList(1)
                badSubjectList{i}(j) = 0;
                connDataTMP{i}{j} = squeeze(mean(connDataTMP2,1,'omitnan'));
                
                connDataElectrodeGroupTMP2 = zeros(numRefElectrodes,numElectrodeGroups,length(freqVals));
                for e=1:numRefElectrodes
                    for b=1:numElectrodeGroups
                        if strcmp(groupType,'rel')
                            connDataElectrodeGroupTMP2(e,b,:) = squeeze(mean(connDataTMP2(e,electrodeGroupList{e,b},:),2,'omitnan'));
                        else
                            connDataElectrodeGroupTMP2(e,b,:) = squeeze(mean(connDataTMP2(e,electrodeGroupList{b},:),2,'omitnan'));
                        end
                    end
                end
                connDataElectrodeGroupTMP{i}{j} = squeeze(mean(connDataElectrodeGroupTMP2,1,'omitnan'));
            else
                badSubjectList{i}(j) = 1;
            end
        end
    end
end

% Remove bad subjects
connData = cell(1,3);
connDataElectrodeGroup = cell(1,3);

for i=1:2
    if pairedDataFlag
        badSubjectPos = find(sum(cell2mat(badSubjectList')));
    else
        badSubjectPos = find(badSubjectList{i});
    end
    x1 = connDataTMP{i};
    x1(badSubjectPos)=[];
    x2 = connDataElectrodeGroupTMP{i};
    x2(badSubjectPos)=[];

    numSubjects = length(x1);
    for j=1:numSubjects
        connData{i}(j,:,:) = x1{j};
        connDataElectrodeGroup{i}(j,:,:) = x2{j};
    end
end
connData{3} = connData{1} - connData{2}; % diff: med-con
connDataElectrodeGroup{3} = connDataElectrodeGroup{1} - connDataElectrodeGroup{2};
end
function cla_Callback(~,~)
        claGivenPlotHandle(axesArray);
        

        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    cla(plotHandles(ii,j));
                end
            end
        end
    end
end