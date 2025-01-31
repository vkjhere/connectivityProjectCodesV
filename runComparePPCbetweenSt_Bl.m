% This GUI helps compare PPC vs freq between stimulus and baseline for a
% given protocol
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
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [520, 20, 190, 20], 'String', 'Reference Electrode Group:');
    referenceDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [710, 20, 120, 20], ...
                                   'String', {'Left Occipital','Right Occipital','Back Occipital','Left Frontal', 'Right Frontal'});
    % Dropdown for FreqRange
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [850, 20, 190, 20], 'String', 'Freq Range:');
    uicontrol('Parent', controlPanel, 'Style', 'edit', 'Position', [900, 20, 50, 20], 'String', '50', ...
              'Tag', 'freqMin');
    uicontrol('Parent', controlPanel, 'Style', 'edit', 'Position', [950, 20, 50, 20], 'String', '60', ...
              'Tag', 'freqMax');
    % Plot Button
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [1150, 20, 80, 20], 'String', 'Plot', ...
              'Callback', @(~, ~) plotConnectivity());

%plot handles _______________________________________
        hConnPPC = getPlotHandles(1,8,[0.05 0.4 0.6 0.45], 0.02, 0.01, 1);
        hTopoRef = getPlotHandles(1,2,[0.66 0.55 0.32 0.3],0.002,0.002,1);
        hTopo  = getPlotHandles(1,3,[0.1 0.05 0.5 0.3],0.02,0,1);
        hConnDistance = getPlotHandles(1,1,[0.7 0.1 0.2 0.25],0.01,0.01,1);
      % Clear Button____________________________________
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [1300, 20, 80, 20], 'String', 'Clear', ...
              'Callback', @(~, ~) cla_Callback());

    % Plot Connectivity Function
    function plotConnectivity()
        % Get selected group, protocol, and analysis choice
        subjectsChoice = groupDropdown.Value;
        protocolName = protocolDropdown.Value;
        refElectrodes = referenceDropdown.Value;
        freqMin = str2double(findobj('Tag', 'freqMin').String);
        freqMax = str2double(findobj('Tag', 'freqMax').String);
        freqRange = [freqMin, freqMax];

        % Map selections to actual strings
        subjectGroups = {'med', 'con'};
        protocols = {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2' };
        refElectrodeGroups = {[14 44 47], [19 49 52], [16 17 18 48], [34 36 3],[62 63 30]};
        freqRanges = {[7 10],[15 25],[60 90]};

        subjectsChoice = subjectGroups{subjectsChoice};
        protocolName = protocols{protocolName};
        refElectrodes = refElectrodeGroups{refElectrodes};

      
        %refElectrodes = [16 17 18 48]; % back occipital
        groupType = 'rel';
        cutoffList = [2 25];
        displaySignificanceFlag = 1;
        pairedSubjectNameList = getPairedSubjectsBK1;            
        subjectNameLists{1} = pairedSubjectNameList(:,1);
        subjectNameLists{2} = pairedSubjectNameList(:,2);
        pairedDataFlag      = 1;
       
analysisChoice1 = 'st';
analysisChoice2 = 'bl';

if ~exist('refElectrodes','var');         refElectrodes = [];           end
if ~exist('groupType','var');             groupType='rel';              end
if ~exist('connMethod','var');            connMethod = 'ppc';           end
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end

if ~exist('freqRangeList','var')
    freqRangeList{1} = freqRange; % takes from GUI
%     freqRangeList{2} = [20 34]; % SG
%     freqRangeList{3} = [35 65]; % FG
end
numFreqRanges = length(freqRangeList);
freqRangeColors = copper(numFreqRanges);

if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100]; % freqLims
    axisRangeList{2} = [0 1]; % yLims for 
    axisRangeList{3} = [0 1]; % cLimsTopo
end
if ~exist('cutoffList','var')
    cutoffList = [2 30]; % elcs & trials
end

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 1;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

saveFolderName = 'savedData';
capType = 'actiCap64_UOL';


[connDataSt, ~,connDataElectrodeGroupSt,~,~, binnedCenters] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice1,cutoffList,pairedDataFlag,saveFolderName,capType);
[connDataBl,freqVals,connDataElectrodeGroupBl,electrodeGroupList,groupNameList, ~] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice2,cutoffList,pairedDataFlag,saveFolderName,capType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];

titleStr{1} = 'Stimulus ';
titleStr{2} = 'Baseline ';
titleStr{3} = 'difference';
freqLims = axisRangeList{1};
yLims = axisRangeList{2};
cLimsTopo = axisRangeList{3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayDataFlag
    % Display electrode positions
    if strcmp(groupType,'rel')  
        montageChanlocs = showElectrodeGroups(hTopoRef,capType,electrodeGroupList(1,:),groupNameList); % Just show for the first electrode
    else
        montageChanlocs = showElectrodeGroups(hTopoRef,capType,electrodeGroupList,groupNameList);
    end
end
%%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%%%%%
freqPosList = cell(1,numFreqRanges);
lineNoiseRange = [46 54];
lineNoiseFreqVals = 46:54;
badFreqPos = intersect(find(freqVals>=lineNoiseRange(1)),find(freqVals<=lineNoiseRange(2)));
for i = 1:numFreqRanges
    freqPosList{i} = setdiff(intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<=freqRangeList{i}(2))),badFreqPos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numElectrodeGroups = length(electrodeGroupList);

for m=1:numElectrodeGroups
    if strcmp(subjectsChoice, 'med')
    dataToPlot{1} = squeeze(connDataElectrodeGroupSt{1}(:,m,:)); % med, St
    dataToPlot{2} = squeeze(connDataElectrodeGroupBl{1}(:,m,:)); % med , Bl
    else
    dataToPlot{1} = squeeze(connDataElectrodeGroupSt{2}(:,m,:)); % con, st
    dataToPlot{2} = squeeze(connDataElectrodeGroupBl{2}(:,m,:)); % con , bl
    end
    
    if displayDataFlag
        % Connectivity as a function of frequency__________________
   displaySettings.colorNames(1,:) = [0 1 0]; % green for St
   displaySettings.colorNames(2,:) =  [1 0 0]; % red for Bl
   displaySettings.colorNames(3,:) =  [0 0 0]; % black for diff

        displayAndcompareData(hConnPPC(m), dataToPlot,freqVals,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,~pairedDataFlag);
        sgtitle(["PPC of " + subjectsChoice + " during " + protocolName + " -- baseline and stimulus wrt elecs " + num2str(refElectrodes)]);

        title(groupNameList{m});
        xlim(hConnPPC(m),freqLims);
        if m==1
        ylabel('PPC'); xlabel('frequency');
        end
        for j=1:2
            text(10,yLims(2)-0.1*j,[titleStr{j} '(' num2str(size(dataToPlot{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hConnPPC(m));
            
        end
    end
end
 % Display topoplots______________________
topoplotDataToReturn = cell(2,numFreqRanges);
    if strcmp(subjectsChoice, 'med')
            xtemp{1} = connDataSt{1}; % st for meditators
            xtemp{2} = connDataBl{1}; % bl
            xtemp{3} = xtemp{1} -xtemp{2}; % St - Bl
    else 
            xtemp{1} = connDataSt{2}; % st 
            xtemp{2} = connDataBl{2}; % bl
            xtemp{3} = xtemp{1} -xtemp{2}; % St - Bl
    end
for i=1:numFreqRanges

    dataToCompare = cell(1,2);
    % Topoplots_________

    for j=1:3 % St, Bl, and difference
            x = squeeze(mean(xtemp{j}(:,:,freqPosList{i}),3, 'omitnan'));
            if useMedianFlag
                data = squeeze(median(x,1,'omitnan'));
            else
                data = squeeze(mean(x,1,'omitnan'));
            end
        
        topoplotDataToReturn{i,j} = data;

        if displayDataFlag
                axes(hTopo(i,j)); %#ok<*LAXES>
                colormap('jet');
                topoplot(data,montageChanlocs,'electrodes','on','plotrad',0.6,'headrad',0.6); colorbar;
                %(.., 'maplimits',cLimsTopo,..)
                if i==1
                title(titleStr{j},'color',displaySettings.colorNames(j,:));
                end
        end

        % Conn versus distance_______________________
        %x = squeeze(mean(connDataElectrodeGroupSt{j}(:,:,freqPosList{i}),3));
          if strcmp(subjectsChoice, 'med')
            dataToCompare{1} = squeeze(mean(connDataElectrodeGroupSt{1}(:,:,freqPosList{i}),3));
            dataToCompare{2} = squeeze(mean(connDataElectrodeGroupBl{1}(:,:,freqPosList{i}),3));
          else
            dataToCompare{1} = squeeze(mean(connDataElectrodeGroupSt{2}(:,:,freqPosList{i}),3));
            dataToCompare{2} = squeeze(mean(connDataElectrodeGroupBl{2}(:,:,freqPosList{i}),3));
          end
        
   end
%     flipdataToCompare = flip(dataToCompare);
    displayAndcompareData(hConnDistance(i), dataToCompare, binnedCenters, displaySettings,yLims,1,useMedianFlag,~pairedDataFlag);
            ylabel('PPC'); xlabel('cos\theta');
        if i==1
        title('PPC vs cos\theta');
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
        claGivenPlotHandle(hConnPPC);
        claGivenPlotHandle(hTopoRef);
        claGivenPlotHandle(hTopo);
        claGivenPlotHandle(hConnDistance);

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