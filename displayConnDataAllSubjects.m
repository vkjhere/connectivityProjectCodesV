% plots Pairwise Phase Consistency (PPC)
% It calls 'psdOnTopOfPPC' to add PSDs on top -scaled down and shifted up version to fit in the plots
% subjectNameLists - 1x2 cell array containing subjectNames for meditators and controls
% protocolName - one of EO1, EC1, G1, M1, G2, EC2, EO2, M2
% analysisChoice - 'st', 'bl' or 'combined'
% refElectrodes - seed electrodes for which connectivity should be measured and averaged
% connMethod - connectivity method
% Option added to return connectivity and topoplot data. Also to simply return these without displaying here.

function [connDataToReturn,topoplotDataToReturn,connData,freqVals,connDataElectrodeGroup,electrodeGroupList,groupNameList,binnedCenters] = displayConnDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,displayDataFlag)

if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('refElectrodes','var');         refElectrodes = [];           end
if ~exist('groupType','var');             groupType='rel';              end
if ~exist('connMethod','var');            connMethod = 'ppc';           end
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('refChoice', 'var');            refChoice = 'none';           end
if ~exist('stRange','var');               stRange = [0.25 1.25];        end
if ~exist('badElectrodeRejectionFlag','var'); badElectrodeRejectionFlag=1;  end
if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [20 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100]; % freqLims
    axisRangeList{2} = [0 1]; % yLims for PPC
    axisRangeList{3} = [0 1]; % cLimsTopo
end
if ~exist('cutoffList','var')
    cutoffList = [2 30]; % elecs & trials
end

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('hAllPlots','var');             hAllPlots = [];               end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 0;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end
numFreqRanges = length(freqRangeList);
freqRangeColors = copper(numFreqRanges);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocol Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
protocolPos = find(strcmp(protocolNameList,protocolName));

if ~strcmp(refChoice,'none')
    protocolPosRef = find(strcmp(protocolNameList,refChoice));
else
    protocolPosRef = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [0.8 0 0.8];      % Purple 
displaySettings.colorNames(2,:) = [0.25 0.41 0.88]; % Cyan
displaySettings.colorNames(3,:) = [0 0 0]; % black
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';
titleStr{3} = 'Med - Con';

freqLims = axisRangeList{1};
%yLims = axisRangeList{2};
yLims = [0 1.1];
cLimsTopo = axisRangeList{3};
capType = 'actiCap64_UOL';

saveFolderName = 'savedData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayDataFlag
    if isempty(hAllPlots)
        electrodeGroupList = getElectrodeGroupsConn(groupType,1);
        numGroups = length(electrodeGroupList);
        hAllPlots.hConn1 = getPlotHandles(1,numGroups,[0.05 0.55 0.6 0.3],0.01,0.01,1);
        hAllPlots.hConn2 = getPlotHandles(numFreqRanges,numGroups,[0.05 0.05 0.6 0.45],0.01,0.01,0);
        hAllPlots.hTopoRef = getPlotHandles(1,2,[0.66 0.55 0.32 0.3],0.002,0.002,1);
        hAllPlots.hTopo  = getPlotHandles(numFreqRanges,3,[0.66 0.05 0.2 0.45],0,0,1);
        hAllPlots.hConn3 = getPlotHandles(numFreqRanges,1,[0.88 0.05 0.1 0.45],0.01,0.01,1);
    else
        hConn1 = hAllPlots.hConn1;
        hConn2 = hAllPlots.hConn2;
        hTopoRef = hAllPlots.hTopoRef;
        hTopo = hAllPlots.hTopo;
        hConn3  = hAllPlots.hConn3;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[connData,freqVals,connDataElectrodeGroup,electrodeGroupList,groupNameList, binnedCenters] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice,cutoffList,pairedDataFlag,saveFolderName,capType);

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
freqValsNoLineNoise = setdiff(freqVals, lineNoiseFreqVals); % removed bad frequency points. No line noise now.
freqValsPosNoLineNoise = find(ismember(freqVals, freqValsNoLineNoise));
%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,~] = psdOnTopOfPPC(subjectNameLists,badEyeCondition,badTrialVersion,stRange, protocolPos, protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffList,pairedDataFlag, hAllPlots, displayDataFlag);

numElectrodeGroups = length(electrodeGroupList);

connDataToReturn = cell(numElectrodeGroups,numFreqRanges);

for i=1:numElectrodeGroups

    dataToPlot{1} = squeeze(connDataElectrodeGroup{1}(:,i,:));
    dataToPlot{2} = squeeze(connDataElectrodeGroup{2}(:,i,:));

    if displayDataFlag
        % Connectivity as a function of frequency__________________
        displayAndcompareData(hConn1(i),dataToPlot,freqVals,displaySettings, yLims, 1,useMedianFlag,~pairedDataFlag);
        title(groupNameList{i});
        xlim(hConn1(i),freqLims);
        if i==1
        ylabel('PPC');
        end
        for j=1:2
            text(10,0.35-0.1*j,[titleStr{j} '(' num2str(size(dataToPlot{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hConn1(i));
            
        end
    end

    % Violin plots__________________________
    for j=1:numFreqRanges
        tmpConn = cell(1,2);
          dataToPlot{1} = squeeze(connDataElectrodeGroup{1}(:,i,:));
          dataToPlot{2} = squeeze(connDataElectrodeGroup{2}(:,i,:));
        for k=1:2
            tmpConn{k} = mean(dataToPlot{k}(:,freqPosList{j}),2);
        end

        connDataToReturn{i,j} = tmpConn;

        if displayDataFlag
            % display violin plots for connectivity
            displaySettings.plotAxes = hConn2(j,i);
            displaySettings.setYLim  = yLims;
            displaySettings.commonYLim = 1;
            displaySettings.yPositionLine = 0.2;
            if ~useMedianFlag
                displaySettings.parametricTest = 1;
                displaySettings.medianFlag = 0;
            else
                displaySettings.parametricTest = 0;
                displaySettings.medianFlag = 1;
            end

            if i==1 && j==3
                displaySettings.showYTicks=1;
                displaySettings.showXTicks=1;
            else
                displaySettings.showYTicks=0;
                displaySettings.showXTicks=0;
            end

             displayViolinPlot(tmpConn,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);

            if i==1
                ylabel(hConn2(j,i),['Conn (' num2str(freqRangeList{j}(1)) '-' num2str(freqRangeList{j}(2)) ' Hz)'],'color',freqRangeColors(j,:));
            end

            % Add lines in ConnVsFreq plots________
            for k=1:2
                line([freqRangeList{j}(k) freqRangeList{j}(k)],yLims,'color',freqRangeColors(j,:),'parent',hConn1(i));
            end
        end
    end
end

 % Display topoplots______________________
topoplotDataToReturn = cell(2,numFreqRanges);

for i=1:numFreqRanges

    dataToCompare = cell(1,2);
    % Topoplots_________
   
    for j=1:3
        if ~pairedDataFlag && j~=3 % not plottting the Med-Con diff (j) for unpaired case
            x = squeeze(mean(connData{j}(:,:,freqPosList{i}),3, 'omitnan'));
        elseif pairedDataFlag
            x = squeeze(mean(connData{j}(:,:,freqPosList{i}),3, 'omitnan'));
        else ; x=[];
        end
            if useMedianFlag
                data = squeeze(median(x,1,'omitnan'));
            else
                data = squeeze(mean(x,1,'omitnan'));
            end
        
        topoplotDataToReturn{i,j} = data;

        if displayDataFlag
                axes(hTopo(i,j)); %#ok<*LAXES>
                topoplot(data,montageChanlocs,'electrodes','on','plotrad',0.6,'headrad',0.6); colorbar;
                %(.., 'maplimits',cLimsTopo,..)
                if i==1
                title(titleStr{j},'color',displaySettings.colorNames(j,:));
                end
        end

        % Conn versus distance__________________
        if j~=3 % not plotting the diff: med-con
        x = squeeze(mean(connDataElectrodeGroup{j}(:,:,freqPosList{i}),3));
        dataToCompare{j} = x;
        end
   end
%     flipdataToCompare = flip(dataToCompare);
    displayAndcompareData(hConn3(i), dataToCompare, 1:numElectrodeGroups, displaySettings,yLims,1,useMedianFlag,~pairedDataFlag);
        if i==1
        title('PPC vs cos\theta');
       end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % if paired dataset, take common bad subjects
    if pairedDataFlag
        badSubjectList{1} = union(find(badSubjectList{1}),find(badSubjectList{2}));
        badSubjectList{2} = badSubjectList{1};
    else 
        for o=1:2
            badSubjectList{o}= find(badSubjectList{o});
        end
    end

% Remove bad subjects ____
connData = cell(1,3);
connDataElectrodeGroup = cell(1,3);

for q=1:2
x1 = connDataTMP{q};
x2 = connDataElectrodeGroupTMP{q};
badSubjectPos = badSubjectList{q};
badSubjectPos = badSubjectPos(badSubjectPos<=length(connDataTMP{q})); % coz sometimes the data of terminal subject might be ubsent in connDataTMP itself. 
x1(badSubjectPos)=[];
x2(badSubjectPos)=[];  

numSubjects = length(x1);

    for p=1:numSubjects
        connData{q}(p,:,:) = x1{p};
        connDataElectrodeGroup{q}(p,:,:) = x2{p};
    end
end

    if pairedDataFlag 
    connData{3} = connData{1} - connData{2}; % diff: med-con
    connDataElectrodeGroup{3} = connDataElectrodeGroup{1} - connDataElectrodeGroup{2};
    end
end