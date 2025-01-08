% groupType - abs or rel. When set to absolute (abs), the electrode groups are fixed (e.g. occipital, central etc). When set to relative (rel), the groups are based on their distance from the seed electrodes

function runDisplayConnDataAllSubjects(groupType)

if ~exist('groupType','var');       groupType='rel';                    end

fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
backgroundColor = 'w'; panelHeight = 0.125;
colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%% Subject Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel1 = uipanel('Title','Subjects','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.025 1-panelHeight 0.15 panelHeight]);

% Comparison - paired or unpaired
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','Comparison','FontSize',fontSizeSmall);
comparisonList = [{'paired'} {'unpaired'}];
hComparison = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 2/3 0.5 1/3],'Style','popup','String',comparisonList,'FontSize',fontSizeSmall);

% Gender - all, male, female
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','text','String','Gender','FontSize',fontSizeSmall);
genderList = [{'all'} {'male'} {'female'}];
hGender = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 1/3 0.5 1/3],'Style','popup','String',genderList,'FontSize',fontSizeSmall);

% Age - all, young, mid
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0 0.5 1/3],'Style','text','String','Age','FontSize',fontSizeSmall);
ageList = [{'all'} {'young'} {'mid'}];
hAge = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 1/3],'Style','popup','String',ageList,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%% Protocol Details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel2 = uipanel('Title','Protocol','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.175 1-panelHeight 0.15 panelHeight]);

% Protocol
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','ProtocolName','FontSize',fontSizeSmall);
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
hProtocol = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 2/3 0.5 1/3],'Style','popup','String',protocolNameList,'FontSize',fontSizeSmall);

% AnalysisChoice
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','text','String','Analysis','FontSize',fontSizeSmall);
analysisChoiceList1 = [{'spontaneous (bl)'} {'stimulus (st)'} {'combined'}];
analysisChoiceList2 = [{'bl'} {'st'} {'combined'}];
hAnalysisChoice = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 1/3 0.5 1/3],'Style','popup','String',analysisChoiceList1,'FontSize',fontSizeSmall);

% RefChoice
refElectrodeList{1} = [14 44 47]; refElectrodeListName{1} = 'LeftOccipital';
refElectrodeList{2} = [19 49 52]; refElectrodeListName{2} = 'RightOccipital';
refElectrodeList{3} = [16 17 18 48]; refElectrodeListName{3} = 'BackOccipital';

refElectrodeList{4} = [34 36 3]; refElectrodeListName{4} = 'LeftFrontal';
refElectrodeList{5} = [62 63 30]; refElectrodeListName{5} = 'RightFrontal';




uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0 0.5 1/3],'Style','text','String','Ref Elecs','FontSize',fontSizeSmall);
hRefChoice = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 1/3],'Style','popup','String',refElectrodeListName,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Bad Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel3 = uipanel('Title','Bad Electrode Condition','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.325 1-panelHeight 0.15 panelHeight]);

% Bad Eye condition
uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','BadEyeCondition','FontSize',fontSizeSmall);
badEyeConditionList1 = [{'eye position (ep)'} {'none (wo)'}]; badEyeConditionList2 = [{'ep'} {'wo'}];
hBadEye = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 2/3 0.5 1/3],'Style','popup','String',badEyeConditionList1,'FontSize',fontSizeSmall);

% Bad Trial Version
uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','text','String','BadTrialVersion','FontSize',fontSizeSmall);
badTrialVersionList = {'v8'};
hBadTrialVersion = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 1/3 0.5 1/3],'Style','popup','String',badTrialVersionList,'FontSize',fontSizeSmall);

% % Bad Electrode Choice
% uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 0 0.5 1/3],'Style','text','String','BadElecChoice','FontSize',fontSizeSmall);
% badElectrodeChoiceList = [{'Reject badElectrodes of protocolName'} {'Reject common badElectrodes of all protocols'} {'Reject badElectrodes of G1'}];
% hBadElectrodeChoice = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 1/3],'Style','popup','String',badElectrodeChoiceList,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%% Freq Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel4 = uipanel('Title','Freq Ranges','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.475 1-panelHeight 0.15 panelHeight]);
freqRangeList0{1} = [7 10];
freqRangeList0{2} = [16 28];
freqRangeList0{3} = [60 98];

numFreqRanges = length(freqRangeList0);
hFreqRangeMin = cell(1,numFreqRanges);
hFreqRangeMax = cell(1,numFreqRanges);

for i=1:numFreqRanges
    uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 1-i/numFreqRanges 0.5 1/numFreqRanges],'Style','text','String',['Freq Range' num2str(i)],'FontSize',fontSizeSmall);
    hFreqRangeMin{i} = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-i/numFreqRanges 0.25 1/numFreqRanges], ...
        'Style','edit','String',num2str(freqRangeList0{i}(1)),'FontSize',fontSizeSmall);
    hFreqRangeMax{i} = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-i/numFreqRanges 0.25 1/numFreqRanges], ...
        'Style','edit','String',num2str(freqRangeList0{i}(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Axis Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel5 = uipanel('Title','Axis Ranges','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.625 1-panelHeight 0.15 panelHeight]);
axisRangeList0{1} = [0 100]; axisRangeName{1} = 'Freq Lims (Hz)';
axisRangeList0{2} = [0 1]; axisRangeName{2} = 'YLims';
axisRangeList0{3} = [0 1]; axisRangeName{3} = 'cLims (topo)';

numAxisRanges = length(axisRangeList0);
hAxisRangeMin = cell(1,numAxisRanges);
hAxisRangeMax = cell(1,numAxisRanges);

for i=1:numAxisRanges
    uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 1-i/numAxisRanges 0.5 1/numAxisRanges],'Style','text','String',axisRangeName{i},'FontSize',fontSizeSmall);
    hAxisRangeMin{i} = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-i/numAxisRanges 0.25 1/numAxisRanges], ...
        'Style','edit','String',num2str(axisRangeList0{i}(1)),'FontSize',fontSizeSmall);
    hAxisRangeMax{i} = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-i/numAxisRanges 0.25 1/numAxisRanges], ...
        'Style','edit','String',num2str(axisRangeList0{i}(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Cutoff Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel6 = uipanel('Title','Cutoffs','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.775 1-panelHeight 0.1 panelHeight]);
cutoffList0 = [2 30]; cutoffNames = [{'Num Elecs'} {'Num Trials'}];

numCutoffRanges = length(cutoffList0);
hCutoffs = cell(1,numCutoffRanges);

for i=1:numCutoffRanges
    uicontrol('Parent',hPanel6,'Unit','Normalized','Position',[0 1-i/numCutoffRanges 0.5 1/numCutoffRanges],'Style','text','String',cutoffNames{i},'FontSize',fontSizeSmall);
    hCutoffs{i} = uicontrol('Parent',hPanel6,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-i/numCutoffRanges 0.5 1/numCutoffRanges], ...
        'Style','edit','String',num2str(cutoffList0(i)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Plot Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel7 = uipanel('Title','Plot','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.875 1-panelHeight 0.1 panelHeight]);

hUseMedianFlag = uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 2/3 1 1/3],'Style','togglebutton','String','Use Median','FontSize',fontSizeMedium);
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium,'Callback',{@rescale_Callback});
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0.5 1/3 0.5 1/3],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback});
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 0 1 1/3],'Style','pushbutton','String','plot','FontSize',fontSizeMedium,'Callback',{@plot_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrodeGroupList = getElectrodeGroupsConn(groupType,1,'actiCap64_UOL');
numGroups = length(electrodeGroupList);
hAllPlots.hConn1 = getPlotHandles(1,numGroups,[0.05 0.55 0.6 0.27],0.01,0.01,1);
hAllPlots.hConn2 = getPlotHandles(numFreqRanges,numGroups,[0.05 0.05 0.6 0.45],0.01,0.03,0);
hAllPlots.hTopoRef = getPlotHandles(1,2,[0.66 0.57 0.32 0.3],0.002,0.002,1);
hAllPlots.hTopo  = getPlotHandles(numFreqRanges,3,[0.66 0.05 0.2 0.45],0,0,1);
hAllPlots.hConn3 = getPlotHandles(numFreqRanges,1,[0.88 0.05 0.1 0.45],0.01,0.01,1);

connMethod = 'ppc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_Callback(~,~)

        %%%%%%%%%%%%%%%%%%%%% Get SubjectLists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        comparisonStr=comparisonList{get(hComparison,'val')};
        if strcmp(comparisonStr,'paired')
            pairedSubjectNameList = getPairedSubjectsBK1;            
            subjectNameLists{1} = pairedSubjectNameList(:,1);
            subjectNameLists{2} = pairedSubjectNameList(:,2);
            pairedDataFlag      = 1;
        else
            [~, meditatorList, controlList] = getGoodSubjectsBK1;
            subjectNameLists{1} = meditatorList;
            subjectNameLists{2} = controlList;
            pairedDataFlag      = 0;
        end

        % Sub-select Subjects based on Demographics
        [subjectNameList,~,~,ageListAllSub,genderListAllSub] = getDemographicDetails('BK1');
        % sub-select based on Gender
        genderStr=genderList{get(hGender,'Value')};
        maleSubjectNameList = subjectNameList(strcmpi(genderListAllSub, 'M'));
        femaleSubjectNameList = subjectNameList(strcmpi(genderListAllSub, 'F'));
        if  strcmp(genderStr,'male')
            subjectNameLists{1} = intersect(subjectNameLists{1},maleSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},maleSubjectNameList,'stable');
        elseif strcmp(genderStr,'female')
            subjectNameLists{1} = intersect(subjectNameLists{1},femaleSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},femaleSubjectNameList,'stable');
        end

        % sub-select based on Age
        ageStr=ageList{get(hAge,'Value')};
        youngSubjectNameList = subjectNameList(ageListAllSub<40);
        midSubjectNameList = subjectNameList(ageListAllSub>=40);

        if  strcmp(ageStr,'young')
            subjectNameLists{1} = intersect(subjectNameLists{1},youngSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},youngSubjectNameList,'stable');
        elseif strcmp(ageStr,'mid')
            subjectNameLists{1} = intersect(subjectNameLists{1},midSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},midSubjectNameList,'stable');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        protocolName = protocolNameList{get(hProtocol,'val')};
        analysisChoice = analysisChoiceList2{get(hAnalysisChoice,'val')};
        refElectrodes = refElectrodeList{get(hRefChoice,'val')};

        badEyeCondition = badEyeConditionList2{get(hBadEye,'val')};
        badTrialVersion = badTrialVersionList{get(hBadTrialVersion,'val')};

        % stRange = [0.25 1.25]; % hard coded for now

        freqRangeList = cell(1,numFreqRanges);
        for ii=1:numFreqRanges
            freqRangeList{ii} = [str2double(get(hFreqRangeMin{ii},'String')) str2double(get(hFreqRangeMax{ii},'String'))];
        end

        axisRangeList = cell(1,numAxisRanges);
        for ii=1:numAxisRanges
            axisRangeList{ii} = [str2double(get(hAxisRangeMin{ii},'String')) str2double(get(hAxisRangeMax{ii},'String'))];
        end

        cutoffList = zeros(1,numCutoffRanges);
        for ii=1:numCutoffRanges
            cutoffList(ii) = str2double(get(hCutoffs{ii},'String'));
        end

        useMedianFlag = get(hUseMedianFlag,'val');

        displayConnDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,1);

    end
    function cla_Callback(~,~)
        claGivenPlotHandle(hAllPlots.hConn1);
        claGivenPlotHandle(hAllPlots.hConn2);
        claGivenPlotHandle(hAllPlots.hTopoRef);
        claGivenPlotHandle(hAllPlots.hTopo);
        claGivenPlotHandle(hAllPlots.hConn3);

        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    cla(plotHandles(ii,j));
                end
            end
        end
    end
    function rescale_Callback(~,~)
        axisLims = [str2double(get(hAxisRangeMin{1},'String')) str2double(get(hAxisRangeMax{1},'String')) str2double(get(hAxisRangeMin{2},'String')) str2double(get(hAxisRangeMax{2},'String'))];
        cLims = [str2double(get(hAxisRangeMin{3},'String')) str2double(get(hAxisRangeMax{3},'String'))];

        rescaleGivenPlotHandle(hAllPlots.hConn1,axisLims);
        rescaleZGivenPlotHandle(hAllPlots.hTopo,cLims);

        function rescaleGivenPlotHandle(plotHandles,axisLims)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    axis(plotHandles(ii,j),axisLims);
                end
            end
        end
        function rescaleZGivenPlotHandle(plotHandles,cLims)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    clim(plotHandles(ii,j),cLims);
                end
            end
        end
    end

end