clearvars    

siteCodes = {'MC','GC','DT'};
% get dataa and effort
[effort,latLongs,errRates] = gofmx_dates3;
dayTimeSeriesTable = readtable('G:\GOM_Spermwhale\TPWS\MC_GC_DT_binsize011000_Group_density2_Pm.xls');
spCode = 'Pm';
outDir = sprintf('E:\\NASData\\ModelData\\%s',spCode);
% startTime = 734274; %<- May 16th, 2010, a Sunday
startTime = 734276; %<- May 18th, 2010, a Tuesday

endTime = datenum([2014,06,01,00,00,00]);
binInt = datenum([0,0,7,0,0,0]);

dateCol = table2array(dayTimeSeriesTable(:,1));
dateVec = datenum(x2mdate(dateCol));
pOccur = [];
siteVec = dayTimeSeriesTable.Site;

for siteNum = 1:length(siteCodes)

    siteIndices = find(strcmp(siteVec,siteCodes{siteNum})>0);
    
    % check for bins outside effort, and set times to NaN so they will be
    % excluded.
    myEffort = effort{siteNum};
    dailyDensities = dayTimeSeriesTable.meanDensity(siteIndices);
    allMonitoredDays = dateVec(siteIndices);
    inBin = [];
    for iE = 1:size(myEffort,1)
        [~,binTemp] = histc(allMonitoredDays,myEffort(iE,:));
        inBin = [inBin;find(binTemp>0)];
        binTemp = [];
    end
    
    outsideEffort = setdiff(1:size(allMonitoredDays,1),inBin);
    if ~isempty(outsideEffort)
        allMonitoredDays(outsideEffort,:) = NaN;
    end
    
    weekVec = startTime:binInt:endTime+1;

    positiveDays = allMonitoredDays(dailyDensities>0);
   
    [nPosWeeks,binPosWeeks] = histc(positiveDays,weekVec);
    [nWeeks,binWeeks] = histc(allMonitoredDays,weekVec);

    pOccur(:,siteNum) = nPosWeeks./nWeeks;
end
%%
cellOut = [cellstr(datestr(weekVec','yyyy-mm-dd')),cellstr(num2str(pOccur,'%.4f,%.4f,%.4f'))];
fileID = fopen(fullfile(outDir,sprintf('ALLSITES_weeklyPOccurrence_%s_jahStart.csv',...
     spCode)),'w');
myHeaders = ['dateStr',siteCodes];
fprintf(fileID,'%s,%s,%s,%s,%s,%s\n',myHeaders{:});
formatSpec = '%s,%s\n';
[nrows,ncols] = size(cellOut);
for row = 1:nrows
    fprintf(fileID,formatSpec,cellOut{row,:});
end
fclose(fileID);
% %  xlswrite(fullfile(outDir,sprintf('ALLSITES_binsize%s_%s_density_jahStart.xls',...
% %        num2str(datevec(binInt)')',columnOrder{spWrite})), [colHeaders;mat3XLScellALL])