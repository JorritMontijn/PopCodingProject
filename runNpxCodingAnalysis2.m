%% aim
%{
"distinct subnetworks in mouse visual system for visual and non-visual signals"

At single trial basis, split which info information limiting noise
(projected unto f') and non-limiting noise (orthogonal directions). Then
look at correlation of those noise types between areas. E.g., lgn diff
corrs correlate with v1, but SC non-diff corrs correlate with v1


Introduce with: spike count correlations are generally low, but important.
Show how to decompose noise correlations on single trials and that
orth>para. Show that inter-areal correlations of shuffles are 0, then show
real data: all red. Trial-to -trial variability in MD coding is highly
correlated despite pairwise noise correlations being low

%}
%% define qualifying areas
cellUseAreas = {...
	'Subiculum',...
	'Posterior complex of the thalamus',...
	'Anterior pretectal nucleus',...
	'Superior colliculus',...
	'Lateral posterior nucleus of the thalamus',...
	'Nucleus of the optic tract',...
	'Dorsal part of the lateral geniculate complex',...
	'Primary visual area',...
	'posteromedial visual area',...
	'Anteromedial visual area',...
	'Retrosplenial area',...
	'Anterior area',...
	};
strFigDir = 'F:\Data\Results\NpxDims\';


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating');
end

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
cellPerfLR = cell(1,intAreas);
cellCorr_MuPara = cell(1,intAreas);
cellCorr_MuOrth = cell(1,intAreas);
cellCorr_MuParaS = cell(1,intAreas);
cellCorr_MuOrthS = cell(1,intAreas);
cellPara = cell(1,intAreas);
cellParaS = cell(1,intAreas);
cellOrth = cell(1,intAreas);
cellOrthS = cell(1,intAreas);
		
%% go through recordings
tic
for intRec=1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Rec
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Rec}));
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellStim(cellfun(@(x) x.structEP.intStimTypes,sThisRec.cellStim) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellStim);
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecOrientation = structStim.Orientation;
	
	%ensure repetitions are full
	[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = label2idx(vecOrientation);
	intRepNum = min(vecCounts);
	indRemove = vecRepetition>intRepNum;
	vecOrientation(indRemove) = [];
	vecStimOnTime(indRemove) = [];
	vecStimOffTime(indRemove) = [];
	
	%remove neurons from other recs
	indQualifyingNeurons = contains({sAggNeuron.Rec},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	
	%% select area 1
	for intArea1=1:intAreas
		strArea1 = cellUseAreas{intArea1};
		indArea1Neurons = contains({sUseNeuron.Area},strArea1,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% get spike times
		cellSpikeTimes1 = {sArea1Neurons.SpikeTimes};
		matSpikeCountsArea1 = getSpikeCounts(cellSpikeTimes1,vecStimOnTime,vecStimOffTime);
		matSpikeCountsShuffled1 = matSpikeCountsArea1(:,randperm(size(matSpikeCountsArea1,2)));
		[intNeurons1,intTrials] = size(matSpikeCountsArea1);
		[vecUniqueOris,dummy,vecOriIdx] = unique(vecOrientation);
		intTrials = numel(vecOrientation);
		intOriNum = numel(vecUniqueOris);
		
		%% compare parallel/orthogonal % real vs shuffled
		%if coding of ori is area's prime task, it ought to minimize
		%parallel variability (V1), if non-visual area it shouldn't
		%{
		%does it depend on # of neurons?
		intIters = 100;
		matMeanParaPerc = zeros(intIters,intNeurons1);
		matMeanParaPerc_S = zeros(intIters,intNeurons1);
		parfor intIter=1:intIters
			for intPopSize=2:intNeurons1
				%real
				[vecNoiseParallel1,vecNoiseOrthogonal1,vecNoiseTotal1] = doCrossValidatedNoiseParaOrtho(matSpikeCountsArea1(randperm(intNeurons1,intPopSize),:),vecOrientation,intTypeCV);
				matMeanParaPerc(intIter,intPopSize) = mean(vecNoiseParallel1 ./ (vecNoiseParallel1 + vecNoiseOrthogonal1));
				
				%shuffled
				[vecNoiseParallel1_S,vecNoiseOrthogonal1_S,vecNoiseTotal1_S] = doCrossValidatedNoiseParaOrtho(matSpikeCountsShuffled1(randperm(intNeurons1,intPopSize),:),vecOrientation,intTypeCV);
				matMeanParaPerc_S(intIter,intPopSize) = mean(vecNoiseParallel1_S ./ (vecNoiseParallel1_S + vecNoiseOrthogonal1_S));
			end
		end
		
		%% plot
		figure
		errorbar(1:intNeurons1,mean(matMeanParaPerc,1),std(matMeanParaPerc,[],1))
		hold on
		errorbar(1:intNeurons1,mean(matMeanParaPerc_S,1),std(matMeanParaPerc_S,[],1))
		hold off
		title(sprintf(strArea1))
		drawnow;
		%}
		
		%% is total variability correlated with population mean activity?
		intTypeCV = 2;
		[vecNoiseParallelCV,vecNoiseOrthogonalCV,vecNoiseTotalCV] = doCrossValidatedNoiseParaOrtho(matSpikeCountsArea1,vecOrientation,intTypeCV);
		vecMeanR = mean(matSpikeCountsArea1,1);
		
		dblR_MuTot = corr(vecMeanR',vecNoiseTotalCV');
		dblR_MuPara = corr(vecMeanR',vecNoiseParallelCV');
		dblR_MuOrth = corr(vecMeanR',vecNoiseOrthogonalCV');
		
		[vecNoiseParallelCV_S,vecNoiseOrthogonalCV_S,vecNoiseTotalCV_S] = doCrossValidatedNoiseParaOrtho(matSpikeCountsShuffled1,vecOrientation,intTypeCV);
		vecMeanR_S = mean(matSpikeCountsShuffled1,1);
		
		dblR_MuTot_S = corr(vecMeanR_S',vecNoiseTotalCV_S');
		dblR_MuPara_S = corr(vecMeanR_S',vecNoiseParallelCV_S');
		dblR_MuOrth_S = corr(vecMeanR_S',vecNoiseOrthogonalCV_S');
		
		%% is para/ortho separation correlated with decoding?
		%real
		%%{
		intStimTypes = numel(unique(vecOrientation));
		dblRandP = 1/intStimTypes;
		intTrials = size(matSpikeCountsArea1,2);
		dblLambda = intTrials/intStimTypes;
		[dblPerformanceLR,vecDecodedIndexCV,dummy1,matWeights,dummy2,matConfusionLR] = doCrossValidatedDecodingLR(matSpikeCountsArea1,vecOrientation,intTypeCV,dblLambda);
		[phat,vecCI] = binofit(dblPerformanceLR*intTrials,intTrials);
		
		%shuffled
		%[dblPerformanceLR_S,vecDecodedIndexCV,dummy1,matWeights,dummy2,matConfusionLR_S] = doCrossValidatedDecodingLR(matSpikeCountsShuffled1,vecOrientation,intTypeCV,dblLambda);
		%compare
		%[p,z] = bino2test(dblPerformanceLR*intTrials,intTrials,dblPerformanceLR_S*intTrials,intTrials);
		%}
		%% save data
		cellPerfLR{intArea1} = [cellPerfLR{intArea1} dblPerformanceLR];
		cellPara{intArea1} = [cellPara{intArea1} mean(vecNoiseParallelCV)];
		cellParaS{intArea1} = [cellParaS{intArea1} mean(vecNoiseParallelCV_S)];
		cellOrth{intArea1} = [cellOrth{intArea1} mean(vecNoiseOrthogonalCV)];
		cellOrthS{intArea1} = [cellOrthS{intArea1} mean(vecNoiseOrthogonalCV_S)];
		
		cellCorr_MuPara{intArea1} = [cellCorr_MuPara{intArea1} dblR_MuPara];
		cellCorr_MuOrth{intArea1} = [cellCorr_MuOrth{intArea1} dblR_MuOrth];
		cellCorr_MuParaS{intArea1} = [cellCorr_MuParaS{intArea1} dblR_MuPara_S];
		cellCorr_MuOrthS{intArea1} = [cellCorr_MuOrthS{intArea1} dblR_MuOrth_S];
	end
end
toc

%% plot
vecMeanDecPerf = cellfun(@mean,cellPerfLR);
[dummy,vecReorderAreas] = sort(vecMeanDecPerf,'descend');

%cellPerfLR = cell(1,intAreas);
%cellCorr_MuPara = cell(1,intAreas);
%cellCorr_MuOrth = cell(1,intAreas);
%cellCorr_MuParaS = cell(1,intAreas);
%cellCorr_MuOrthS = cell(1,intAreas);
figure
subplot(2,3,1:3)
%scatter(cell2vec(cellCorr_MuPara)-cell2vec(cellCorr_MuParaS),cell2vec(cellCorr_MuOrth)-cell2vec(cellCorr_MuOrthS),[],cell2vec(cellPerfLR))
cla;
hold on
for intArea=1:intAreas
	intPlotArea = vecReorderAreas(intArea);
	intRecs = numel(cellCorr_MuPara{intPlotArea});
	%plot(repmat(intArea+[-0.5 0.5],[intRecs 1])',[cellCorr_MuPara{intPlotArea}; cellCorr_MuOrth{intPlotArea}])
	%plot(repmat(intArea+[-0.5 0.5],[intRecs 1])',[cellCorr_MuPara{intPlotArea} - cellCorr_MuOrth{intPlotArea}; cellCorr_MuParaS{intPlotArea} - cellCorr_MuOrthS{intPlotArea}])
	plot(repmat(intArea+[-0.2 0.2],[intRecs 1])',[cellOrth{intPlotArea} ./ cellPara{intPlotArea}; cellOrthS{intPlotArea} ./ cellParaS{intPlotArea}])
	
	%cellCorr_MuOrth{intArea}
	%cellCorr_MuOrthS{intArea}
end
hold off
set(gca,'xtick',1:intAreas,'xticklabel',cellUseAreas(vecReorderAreas))
ylabel('Fraction of variability orthogonal to f''')
xtickangle(15)
title('|Orth|/|Para| vs |Orth_S|/|Para_S|, ordered by decoding performance')
fixfig;
maxfig;

subplot(2,3,4)
scatter(cell2vec(cellOrth),cell2vec(cellPerfLR))
xlabel('f''-orthogonal magnitude')
ylabel('LR decoding performance')
fixfig;

subplot(2,3,5)
scatter(cell2vec(cellPara),cell2vec(cellPerfLR))
xlabel('f''-parallel magnitude')
ylabel('LR decoding performance')
fixfig;

subplot(2,3,6)
scatter(cell2vec(cellOrth) ./ cell2vec(cellPara),cell2vec(cellPerfLR))
xlabel('Orth/para fraction')
ylabel('LR decoding performance')
fixfig;

%% TO DO
%normalize magnitude for pop size
%compare with shuffled f' rather than shuffled trials
error('to do!')

%% save fig
strFigFile = sprintf('OrthParaMagnitudes_%s',getDate);
maxfig;
drawnow;
export_fig([strFigDir strFigFile '.tif']);
export_fig([strFigDir strFigFile '.pdf']);
		