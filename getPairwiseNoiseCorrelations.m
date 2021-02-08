function [matNoiseCorrs,matSignalCorrs] = getPairwiseNoiseCorrelations(matSpikeCounts,vecOrientation)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
		
	[intNeurons,intTrials] = size(matSpikeCounts);
	[vecUniqueOris,dummy,vecOriIdx] = unique(vecOrientation);
	intTrials = numel(vecOrientation);
	intOriNum = numel(vecUniqueOris);
	
	matNoiseCorrs = nan(intNeurons,intNeurons,intOriNum);
	matSignalCorrs = nan(intNeurons,intNeurons);
	matMeanResp = nan(intNeurons,intOriNum);
	for intOriIdx=1:intOriNum
		%stim
		vecTrialsOri = find(vecOriIdx==intOriIdx);
		intRepNum = numel(vecTrialsOri);
		matAct = matSpikeCounts(:,vecTrialsOri);
		matMeanResp(:,intOriIdx) = mean(matAct,2);
		
		%noise corrs
		matNoiseCorrs(:,:,intOriIdx) = corr(matAct');
	end
	matSignalCorrs = corr(matMeanResp');
end

