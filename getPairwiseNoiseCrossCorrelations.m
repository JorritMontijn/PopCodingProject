function matNoiseCrossCorrs = getPairwiseNoiseCrossCorrelations(matSpikeCountsArea1,matSpikeCountsArea2,vecOrientation)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
		
	[intNeurons1,intTrials] = size(matSpikeCountsArea1);
	[intNeurons2,intTrials] = size(matSpikeCountsArea2);
	[vecUniqueOris,dummy,vecOriIdx] = unique(vecOrientation);
	intTrials = numel(vecOrientation);
	intOriNum = numel(vecUniqueOris);
	
	matNoiseCrossCorrs = nan(intNeurons1,intNeurons2,intOriNum);
	for intOriIdx=1:intOriNum
		%stim
		vecTrialsOri = find(vecOriIdx==intOriIdx);
		intRepNum = numel(vecTrialsOri);
		matAct1 = matSpikeCountsArea1(:,vecTrialsOri);
		matAct2 = matSpikeCountsArea2(:,vecTrialsOri);
		
		%noise corrs
		matNoiseCrossCorrs(:,:,intOriIdx) = corr(matAct1',matAct2');
	end
end

