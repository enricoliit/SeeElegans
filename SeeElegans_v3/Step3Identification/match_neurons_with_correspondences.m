function [tpr, tot] = match_neurons_with_correspondences(neurons_identified, correspondences)
    % Iterate through each neuron name
    tot = length(neurons_identified.name);
    tp = 0;
    for i = 1:length(neurons_identified.name)
        currentName = neurons_identified.name{i}{1};
        % Search through the correspondences
        for j = 1:size(correspondences, 1)
            if ismember(currentName, correspondences{j, 2})
                % If a match is found, print the name and associated ids
                neuronId = neurons_identified.id(i);
                correspondenceId = correspondences{j, 1};
                if neuronId == correspondenceId
                    tp = tp + 1;
                    fprintf('<strong>Match found: %s - Neuron ID: %d, Correspondence ID: %d</strong>\n', ...
                        currentName, neuronId, correspondenceId);
                else
                    fprintf('Match found: %s - Neuron ID: %d, Correspondence ID: %d\n', ...
                        currentName, neuronId, correspondenceId);
                end
            end
        end
    end
    tpr = tp/tot;
end