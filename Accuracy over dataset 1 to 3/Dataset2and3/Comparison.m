factor=1

%%
load("neurons_by_hand.mat");
load('neurons_reconstructed_max_2_e.mat');

s1=size(hand_e,1);
s2=numel(neurons_reconstructed_max_2);
guessed=zeros(1,s1);
taken=zeros(1,s2);
sizes=12;
TP=0;
FN=0;
FP=0;
minimum_distance=sizes.*factor;
for i=1:s1
    current_neuron=hand_e(i,[3 2 4]);
    best=3000;
    for j=1:s2
        if ~taken(j)
            current_guess=neurons_reconstructed_max_2(j).coords(1,1:3);
            current_distance=norm(current_neuron-current_guess);
            best=min([best current_distance]);
            if best==current_distance
                best_id=j;
            end
        end
    end
    if best<minimum_distance
        taken(best_id)=1;
        guessed(i)=1;
        TP=TP+1;
    end
end

disp(['TP = ' num2str(sum(guessed)) ';']);
disp(['FN = ' num2str(sum(guessed==0)) ';']);
disp(['FP = ' num2str(sum(taken==0)) ';']);


%%

load("neurons_by_hand.mat");
load('neurons_reconstructed_max_2_c.mat');

s1=size(hand_c,1);
s2=numel(neurons_reconstructed_max_2);
guessed=zeros(1,s1);
taken=zeros(1,s2);
sizes=21;

TP=0;
FN=0;
FP=0;
minimum_distance=sizes.*factor;


for k = 1:size(hand_c,1)
    for h = 1:numel(neurons_reconstructed_max_2)
        current_neuron=hand_c(k,[3 2 4]);
        current_guess=neurons_reconstructed_max_2(h).coords(1,1:3);
        C(h,k) = vecnorm(current_neuron-current_guess);
    end
end
C(C>minimum_distance)=Inf;
C=C.^2;
highest_cost=max(C(C~=Inf));
UnmatchedCost=double(highest_cost*1.05);
% lowest cost solution
UnmatchedCost(isempty(UnmatchedCost))=0;
links=matchpairs(C,minimum_distance);
links=links(:,[2 1]);
for i=1:s1
    current_neuron=hand_c(i,[3 2 4]);
    best=3000;
    for j=1:s2
        if ~taken(j)
            current_guess=neurons_reconstructed_max_2(j).coords(1,1:3);
            current_distance=norm(current_neuron-current_guess);
            best=min([best current_distance]);
            if best==current_distance
                best_id=j;
            end
        end
    end
    if best<minimum_distance
        taken(best_id)=1;
        guessed(i)=1;
        TP=TP+1;
    end
end
disp(' ');
disp(['TP = ' num2str(sum(guessed)) ';'])
disp(['FN = ' num2str(sum(guessed==0)) ';'])
disp(['FP = ' num2str(sum(taken==0)) ';'])







