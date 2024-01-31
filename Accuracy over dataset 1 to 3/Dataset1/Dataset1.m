clear all

%%
load("neurons_by_hand_1.mat");
load("neurons_reconstructed_max_2.mat");
load('Dataset1.mat')

cdim=1;
swscd=size(wholestack,2);
for i=1:numel(neurons_reconstructed_max_2)
    neurons_reconstructed_max_2(i).coords(:,cdim)=-1.*neurons_reconstructed_max_2(i).coords(:,cdim)+swscd;
end

factor=1;
s1=size(neurons_by_hand,1);
hande_num=s1;
s2=numel(neurons_reconstructed_max_2);
guessed=zeros(1,s1);
taken=zeros(1,s2);
sizes=12;
TP=0;
FN=0;
FP=0;
minimum_distance=sizes.*factor;
for i=1:s1
    current_neuron=neurons_by_hand(i,[3 2 4]);
    best=3000;
    for j=1:s2
        if ~taken(j)
            current_guess=neurons_reconstructed_max_2(j).coords(1,[1 2 3]);
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
TP=sum(guessed);
FN=sum(guessed==0);
FP=sum(taken==0);
disp(['TP = ' num2str(TP) ';']);
disp(['FN = ' num2str(FN) ';']);
disp(['FP = ' num2str(FP) ';']);

% disp(['TP = ' num2str(100.*TP/hande_num) ';']);
% disp(['FN = ' num2str(100.*FN/hande_num) ';']);
% disp(['FP = ' num2str(100.*FP/hande_num) ';']);

ntable_ds1(1,1:4)={'Dataset1 (TP,FN,FP)',num2str(100.*TP./hande_num),num2str(100.*FN./hande_num,4),num2str(100.*FP./hande_num,4)};


figure
I=sum(wholestack(:,end:-1:1,:,1),3)';
imagesc(log(I))
colormap("gray")
axis equal
axis off

hold on
for i=1:s1
    colr='b';
    if guessed(i)==0
        colr='y';
    end
    hold on
    viscircles([neurons_by_hand(i,2),neurons_by_hand(i,3)],6,'LineWidth',0.01,'Color',colr);
end

for i=1:s2
    colr='g';
    if taken(i)==0
        colr='r';
    end
    hold on
    viscircles([neurons_reconstructed_max_2(i).coords(:,2), neurons_reconstructed_max_2(i).coords(:,1)],6,'LineWidth',0.01,'Color',colr);
end

savefig('Dataset1.fig');
print(gcf,'Dataset1.png','-dpng','-r300');
print(gcf,'Dataset1.svg','-dsvg');
