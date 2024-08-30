%% feature_definition

%% 
% 5 features selected by visual inspection for each subject
f1=[9 9;9,10;9,20;9 8; 1 9];
f2=[13 5; 13 6; 15 5;13 4; 15 6];
f3=[11 6; 7 6; 9 12;9 11; 3 6];
f4=[11 5; 7 4; 7 5; 11 4;2 5];
f5=[7 6; 12 6;2 6; 7 5;3 6];
f6=[8 5; 8 6; 7 6;7 5; 11 5];
f7=[2 12; 4 12; 14 6; 13 6; 11 5];
f8=[12 5; 16 5; 8 5;16 6;10 5];
% arranged within a structure
feature_selected=struct('s1',f1,'s2',f2,'s3',f3,'s4',f4,'s5',f5,'s6',f6,'s7',f7,'s8',f8);
subject = {'s1','s2','s3','s4','s5','s6','s7','s8'};
%Definition of a basic matrix for visual representation
base_mat=-1*ones(16,23);
%creation of a structure: F_5; containing for each subject a matrix for visual representation of the chosen features
F_5=struct;
for zi=1:8    %for each subject
    curr_mat=base_mat;
    index=[feature_selected.(subject{zi})];  %feature extraction from the previously created structure
    for ii=1:3   %selection of the first 3 chosen features
        r=index(ii,1);
        c=index(ii,2);
        curr_mat(r,c)=0;
    end
    for ii=4:5   %selection of remaining features
        r=index(ii,1);
        c=index(ii,2);
        curr_mat(r,c)=+1;
    end 
    [F_5.(subject{zi})] = curr_mat;  %inserimento all'interno della struttura F_5 
    base_mat=-1*ones(16,23); %Definition of a basic matrix for visual representation
end
%saving structure F_5 in feature_selected.mat
save('feature_selected.mat','F_5')