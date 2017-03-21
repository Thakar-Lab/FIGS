function [ Results_FCM] = Process_FCM(U,cluster_association_user,Genes)

% Developed by Atif Khan, March 21, 2017. 
% Process FCM results to associate genes to the clusters based on user
% defined association criteria. 

size_U=size(U); % Get the number of genes and the clustes
Data_Genes = upper(Genes); % Convert all gene symbols to upper case

% Get the statistics for the fuzzy membership matrix 
sortedU=sort(U,'descend');
maxU=max(U);
meanU = mean(U);
stdU=std(U);
y=cluster_association_user; % Get user defined association criteria
    
if y==1;  % for case 1: max membership value, hard or crisp clusters. 
for i=1:size_U(1,1)
index1 = find(U(i,:) ==maxU);
indexed = index1;
index{1,i}(1,:) = indexed;
index{1,i}(2,:)= U(i,indexed);
[a_sorted, a_order] = sort(index{1,i}(2,:),'descend');
index{1,i}(1,:) = index{1,i}(1,a_order);
index{1,i}(2,:) = a_sorted;
genes_index{i}= Data_Genes(index{1,i}(1,:),1);
end;

elseif y==2;  % Case 2: associate genes to the clusters based on their membership value distribution
for i=1:size_U(1,1)
index1 = find(U(i,:) >=(meanU+stdU));
indexed = index1;
index{1,i}(1,:) = indexed;
index{1,i}(2,:)= U(i,indexed);
[a_sorted, a_order] = sort(index{1,i}(2,:),'descend');
index{1,i}(1,:) = index{1,i}(1,a_order);
index{1,i}(2,:) = a_sorted;
genes_index{i}= Data_Genes(index{1,i}(1,:),1);
end;

elseif y==3;  % Case 3: mean(max membership value))2
for i=1:size_U(1,1)
index1 = find(U(i,:) >= (mean(maxU)/2));
indexed = index1;
index{1,i}(1,:) = indexed;
index{1,i}(2,:)= U(i,indexed);
[a_sorted, a_order] = sort(index{1,i}(2,:),'descend');
index{1,i}(1,:) = index{1,i}(1,a_order);
index{1,i}(2,:) = a_sorted;
genes_index{i}= Data_Genes(index{1,i}(1,:),1);
end;

else    

user_threshold = cluster_association_user;  % Case 4: Use user defined fuzzy assiciation value

for i=1:size_U(1,1)
index1 = find(U(i,:) >= user_threshold);
indexed = index1;
index{1,i}(1,:) = indexed;
index{1,i}(2,:)= U(i,indexed);
[a_sorted, a_order] = sort(index{1,i}(2,:),'descend');
index{1,i}(1,:) = index{1,i}(1,a_order);
index{1,i}(2,:) = a_sorted;
genes_index{i}= Data_Genes(index{1,i}(1,:),1);
end;

end;

% Convert indexes of all the clusters to corresponding gene symbols
s = struct('a',genes_index);
for i=1:size_U(1,1)
T=s(i).a;
size_T=size(T);
size_T=size_T(1,1);
FCM_Genes(1:size_T,i)=T;
end;

% Get and store the size of each cluster
for i=1:size_U(1,1)
FCM_size(i) = sum(~cellfun('isempty',FCM_Genes(:,i)));
end;

% Find and store the most fuzzy gene
yourvector=FCM_Genes(:);
yourvector=yourvector(~cellfun('isempty',yourvector));
yourvector_unique=unique(yourvector);
yourvector_unique_size=size(yourvector_unique);

if isempty(yourvector_unique)
    Most_Fuzzy_Genes='No Fuzzy Gene Found'; % In case there is no fuzzy gene found. 
else 
    for i=1:yourvector_unique_size(1,1)
    matched=strmatch(yourvector_unique(i,1),yourvector(1:size(yourvector),1),'exact');
    counted(i,1:2)=size(matched);
    clear matched;
    end;
    yourvector_unique(:,2)=num2cell(counted(:,1));
    Most_Fuzzy_Genes=sortrows(yourvector_unique,-2);
end; 

% Process the overlap among the clusters  
for i=1:size_U(1,1)
    for j=1:size_U(1,1)
    Cluss_Match = intersect(FCM_Genes(1:FCM_size(1,i),i),FCM_Genes(1:FCM_size(1,j),j));
    size_clus=size(Cluss_Match);
    size_clus=size_clus(1,1);
    Cluss_Matched(i,j)=size_clus;
    clear size_clus;
    end;
end;
Cluss_Matched(logical(eye(size(Cluss_Matched)))) = 0;

% Store all the results in structure 
Results_FCM.Fuzzy_Genes=FCM_Genes;
Results_FCM.FCM_size=FCM_size;
Results_FCM.Most_Fuzzy_Genes=Most_Fuzzy_Genes;
Results_FCM.U_of_Genes=index;
Results_FCM.Cluss_Clus_Overlap=Cluss_Matched;
Results_FCM.genes_index=index;

% Save fuzzy gene-sets in excel file 
FCM_GeneSet_=FCM_Genes;
T = cell2table(FCM_GeneSet_);
writetable(T,'FIGS_Results.csv');
end







