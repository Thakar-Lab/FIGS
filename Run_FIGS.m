function Run_FIGS

%    Developed by Atif Khan, March 21, 2017. 
%    #### RUN_FIGS Summary ####
%    Load data through GUI
%    GUI will help you locate and upload your data file 
%    First Column of data file must contain the gene names 
%    First Row may contain gene names (optional) 
%    Loading Pathways for enrichment is optional. You will be provided an
%    option to upload your own pathway file for enrichment 
%    Each row of Pathway file will be treated as indivisual pathway. 
%    First element of each row Pathway file will be considered as pathway name.  

% Get the data file from user specified location 
    [filename, pathname]=uigetfile;      % open GUI for lloading file 
    if ~ischar(filename)
    errordlg('No File Uploaded. Re-run program and upload file','File Error');
    return;
    end;
    Path=strcat(pathname,filename);   % Get full path and file name 
    Data_all=readtable(Path);   % Store data into table 
    Genes=Data_all{:,1};    % Get all genes names from the table 
    Data=table2array(Data_all(:,2:end));    % Get data (pair-wise distance between the genes)
    clear Path filename  pathname Data_all   % Clear non resuable variables 

    
% Display the size of the data file after successful uplaod. 
    Gene_Size=length(Genes);
    Message=['*** Successfully uploaded ' num2str(Gene_Size)  ' genes and ' num2str(Gene_Size) ' x ', num2str(Gene_Size), ' Mutual Information Matrix ***  '];
    fprintf(Message);
    fprintf('\n');
 
% Ask users if they want to provide their pathways for enrichment analysis. 
    PW_choice = questdlg('Would you like a upload your pathways for enrichments? (You can upload .txt file where each row corresponds to a pathway. First element of each row will be treated as pathway name. Please assign brief pathway names for better representation.  ', ...
        'Pathways for Enrichment', ...
        'Yes','No','No');
% Handle response
    switch PW_choice
        case 'Yes'     % Get the data for pathways 
        [filename, pathname]=uigetfile;      % open GUI for lloading file 
        Path=strcat(pathname,filename);   % Get full path and file name 
        PW_all=readtable(Path,'ReadVariableNames',false);   % Store data into table 
        PW.PW_Name=PW_all{:,1};    % Get all genes names from the table 
        size_PW=size(PW_all);

        for i=1:size_PW(1,1);
        PW_temp=table2cell(PW_all(i,2:end));
        %# find empty cells
        emptyCells = cellfun(@isempty,PW_temp);
        %# remove empty cells
        PW_temp(emptyCells) = [];
        PW.PWs{1,i}=PW_temp;
        PW.PW_Size(1,i)=length(PW_temp);
        end;
        clear Path filename pathname PW_temp emptyCells % Clear non resuable variables     

        case 'No'
            disp('You decided not to provide your pathways for enrichment.')        
    end

% Pre-Processing of Gene-Gene Mutal Information Matrix 
    Norm_Data=mat2gray(Data);    % Normalize data to 0-1 usign min-max normalization 
    Processed_Data=1-(Norm_Data);    % Convert similarity measure to dissimilarity measure 

% Begining of Fuzzy C-Means Clustering Processing 
% Get the number of clusters and fuzziness between the clusters. 
    prompt = {'Enter number of clusters:','Enter fuzziness between the clusters (>1):'};
    dlg_title = 'Fuzzy C-Means Parameters';
    num_lines = 1;
    defaultans = {'10','1.1'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    no_of_clusters =answer(1,1);
    no_of_clusters=str2double(no_of_clusters);
    fuzziness=str2double(answer(2,1));

% Get the clsuter association criteria. 
    Clus_ass_choice = menu('Choose cluster association criteria, MV=Membership Value','Hard (max MV)','(Mean+Std) of MV','Mean(Max MV)/2', 'Enter your own threshold');
    if Clus_ass_choice==4
    prompt = {'Enter your own membership value threshold between (0-1)'};
    dlg_title = 'Cluster Association Threshold';
    num_lines = 1;
    defaultans = {'0.5'};
    clus_user_MV = inputdlg(prompt,dlg_title,num_lines,defaultans);
    Clus_ass_choice=str2double(clus_user_MV(1,1));
    else
    end; 

% Get the Ward hierarchical clusters from the dissimilarity matrix
    Ward_link = linkage(Processed_Data,'ward','euclidean');
    Ward_clusters = cluster(Ward_link,'maxclust', no_of_clusters);
    size_Data=size(Processed_Data);
    for i=1:no_of_clusters
    index1 = find(Ward_clusters == i);
    indexed = index1;
    index{1,i}(1,:) = indexed;
    index{1,i}(2,:)= Ward_clusters(indexed,1);
    [a_sorted, a_order] = sort(index{1,i}(2,:),'descend');
    index{1,i}(1,:) = index{1,i}(1,a_order);
    index{1,i}(2,:) = a_sorted;
    Ward_genes{i}= Genes(index{1,i}(1,:),1);
    end;
    clear index index1 indexed a_sorted a_order;  
    
% Get the centroids from Ward clusters. 
    for i=1:no_of_clusters
    index1 = find(Ward_clusters == i);
    indexed = index1;
    Ward_Transformed(i,:)=0;
    Ward_Transformed(i,index1)=1;
    end;
    
% Initialize fuzzy partitions by using Ward based centroids 
     col_sum = sum(Ward_Transformed);
     Ward_centers = Ward_Transformed./col_sum(ones(no_of_clusters, 1), :);
     
% Perform Fuzzy C-Means Clusterign on Ward Initialzed Centers 
     options=[fuzziness NaN NaN 1];
     [center, U, obj_fcm] = FIGS_FCM(Processed_Data, no_of_clusters, options, Ward_centers);
     Results_FCM = Process_FCM(U,Clus_ass_choice,Genes);
     Results_FCM.Genes=Genes;
     save('Results_FIGS', 'Results_FCM','Ward_genes');

% Get and allign all pathway names. 
      switch PW_choice
      case 'Yes'     % Check if user has provided pathways 
        s = struct('a',PW.PWs);
        for i=1:length(PW.PWs)  % Put all pathways in single object
        T=s(i).a;
        size_T=size(T);
        size_T=size_T(1,2);
        PW_Genes(i,1:size_T)=T;
        end;
        PW.PW_Genes=(PW_Genes)';

        ALL_PW_Genes = PW_Genes(:);
        ALL_PW_Genes = ALL_PW_Genes(~cellfun('isempty',ALL_PW_Genes));
        PW.Unique_Pathway_Genes = unique(ALL_PW_Genes);
        PW.Size_Uni_PW_Genes=size(PW.Unique_Pathway_Genes);
        clear ALL_PW_Genes T size_T PW_Genes s i ;
      
% Perform Hypergeometric test 
      Results_HypGeoTest = HypGeoTest(Results_FCM, PW);
% Dsiplay overp as an interactive circular graph
      Display_Circular_Graph(Results_FCM,PW,Results_HypGeoTest);
      end; 

    Message=['*** Successfully scored resutls in FIGS_Results.cvs file ***  '];
    fprintf(Message);
    fprintf('\n');

end

