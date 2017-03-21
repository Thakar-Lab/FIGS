function [ output_args ] = Display_Circular_Graph(Results_FCM,PW,Results_HypGeoTest)

% Developed by Atif Khan, March 21, 2017. 

% Sort the p-values from hypergeometric test 
[a b]=sort(Results_HypGeoTest.hyp_geo_test,'ascend');
p=0.05;   % p-value significance level 

% Get the names of most significantly enriched pathways for each FCM
% cluster 
for i=1:length(PW.PW_Size)
for j=1:length(Results_FCM.FCM_size)
if(a(i,j)<p)
Most_Sig_PW_Enriched{i,j}= PW.PW_Name(b(i,j),1);  
else
Most_Sig_PW_Enriched{i,j}=[];    % empty if no significant pathway found 
end;
end;
end;

Label=Most_Sig_PW_Enriched(1,:);            %# Get the list of only most significant pathway
emptyIndex = cellfun(@isempty,Label);       %# Find indices of empty cells
Label(emptyIndex) = {' \color{blue}'};      %# Fill empty cells with 0

    % Loop to concatenate cluster number and the corresponding PW enriched
    for i=1:length(Results_FCM.FCM_size)
    Col_Name(1,i)=cellstr(['C',num2str(i),'  ',char(Label{1,i})]);
    end;
    Col_Name=regexprep(Col_Name,'_',' ');  % replace '_' with space for better representation

% Make circular graph that displays overlap among FCM clusters labelled
% with most significantly enriched pathwyas 
figure;
myColorMap = lines(length(Results_FCM.Cluss_Clus_Overlap));
circularGraph(Results_FCM.Cluss_Clus_Overlap,'Colormap',myColorMap,'Label',Col_Name);
set(gcf,'NumberTitle','off');
set(gcf,'Name',['Overlap between FIGS gene-sets     ', char(169),'Thakar-Lab   @URMC']);
set(findobj(gcf,'Type','text'),'FontSize',10);
%set(gcf, 'MenuBar', 'None'); % Allow user to copy/save resutls  


 X = Results_HypGeoTest.hyp_geo_test;
 %# convert matrix of numbers to cell array of strings (right aligned)
 XX = reshape(strtrim(cellstr(num2str(X(:)))), size(X));
 for i=1:length(Results_FCM.FCM_size)
 %# find cells matching condition (p-values <0.05)
 idx = ( X(:,i) < p );
 %# use HTML to style these cells in red color 
 XX(idx,i) = strcat(...
 '<html><span style="color: #FF0000; font-weight: bold;">', ...
 XX(idx,i), ...
 '</span></html>');
 end;
    % Name the columns of table with FCM cluster names  
    for i=1:length(Results_FCM.FCM_size)
    Table_Col_Name(i)=cellstr(['C',num2str(i)]);
    end;
    
% Make table that dispalys composition of each FCM cluster and also gives
% the p-values of enrichment with corresponding pathways    
    
f = figure;
h = uitable(f(1), 'Units','normalized', 'ColumnWidth','auto','Position',[0 0.5 1 .5]);  % First Table of P Values 
set(h(1), 'Data',XX,'rowname',PW.PW_Name,'columnname',Table_Col_Name)
% Second Table of genes 
h(2) = uitable(f,'columnname',Table_Col_Name, 'ColumnWidth','auto',...
                     'data',Results_FCM.Fuzzy_Genes, ...
                     'units','normalized', ...
                     'pos',[0 0 1 0.5]);
   
% Set the title of the GUI                  
set(gcf,'NumberTitle','off'); 
set(gcf,'Name',['FIGS: Genes sets    ', char(169),'Thakar-Lab   @URMC']);
%set(gcf, 'MenuBar', 'None');  % Allow user to copy/save resutls  

end

