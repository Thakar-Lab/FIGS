function [ Results_HypGeoTest ] = HypGeoTest( Results_FCM, PW )

% Developed by Atif Khan, March 21, 2017. 

% Perform Hypergeometric Test 
for i=1:length(Results_FCM.FCM_size)
       for j=1:length(PW.PW_Size)
       FCM_Match_1 = intersect(Results_FCM.Fuzzy_Genes(1:Results_FCM.FCM_size(i),i),PW.PWs{1,j});
       size_of_match(j,1:2)=size(FCM_Match_1);
       FCM_Match(i,j)=size_of_match(j,1);   % This gives the common genes in bw FCM clusters and each PW
      
       % Hypergeometric test processing
       x_act=intersect(Results_FCM.Fuzzy_Genes(1:Results_FCM.FCM_size(i),i),PW.PWs{1,j});
       x=size(x_act);
       K_act=intersect(PW.PWs{1,j},Results_FCM.Genes);
       K=size(K_act);
       N=Results_FCM.FCM_size(i); % Size of fuzzy cluster
       hyp_geo = hygecdf(x(1,1),length(Results_FCM.Genes),K(1,1),N,'upper');
       hyp_geo_test(i,j)=hyp_geo(1,1);
       end;
       hyp_geo_test_size=size(hyp_geo_test);
end;

hyp_geo_test(hyp_geo_test == 0) = NaN;  % Convert indices of 0 probability to NaN. 

 % Find the total number of significantly enriched PWs for each FCM cluster
 for i=1:length(Results_FCM.FCM_size)
 Sig_Clus(i) = numel(hyp_geo_test(hyp_geo_test(i,1:hyp_geo_test_size(1,2)) <0.05));
 end;

FCM_Match=transpose(FCM_Match);
hyp_geo_test=transpose(hyp_geo_test); 

Results_HypGeoTest.FCM_Match=FCM_Match;
Results_HypGeoTest.hyp_geo_test=hyp_geo_test;
Results_HypGeoTest.Sig_Clus=Sig_Clus;
% Save the results
save('Results_HypGeoTest', 'Results_HypGeoTest');

end

