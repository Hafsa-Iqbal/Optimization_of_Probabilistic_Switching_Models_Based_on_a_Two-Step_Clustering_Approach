function labels_in_level = level_info (HirCluster, inial_number_clusters)
labels_in_level{1,1} = 1:inial_number_clusters;
labels_in_level{1,2} = [];
current_level = 2;
list_of_levels = [1,2];
current_number_clusters = inial_number_clusters;
max_number_clusters(1) =  inial_number_clusters;
for i = 1:size(HirCluster,1)
    current_number_clusters = current_number_clusters + 1; 
    
    if HirCluster(i,1) < inial_number_clusters && HirCluster(i,2) < inial_number_clusters
        max_number_clusters(current_level) = current_number_clusters;
        labels_in_level{1,current_level} = [labels_in_level{1,current_level} current_number_clusters];

    else
        first_belong = find( ((HirCluster(i,1) > max_number_clusters) == 0), 1);
        second_belong = find( ((HirCluster(i,2) > max_number_clusters) == 0), 1);
        
        current_level = max([first_belong, second_belong]) + 1;
        if  ~isempty(setdiff(current_level, list_of_levels))
            labels_in_level{1,current_level} = [];
            list_of_levels = [list_of_levels current_level];
        end
        
        max_number_clusters(current_level) = current_number_clusters;
        labels_in_level{1,current_level} = [labels_in_level{1,current_level} current_number_clusters]; 
    end 
end