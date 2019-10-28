List_directory = dir("IFNG_project/Image_analysis_kinetic/Data/") ;
List_directory = struct2table(List_directory) ; 
List_directory = List_directory.name ; 
List_directory = (List_directory(3:size(List_directory,1)));
Output_directory =  "IFNG_project/Image_analysis_kinetic/Results/"

parfor number_lame=1:size(List_directory,1) 
    
    %mkdir(strcat(Output_directory,List_directory{number_lame}));
    List_samples = dir(strcat("/full/path/to/home/",List_directory{number_lame}));
    List_samples = struct2table(List_samples) ; 
    List_samples = List_samples.name ; 
    List_samples = (List_samples(3:size(List_samples,1)));
    
    for k=1:size(List_samples)
        
        Temp_image = LoadImage(char(strcat("/full/path/to/home//",List_directory{number_lame},"/",num2str(k),"/")));
        
        %%Extracting the nuclei directly
        Nuclei = Nuclei_identification(imadjust(Temp_image(:,:,1,2)));

        %%Extracting STAT1 signal
        se = strel('disk',4);

        STAT = imbinarize(imadjust(Temp_image(:,:,1,1)));
        STAT = imclose(STAT,se);
        STAT = imfill(STAT, 'holes');

        %Extracting Cell enveloppe
        Cell_enveloppe = Nuclei_identification(STAT);
        
        %Extracting score

         Overlapping_score = [];

         STAT_signal = Temp_image(:,:,1,1);
         Nucleus_signal = Temp_image(:,:,1,2);

        for i=1:size(Cell_enveloppe)
    
            temp_enveloppe = Cell_enveloppe{i};  
            temp_enveloppe_ROI = poly2mask(temp_enveloppe(:,1),temp_enveloppe(:,2),size(Temp_image,1),size(Temp_image,2));
            global_STAT_signal = STAT_signal(temp_enveloppe_ROI);
            temp_nuclei_signal = Nucleus_signal(temp_enveloppe_ROI);
            temp_score = corr(global_STAT_signal,temp_nuclei_signal,"Type","P");
            Overlapping_score = [Overlapping_score temp_score];
            
        end

        writetable(array2table(Overlapping_score.'), strcat(Output_directory,List_directory{number_lame},"_",num2str(k),".txt"),"delimiter","\t")

    end

end


