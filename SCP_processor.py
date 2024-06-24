import pandas as pd
import numpy as np

class SCP_processor:
    def sumIDs(IDMatrix):
        """_summarize the ID matrix infor into ID summary_
        

        Args:
            IDMatrix (_type_): _protein or pepetides matrix_
            0 Symbol/Annotated Sequence 	run1 	run2 	run3 
            1 P023D12	MS2 	MBR 	NaN 
            2 P1222	NaN 	ID 	NaN 
        ID: means we don't know the ID mode

        Returns:
            _type_: _description_
                                        names  MS2_IDs  MBR_IDs  Total_IDs
    0            10ng_QC_1_channel2 Intensity      NaN      NaN       3650
    1            10ng_QC_2_channel1 Intensity      NaN      NaN       3604
    ....
        """
        # removes the columns that don't have ID data
        columns = [col for col in IDMatrix.columns if not any(
            substring in col for substring in [
                'Symbol', 'Annotated Sequence'])]
        #put each ID_Modes into a list
        returnNames = []
        MS2_ID = []
        MBR_ID = []
        total_ID = []
        for eachColumn in columns:
            MS2_ID.append(len(IDMatrix[eachColumn][IDMatrix[eachColumn] == "MS2"])) #PD differentiates
            MBR_ID.append(len(IDMatrix[eachColumn][IDMatrix[eachColumn] == "MBR"])) #PD differentiates
            total_ID_each = len(IDMatrix[eachColumn][IDMatrix[eachColumn] == "ID"]) #some don't so we count total directly
            if total_ID_each == 0: #otherwise we sum
                total_ID_each = len(IDMatrix[eachColumn][
                    IDMatrix[eachColumn] == "MS2"]) + len(IDMatrix[
                    eachColumn][IDMatrix[eachColumn] == "MBR"])
            total_ID.append(total_ID_each)

        return pd.DataFrame({'names': columns,
                            'MS2_IDs': MS2_ID,
                            'MBR_IDs': MBR_ID,
                            'Total_IDs': total_ID})

    def generate_column_from_name_mapping(columns, partial_column_name_mapping):
        #input is column names, and a dictionary with what you want each column (key) to be renamed to (value)
        column_name_mapping = {}
        for col in columns:
            for key, value in partial_column_name_mapping.items():
                if key in col: #in the case of PD, we are looking for a pattern within the column name
                    column_name_mapping[col] = value
                    break
        return column_name_mapping

    def generate_column_to_name_mapping(columns, partial_column_name_mapping):
        #input is column names, and a dictionary with what you want each column (key) to be renamed to (value)
        column_name_mapping = {}
        for col in columns:
            for key, value in partial_column_name_mapping.items():
                if key == col:  #after we get away from PD's weirdness, then we want exact matches,
                                #so we don't get for example 1-11 when we look for 1-1 for file Identifiers or filenames
                    column_name_mapping[col] = value
                    break
        return column_name_mapping


    def combine_IDs(all_matrix, MS2_matrix):
        # make IDs into MBR
        if "Annotated Sequence" in all_matrix.columns:
            name = "Annotated Sequence"
        elif "Symbol" in all_matrix.columns:  
            name = "Symbol"
        id_cols = all_matrix.columns.tolist()
        id_cols.remove(name)
        
        
        # for eachColumn in id_cols:
        #     if len(all_matrix[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2"))]) > 0:
        #         all_matrix.loc[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2")), [eachColumn]] = MS2_matrix[[eachColumn]]

        all_keys = pd.merge(all_matrix[name], MS2_matrix[name],how="outer")
        # print(all_keys)

        all_matrix = pd.merge(all_matrix, all_keys, how="right").replace("ID","MBR")
        MS2_matrix = pd.merge(MS2_matrix, all_keys, how="right")

        # print(all_matrix)
        for eachColumn in id_cols:
            # print(len(all_matrix[eachColumn]))
            # print(len(MS2_matrix[eachColumn]))
            if len(all_matrix[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2"))]) > 0:
                all_matrix.loc[(all_matrix[name].isin(MS2_matrix[name]) & (MS2_matrix[eachColumn] == "MS2")), [eachColumn]] = "MS2"

        # print(all_matrix)

        return all_matrix #noticed this changed

    def read_file(queue_id=None, queue_info= None, processor_info = None,
                input1=None, input2=None,input3=None, input4=None, input5=None,
                process_app = None, file_id = 1):
        """_Read data from data manager API or through local files or read directly
        in the webapp_
        Args:
            queue_id (_int_): _processing queue id_
            queue_info (_dict_): _queue info from the API_
            processor_info (_dict_): _processor info from the API_        
            input1 (_str_): _input file 1_ 
            input2 (_str_): _input file 2_
            input3 (_str_): _input file 3_
            input4 (_str_): _input file 4_
            input5 (_str_): _input file 5_
            process_app (_str_): _process app name_
        Returns:
            _dict_: _dictionary containing data all data        
        """

        """ Input files are as followws
        App       Input file
                1                         2                       3         4         5
        PD        _Proteins                 _PeptideGroups                              _InputFiles
        Fragpipe  combined_protein          combined_peptide
        DIANN     diann-output.pg_matrix    diann-output.pr_matrix  protein   peptide   filelist_diann.txt
        """

        min_unique_peptides = 1

        #getting files from data system
            # getting files from data system (webapp)
        if queue_id is not None and processor_info is None:
            # Method 1 pull data directly (used by the webapp)
            process_app = DataAnalysisQueue.objects.filter(
                pk=queue_id).first().processing_app.name
            input1= DataAnalysisQueue.objects.filter(
                pk=queue_id).first().output_file_1
            input2= DataAnalysisQueue.objects.filter(
                pk=queue_id).first().output_file_2  
            input3= DataAnalysisQueue.objects.filter(
                pk=queue_id).first().output_file_3
            input4= DataAnalysisQueue.objects.filter(
                pk=queue_id).first().output_file_4  
            input5= DataAnalysisQueue.objects.filter(
                pk=queue_id).first().output_file_5
        elif queue_info is not None and processor_info is not None:
        # Method 2 pull data from the data system (used by jupyter notebook)
            process_app = processor_info["name"]
            input1= queue_info["output_file_1"]
            input2= queue_info["output_file_2"]  
            input3= queue_info["output_file_3"]
            input4= queue_info["output_file_4"]  
            input5= queue_info["output_file_5"]
        # method 3 feed data directly (through local file paths)
        else:
            analysis_file = input1

        if "FragPipe" in process_app:     # fragpipe results
            # read data
            peptide_table = pd.read_table(input2,low_memory=False)
            protein_table = pd.read_table(input1,low_memory=False)

            #remove 
            #protein remove = 

            # pep_info_columns = ["Index", "Gene", "ProteinID", "SequenceWindow", "Start", "End", "MaxPepProb", "ReferenceIntensity"]
            pep_info_columns = ["Index", "Gene", "ProteinID", "MaxPepProb", "ReferenceIntensity"]


            prot_info_columns = ["Index", "NumberPSM", "MaxPepProb", "ReferenceIntensity"]

            # ALL
            ## Proteins abundance table
            protein_table.rename(columns={'Gene': 'Symbol'},inplace=True)
            prot_abundance = protein_table.drop(prot_info_columns,axis=1)
            prot_info_columns.append("Symbol")
            prot_other_info = protein_table.loc[:,prot_info_columns]
            ## Peptide abundance table
            peptide_table.rename(columns={'Peptide': 'Annotated Sequence'}, inplace=True)
            pep_abundance = peptide_table.drop(pep_info_columns,axis=1)
            pep_info_columns.append("Annotated Sequence")
            pep_other_info = peptide_table.loc[:,pep_info_columns]

            
            # Proteins ID table
            all_ID_cols = prot_abundance.columns
            prot_ID_MS2 = protein_table.loc[:, all_ID_cols]
            
            # Proteins ID table
            all_ID_cols = pep_abundance.columns
            pep_ID_MS2 = peptide_table.loc[:, all_ID_cols]

            # remove "Spectral Count", " MaxLFQ Intensity" or " Intensity" from names
            run_name_list = [name for name in all_ID_cols.drop("Annotated Sequence")]
            
            run_name_list = pd.DataFrame({"Run Names": run_name_list})
            run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

            # print(pep_ID_MS2.columns)

            for item in [prot_abundance,pep_abundance,pep_ID_MS2,prot_ID_MS2]:
                # Generate a new column name mapping using the function
                fileid_mapping = generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))
                item.rename(columns = fileid_mapping,inplace=True)
            

            # get ID matrix tables
            prot_ID = prot_abundance.copy()
            cols = [col for col in prot_ID.columns if col != 'Symbol']
            for col in cols:
                if prot_ID[col].dtype != 'object': # Check if not a string column
                    prot_ID[col].replace(0, np.nan, inplace=True)
                    # Replace all numerical values to ID
                    prot_ID[col] = prot_ID[col].astype(str).str.replace("\d+\.\d+", "ID", regex=True)
            pep_ID = pep_abundance.copy()
            cols = [col for col in pep_ID.columns if col != 'Annotated Sequence	']
            for col in cols:
                if pep_ID[col].dtype != 'object': # Check if not a string column
                    pep_ID[col].replace(0, np.nan, inplace=True)
                    # Replace all numerical values to ID
                    pep_ID[col] = pep_ID[col].astype(str).str.replace("\d+\.\d+", "ID", regex=True)
            #Rename protein numbers
            
            cols = [col for col in prot_ID_MS2.columns if col != 'Symbol']
            for col in cols:
                if prot_ID_MS2[col].dtype != 'object': # Check if not a string column
                    prot_ID_MS2[col].replace(0, np.nan, inplace=True)
                    # Replace all numerical values to ID
                    prot_ID_MS2[col] = prot_ID_MS2[col].astype(str).str.replace("\d+\.\d+", "MS2", regex=True)
            cols = [col for col in pep_ID_MS2.columns if col != 'Annotated Sequence	']
            for col in cols:
                if pep_ID_MS2[col].dtype != 'object': # Check if not a string column
                    pep_ID_MS2[col].replace(0, np.nan, inplace=True)
                    # Replace all numerical values to ID
                    pep_ID_MS2[col] = pep_ID_MS2[col].astype(str).str.replace("\d+\.\d+", "MS2", regex=True)

            # print(run_name_list)
            pep_ID = combine_IDs(pep_ID, pep_ID_MS2)
            prot_ID = combine_IDs(prot_ID, prot_ID_MS2)

            prot_other_info["Source_File"] = input1
            pep_other_info["Source_File"] = input2


        # elif "DIANN" in process_app:
        #     # read in DIANN output files
        #     peptide_table = pd.read_table(input2,low_memory=False)
        #     protein_table = pd.read_table(input1,low_memory=False)
        #     protein_table_MS2 = pd.read_table(input3,low_memory=False)
        #     peptide_table_MS2 = pd.read_table(input4,low_memory=False)
        #     prot_other_info = pd.DataFrame({"Protein": protein_table["Protein.Ids"], "Protein.Group": protein_table["Protein.Group"]})
        #     pep_other_info = pd.DataFrame({"Mapped Proteins": peptide_table["Protein.Group"], "Modified.Sequence": peptide_table["Modified.Sequence"]})

        #     prot_other_info["Source_File"] = "None"
        #     pep_other_info["Source_File"] = "None"

        #     # meta_table = pd.read_csv(input5, sep=' ', header=None, names=["File Name"])
        #     # filter Contaminant
        #     protein_table= protein_table[~protein_table['Protein.Group'].str.contains(
        #         "contam", na=False)]
        #     peptide_table= peptide_table[~peptide_table['Protein.Group'].str.contains(
        #         "contam", na=False)]
        #     protein_table_MS2= protein_table_MS2[~protein_table_MS2['Protein.Group'].str.contains(
        #         "contam", na=False)]
        #     peptide_table_MS2= peptide_table_MS2[~peptide_table_MS2['Protein.Group'].str.contains(
        #         "contam", na=False)]
        #     prot_other_info= prot_other_info[~prot_other_info['Protein'].str.contains(
        #         "contam", na=False)]
        #     pep_other_info= pep_other_info[~pep_other_info['Mapped Proteins'].str.contains(
        #         "contam", na=False)]
            
        #     prot_other_info.rename(columns={'Protein': 'Symbol'}, inplace=True)
        #     pep_other_info.rename(columns={'Modified.Sequence': 'Annotated Sequence'}, inplace=True)
        #     # Replace backslashes with forward slashes if data comes from Windows
        #     # meta_table['File Name'] = meta_table['File Name'].str.replace('\\', '/', regex=False)
        #     # # Apply a lambda function to extract file names without extensions
        #     # meta_table['File Name'] = meta_table['File Name'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
        #     # run_name_list = meta_table['File Name'].tolist()
        #     # run_name_list = pd.DataFrame({"Run Names": run_name_list})
        #     # run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))


        #     # Get the file names from the meta table
        #     protein_path_cols = protein_table_MS2.filter(regex='\\\\|Protein.Ids').columns

        #     ## Proteins
        #     prot_abundance = protein_table.loc[:, protein_path_cols]
        #     prot_abundance_MS2 = protein_table_MS2.loc[:, protein_path_cols]
        #     # Rename Columns to remove file path
        #     file_path_cols = protein_table.filter(regex='\\\\').columns
        #     prot_abundance.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in prot_abundance.columns]
        #     prot_abundance = prot_abundance.rename(columns={'Protein.Ids': 'Symbol'})
        #     prot_abundance["Symbol"] =  prot_abundance["Symbol"].str.replace(";.*","",regex = True)
        #     prot_abundance_MS2.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in prot_abundance_MS2.columns]
        #     prot_abundance_MS2 = prot_abundance_MS2.rename(columns={'Protein.Ids': 'Symbol'})   
        #     prot_abundance_MS2["Symbol"] =  prot_abundance_MS2["Symbol"].str.replace(";.*","",regex = True)

        #     ## Peptides
        #     peptide_path_cols = peptide_table_MS2.filter(regex='\\\\|Modified.Sequence').columns
        #     pep_abundance = peptide_table.loc[:, peptide_path_cols]
        #     pep_abundance_MS2 = peptide_table_MS2.loc[:, peptide_path_cols]
        #     # Rename Columns to remove file path
        #     file_path_cols = peptide_table.filter(regex='\\\\').columns
        #     pep_abundance.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in pep_abundance.columns]
        #     pep_abundance = pep_abundance.rename(columns={'Modified.Sequence': 'Annotated Sequence'})
        #     pep_abundance_MS2.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in pep_abundance_MS2.columns]
        #     pep_abundance_MS2 = pep_abundance_MS2.rename(columns={'Modified.Sequence': 'Annotated Sequence'})

        #     run_name_list = pd.DataFrame(data={"Run Names": [os.path.splitext(os.path.basename(x))[0] for x in file_path_cols]})
        #     run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

        #     for item in [prot_abundance,pep_abundance,prot_abundance_MS2,pep_abundance_MS2]:
        #         # Generate a new column name mapping using the function
        #         fileid_mapping = generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))
        #         item.rename(columns = fileid_mapping,inplace=True)


        #     #convert to str for IDs matrix
        #     pep_ID = pep_abundance.copy()
        #     cols = [col for col in pep_ID.columns if col != 'Symbol']
        #     for col in cols:
        #         if pep_ID[col].dtype != 'object': # Check if not a string column
        #             pep_ID[col].replace(0, np.nan, inplace=True)
        #             # Replace all numerical values to ID
        #             pep_ID[col] = pep_ID[col].astype(str).str.replace("\d+\.\d+", "ID", regex=True)
        #     pep_ID_MS2 = pep_abundance_MS2.copy()
        #     cols = [col for col in pep_ID_MS2.columns if col != 'Symbol']
        #     for col in cols:
        #         if pep_ID_MS2[col].dtype != 'object': # Check if not a string column
        #             pep_ID_MS2[col].replace(0, np.nan, inplace=True)
        #             # Replace all numerical values to ID
        #             pep_ID_MS2[col] = pep_ID_MS2[col].astype(str).str.replace("\d+\.\d+", "MS2", regex=True)
        #     prot_ID = prot_abundance.copy()
        #     cols = [col for col in prot_ID.columns if col != 'Symbol']
        #     for col in cols:
        #         if prot_ID[col].dtype != 'object': # Check if not a string column
        #             prot_ID[col].replace(0, np.nan, inplace=True)
        #             # Replace all numerical values to ID
        #             prot_ID[col] = prot_ID[col].astype(str).str.replace("\d+\.\d+", "ID", regex=True)
        #     prot_ID_MS2 = prot_abundance_MS2.copy()
        #     cols = [col for col in prot_ID_MS2.columns if col != 'Symbol']
        #     for col in cols:
        #         if prot_ID_MS2[col].dtype != 'object': # Check if not a string column
        #             prot_ID_MS2[col].replace(0, np.nan, inplace=True)
        #             # Replace all numerical values to ID
        #             prot_ID_MS2[col] = prot_ID_MS2[col].astype(str).str.replace("\d+\.\d+", "MS2", regex=True)

        #     pep_ID = combine_IDs(pep_ID, pep_ID_MS2)
        #     prot_ID = combine_IDs(prot_ID, prot_ID_MS2)       
            
        # elif "PD" in process_app:
        #     peptide_table = pd.read_table(input2,low_memory=False)
        #     protein_table = pd.read_table(input1,low_memory=False)
            
        #     # filter Contaminant
        #     protein_table= protein_table[(protein_table[
        #         "Protein FDR Confidence: Combined"] == "High") &
        #                     ((protein_table["Master"] == "IsMasterProtein") | 
        #                      (protein_table["Master"] == "Master")) & 
        #                     (protein_table["Contaminant"] == False)]

        #     protein_table.rename(
        #         columns={'# Peptides': 'number of peptides'}, inplace=True)
        #     protein_table=protein_table.query(
        #         "`number of peptides` >= @min_unique_peptides")
        #     peptide_table= peptide_table[(peptide_table[
        #         'Contaminant'] == False) & (peptide_table["Confidence"]== "High")]

        #     meta_table = pd.read_table(input5,low_memory=False)
        #     #filter rows in meta table on File ID column if it is NaN
        #     meta_table = meta_table[meta_table['File ID'].notna()]

        #     # Replace single backslashes with forward slashes in the 'file_paths' column
        #     meta_table['File Name'] = meta_table['File Name'].str.replace('\\', '/', regex=False)
        #     # Apply a lambda function to extract file names without extensions
        #     meta_table['file_names'] = meta_table['File Name'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
        #     file_path_name_dict = dict(zip(meta_table['File ID'], meta_table['file_names']))
        #     run_name_list = pd.DataFrame({"Run Names": file_path_name_dict.values()})
        #     run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))
            
        #     #format the read in table into three different tables: abundance, id and other_info
        #     prot_abundance = protein_table.filter(regex='Abundance:|Symbol')
        #     prot_ID = protein_table.filter(regex='Found in Sample:|Symbol')
        #     prot_other_info = protein_table.loc[:, ~protein_table.columns.str.contains('Found in Sample:|Abundance:')]
            

        #     pep_abundance = peptide_table.filter(regex='Abundance:|Annotated Sequence')
        #     pep_ID = peptide_table.filter(regex='Found in Sample:|Annotated Sequence')
        #     pep_other_info = peptide_table.loc[:, ~peptide_table.columns.str.contains('Found in Sample:|Abundance:')]

        #     prot_other_info["Source_File"] = input1
        #     pep_other_info["Source_File"] = input2

        #     #change column names to file/run names to our fileID

        #     new_dict = {"Abundance: " + key + ":": value for key, value in file_path_name_dict.items()}
        #     for item in [prot_abundance,pep_abundance]:
        #         # Generate a new column name mapping using the function
        #         column_name_mapping = generate_column_from_name_mapping(item.columns, new_dict)
        #         #TODO solving  A value is trying to be set on a copy of a slice from a DataFrame
                
        #         item.rename(columns = column_name_mapping, inplace = True)
        #         #use generate_column_to_name_mapping function because we don't want partial matches as in QC_HeLa.raw and QC_HeLa_20230727235101.raw
        #         fileid_mapping = generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))
        #         item.rename(columns = fileid_mapping,inplace=True)
            

        #     new_dict = {"Found in Sample: " + key + ":": value for key, value in file_path_name_dict.items()}
        #     for item in [pep_ID,prot_ID]:
        #         # Generate a new column name mapping using the function
        #         column_name_mapping = generate_column_from_name_mapping(item.columns, new_dict)

        #         item.rename(columns = column_name_mapping, inplace = True)
        #         #use generate_column_to_name_mapping function because we don't want partial matches
        #         fileid_mapping = generate_column_to_name_mapping(item.columns, dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"])))

        #         item.rename(columns = fileid_mapping,inplace=True)

        #     # replace "High" to MS2 "Peak Found" to MBR, the rest become np.NaN
        #     replacements = {'High': 'MS2', 'Peak Found': 'MBR', "Medium": np.NaN, "Low": np.NaN, "Not Found": np.NaN}
        #     for column in run_name_list["Run Identifier"]:
        #         if column in pep_ID.columns:
        #             pep_ID[column] = pep_ID[column].replace(to_replace=replacements)
        
        #         if column in prot_ID.columns:
        #             prot_ID[column] = prot_ID[column].replace(to_replace=replacements)
        
        # get ID summary by parsing ID Matrix
        protein_ID_summary = sumIDs(prot_ID)
        peptide_ID_summary = sumIDs(pep_ID)
        

        #sets the processing app in run_name_list
        run_name_list["Processing App"] = process_app
        run_name_list["Analysis Name"] = analysis_file


        return {'run_metadata': run_name_list,
                'protein_other_info': prot_other_info,
                'peptide_other_info': pep_other_info,
                'protein_abundance': prot_abundance,
                'protein_ID_matrix': prot_ID,
                'protein_ID_Summary': protein_ID_summary,
                'peptide_abundance': pep_abundance,
                'peptide_ID_matrix': pep_ID,
                'peptide_ID_Summary': peptide_ID_summary,

                }  

    def read_files(queue_ids = None, queue_info = None, processor_info = None, grouped_input_files = []):
        '''
        Creates a list of data objects
        
        Input
        grouped_input_files
        [
            #File 0
        {input1:
        input2:
        input3:   
        input4:
        input5:
        process_app:
        },
            #File 1
        {input1:
        input2:
        input3:   
        input4:
        input5:
        process_app:
        },
        ...
        ]
        '''


        data_objects = []

        i = 0
        for eachGroup in grouped_input_files:
            if queue_ids is not None:
                pass
            else:
                process_app = eachGroup["process_app"]
                input1= eachGroup["input1"]
                input2= eachGroup["input2"]  
                input3= eachGroup["input3"]
                input4= eachGroup["input4"]  
                input5= eachGroup["input5"]

            current_data_object = read_file(input1=input1,input2=input2,
                                            input3=input3,input4 = input4,
                                            input5=input5, process_app=process_app,file_id = i)
            data_objects.append(current_data_object)
            
            i = i + 1

        
        return data_objects


    def outer_join_data_objects(data_objects):
        '''
        Takes in a list of data objects as given by read_files and converts them to a single data object as given by read_files,
        protein info continues to show what was found on each original file, and so forth.
        '''

        first_file = True
        for eachDataObject in data_objects:
            print("***")
            if first_file:
                first_file = False
                final_data_object = eachDataObject
            else:
                final_data_object['run_metadata'] = pd.concat([final_data_object['run_metadata'],eachDataObject['run_metadata']]).reset_index(drop=True)
                final_data_object['protein_other_info'] = pd.concat([final_data_object['protein_other_info'],eachDataObject['protein_other_info']]).reset_index(drop=True)
                final_data_object['peptide_other_info'] = pd.concat([final_data_object['peptide_other_info'],eachDataObject['peptide_other_info']]).reset_index(drop=True)
                final_data_object['protein_ID_Summary'] = pd.concat([final_data_object['protein_ID_Summary'],eachDataObject['protein_ID_Summary']]).reset_index(drop=True)
                final_data_object['peptide_ID_Summary'] = pd.concat([final_data_object['peptide_ID_Summary'],eachDataObject['peptide_ID_Summary']]).reset_index(drop=True)
                duplicates_found = False
                
                #loop through to see if there are any duplicate files
                for eachCol in final_data_object['protein_abundance'].loc[:, final_data_object['protein_abundance'].columns!='Symbol'].columns:
                    if eachCol in eachDataObject['protein_abundance'].columns:
                        duplicates_found = True
                    else:
                        pass
                for eachCol in final_data_object['protein_ID_matrix'].loc[:, final_data_object['protein_ID_matrix'].columns!='Symbol'].columns:
                    if eachCol in eachDataObject['protein_ID_matrix'].columns:
                        duplicates_found = True
                    else:
                        pass
                for eachCol in final_data_object['peptide_abundance'].loc[:, final_data_object['peptide_abundance'].columns!='Annotated Sequence'].columns:
                    if eachCol in eachDataObject['peptide_abundance'].columns:
                        duplicates_found = True
                    else:
                        pass
                for eachCol in final_data_object['peptide_ID_matrix'].loc[:, final_data_object['peptide_ID_matrix'].columns!='Annotated Sequence'].columns:
                    if eachCol in eachDataObject['peptide_ID_matrix'].columns:
                        duplicates_found = True
                    else:
                        pass     
                if duplicates_found:
                    print("Error: files analyzed twice present!!!")
                    quit()
                    print("@#afio2q3")
                else:
                    #merge keeping all proteins
                    # print("!!!!")
                    final_data_object['protein_abundance'] = pd.merge(final_data_object['protein_abundance'],eachDataObject['protein_abundance'],how="outer")
                    final_data_object['protein_ID_matrix'] = pd.merge(final_data_object['protein_ID_matrix'],eachDataObject['protein_ID_matrix'],how="outer")
                    final_data_object['peptide_abundance'] = pd.merge(final_data_object['peptide_abundance'],eachDataObject['peptide_abundance'],how="outer")
                    final_data_object['peptide_ID_matrix'] = pd.merge(final_data_object['peptide_ID_matrix'],eachDataObject['peptide_ID_matrix'],how="outer")
                    
        return final_data_object

    def calculate_missing_values_MS2(data_object,
                                missing_value_thresh=33,
                                is_protein=True,
                                ignore_nan=False):
        """_Filter out proteins/peptides with missing values rate above the
        threshold_

        Args:
            data_object (_panada_): _dataframe contain data for one experimental
            condition_
            missing_value_thresh (int, optional): _description_. Defaults to 33.
            analysis_program (str, optional): _description_.
            ignore_nan: if filter intensity again with Nan threadshold, this 
            helps with the calcualting stdev step.

        Returns:
            _data_object_: _dictionary containing data for one experimental
            'abundances':        Symbol  3_TrypsinLysConly_3A4_channel2 3_TrypsinLysConly_3BC_channel1
    0     A0A096LP49                            0.00                                        10
    1     A0A0B4J2D5                        89850.26                                      3311
    2         A0AVT1                        83055.87                                    312312
        """
        if is_protein:
            name = "Symbol"
            matrix_name = "protein_ID_matrix"
            other_info_name = "protein_other_info"
            abundance_name = "protein_abundance"
            
        else:
            name = "Annotated Sequence"
            matrix_name = "peptide_ID_matrix"
            other_info_name = "peptide_other_info"
            abundance_name = "peptide_abundance"

        #initializes number of missing values to zero
        protein_columns = data_object[matrix_name].assign(missingValues=0)

        i = 0
        # found all the proteins/peptides with missing values rate below
        # the threshold, pep_columns contains the remaining protein/peptide
        # in a pandas dataframe with $names as its column name
        print(data_object[matrix_name].columns)
        for each_column in data_object[matrix_name].loc[:, ~data_object[matrix_name].columns.str.contains(name)].columns:
            # replace "nan" to np.nan
            protein_columns = protein_columns.replace({"nan": np.nan,"NA":np.nan}) 
            protein_columns.loc[protein_columns[each_column] != "MS2", #ID/MBR are still missing values if you are only considering MS2
                                "missingValues"] += 1

            i += 1

        protein_columns = protein_columns.assign(missingValuesRate=(
            protein_columns["missingValues"] / i) * 100)
        
        returnMatrix = pd.DataFrame({"Missing Values Rate": protein_columns["missingValuesRate"], name: protein_columns[name]})
        
        
        
        return returnMatrix

    def filter_by_missing_values(data_object,
                                missing_value_thresh=33,
                                is_protein=True,
                                ignore_nan=False):
        """_Filter out proteins/peptides with missing values rate above the
        threshold_

        Args:
            data_object (_panada_): _dataframe contain data for one experimental
            condition_
            missing_value_thresh (int, optional): _description_. Defaults to 33.
            analysis_program (str, optional): _description_.
            ignore_nan: if filter intensity again with Nan threadshold, this 
            helps with the calcualting stdev step.

        Returns:
            _data_object_: _dictionary containing data for one experimental
            'abundances':        Symbol  3_TrypsinLysConly_3A4_channel2 3_TrypsinLysConly_3BC_channel1
    0     A0A096LP49                            0.00                                        10
    1     A0A0B4J2D5                        89850.26                                      3311
    2         A0AVT1                        83055.87                                    312312
        """
        if is_protein:
            name = "Symbol"
            matrix_name = "protein_ID_matrix"
            other_info_name = "protein_other_info"
            abundance_name = "protein_abundance"
            
        else:
            name = "Annotated Sequence"
            matrix_name = "peptide_ID_matrix"
            other_info_name = "peptide_other_info"
            abundance_name = "peptide_abundance"

        #initializes number of missing values to zero
        protein_columns = data_object[matrix_name].assign(missingValues=0)

        i = 0
        # found all the proteins/peptides with missing values rate below
        # the threshold, pep_columns contains the remaining protein/peptide
        # in a pandas dataframe with $names as its column name
        for each_column in data_object[matrix_name].loc[
                :, ~data_object[matrix_name].columns.str.contains(
                    name)].columns:
            # replace "nan" to np.nan
            protein_columns = protein_columns.replace({"nan": np.nan}) 
            #find missing values and increment those rows (a row is a protein/peptide) total number of missing values

            protein_columns.loc[(protein_columns[each_column] != "MS2")
                                &(protein_columns[each_column] != "MBR")
                                &(protein_columns[each_column] != "ID"), #this is more robust than using nan's in case something fails to convert
                                "missingValues"] += 1

            i += 1

        protein_columns = protein_columns.assign(missingValuesRate=(
            protein_columns["missingValues"] / i) * 100)
        
        protein_columns = protein_columns.query(
            "missingValuesRate < @missing_value_thresh")
        
        protein_columns = protein_columns.loc[:,
                                    protein_columns.columns.str.contains(name)]

        # filter the data_object with the remaining proteins/peptides names
        data_object[abundance_name] = protein_columns.merge(
            data_object[abundance_name])
        data_object[matrix_name] = protein_columns.merge(
            data_object[matrix_name])
        data_object[other_info_name] = protein_columns.merge(
            data_object[other_info_name])
        # In case there is mismatch between ID table and abundance table,
        # mannually remove the row with all NaN values
        # keep rows in data_object[abundance_name] where at least two values are 
        # not NaN(do this to all rows except the first row), otherwise can't
        # calculate the stdev
        if ignore_nan:
            data_object[abundance_name] = data_object[abundance_name].dropna(
                thresh=2, subset=data_object[abundance_name].columns[1:])
            # This will cause the veen diagram to be different from R program
        
        return data_object

    def filter_by_missing_values_MS2(data_object,
                                missing_value_thresh=33,
                                is_protein=True,
                                ignore_nan=False):
        """_Filter out proteins/peptides with missing values rate above the
        threshold_

        Args:
            data_object (_panada_): _dataframe contain data for one experimental
            condition_
            missing_value_thresh (int, optional): _description_. Defaults to 33.
            analysis_program (str, optional): _description_.
            ignore_nan: if filter intensity again with Nan threadshold, this 
            helps with the calcualting stdev step.

        Returns:
            _data_object_: _dictionary containing data for one experimental
            'abundances':        Symbol  3_TrypsinLysConly_3A4_channel2
    0     A0A096LP49                            0.00
    1     A0A0B4J2D5                        89850.26
    2         A0AVT1                        83055.87
        """
        if is_protein:
            name = "Symbol"
            matrix_name = "protein_ID_matrix"
            other_info_name = "protein_other_info"
            abundance_name = "protein_abundance"
            
        else:
            name = "Annotated Sequence"
            matrix_name = "peptide_ID_matrix"
            other_info_name = "peptide_other_info"
            abundance_name = "peptide_abundance"

        protein_columns = data_object[matrix_name].assign(missingValues=0)

        i = 0
        # found all the proteins/peptides with missing values rate below
        # the threshold, pep_columns contains the remaining protein/peptide
        # in a pandas dataframe with $names as its column name
        for each_column in data_object[matrix_name].loc[:, ~data_object[matrix_name].columns.str.contains(name)].columns:
            # replace "nan" to np.nan
            protein_columns = protein_columns.replace({"nan": np.nan}) 
            protein_columns.loc[protein_columns[each_column] != "MS2", #ID/MBR are still missing values if you are only considering MS2
                                "missingValues"] += 1

            i += 1

        protein_columns = protein_columns.assign(missingValuesRate=(
            protein_columns["missingValues"] / i) * 100)
        
        protein_columns = protein_columns.query(
            "missingValuesRate < @missing_value_thresh")
        
        protein_columns = protein_columns.loc[:,
                                    protein_columns.columns.str.contains(name)]

        # filter the data_object with the remaining proteins/peptides names
        data_object[abundance_name] = protein_columns.merge(
            data_object[abundance_name])
        data_object[matrix_name] = protein_columns.merge(
            data_object[matrix_name])
        data_object[other_info_name] = protein_columns.merge(
            data_object[other_info_name])
        # In case there is mismatch between ID table and abundance table,
        # mannually remove the row with all NaN values
        # keep rows in data_object[abundance_name] where at least two values are 
        # not NaN(do this to all rows except the first row), otherwise can't
        # calculate the stdev
        if ignore_nan:
            data_object[abundance_name] = data_object[abundance_name].dropna(
                thresh=2, subset=data_object[abundance_name].columns[1:])
            # This will cause the veen diagram to be different from R program
        
        return data_object


    def NormalizeToMedian(abundance_data, apply_log2=False):
        """_Normalizes each column by multiplying each value in that column with
        the median of all values in abundances (all experiments) and then dividing
        by the median of that column (experiment)._
        we find applying log2 transform first gives more robust results for PCA etc.
        See https://pubs.acs.org/doi/10.1021/acsomega.0c02564
        Args:
            abundance_data (_pd_): _description_
            apply_log2 (_bool_,): _apply log2 to all result_.
        Returns:
            _type_: _description_
            format:
            'abundances':        Symbol  3_TrypsinLysConly_3A4_channel2
            A0A096LP49                    0.000000e+00
        """
        # all the columns/sample list
        columns = [col for col in abundance_data.select_dtypes(include=[
                np.number])]
        data_matrix = abundance_data[columns].values
        # replace 0 with nan
        data_matrix[data_matrix == 0] = np.nan
        medianOfAll = np.nanmedian(data_matrix)
        
        #normalize all median, all median/current run all protein median
        # apply log2 to all the values if apply_log2 is True
        if apply_log2:    
            for each_column in columns:
                abundance_data[each_column] = (
                    np.log2(medianOfAll) * np.log2(abundance_data[each_column]) /
                    np.log2(np.nanmedian(abundance_data[
                        each_column].replace(0, np.nan))))
        else:
            for each_column in columns:
                abundance_data[each_column] = (
                    medianOfAll * abundance_data[each_column] /
                    np.nanmedian(abundance_data[
                        each_column].replace(0, np.nan)))
        #TODO divide by zero error encountered in log2, temporarily set to 0
        abundance_data = abundance_data.replace([np.inf, -np.inf], 0)

        return abundance_data

    def calculate_cvs(abundance_data):
        """_Calculate mean, stdev, cv for withn each protein/peptide abundance_

        Args:
            data_object (_type_): _full data frame_

        Returns:
            _type_: _df with Symbol mean, stdev, cv for each protein/peptide_
        """
        if 'Symbol' in abundance_data.columns:
            name = "Symbol"
        if 'Annotated Sequence' in abundance_data.columns:
            name = "Annotated Sequence"
        abundance_data = abundance_data.assign(
            intensity=abundance_data.loc[:, ~abundance_data.columns.str.contains(
                name)].mean(axis=1, skipna=True),
            stdev=abundance_data.loc[:, ~abundance_data.columns.str.contains(
                name)].std(axis=1, skipna=True),
            CV=abundance_data.loc[:, ~abundance_data.columns.str.contains(name)].std(
                axis=1, skipna=True) / abundance_data.loc[
                :, ~abundance_data.columns.str.contains(name)].mean(
                axis=1, skipna=True) * 100)

        abundance_data = abundance_data.loc[:, [
                name, "intensity", "stdev", "CV"]]
        
        return abundance_data


    def t_test_from_summary_stats(m1, m2, n1, n2, s1, s2, equal_var=False):
        """_Calculate T-test from summary using ttest_ind_from_stats from
        scipy.stats package_

        Args:
            m1 (_type_): _mean list of sample 1_
            m2 (_type_): mean list of sample 2_
            n1 (_type_): sample size list of sample 1_
            n2 (_type_): sample size list of sample 2_
            s1 (_type_): standard deviation list of sample 1_
            s2 (_type_): standard deviation list of sample 2_
            equal_var (_type_, optional): False would perform Welch's
            t-test, while set it to True would perform Student's t-test. Defaults
            to False.

        Returns:
            _type_: _list of P values_
        """

        p_values = []
        for i in range(len(m1)):
            _, benjamini = ttest_ind_from_stats(
                m1[i], s1[i], n1[i], m2[i], s2[i], n2[i], equal_var=equal_var)
            p_values.append(benjamini)

        return p_values

    def impute_knn(abundance_data, k=5):
        """_inpute missing value from neighbor values_

        Args:
            abundance_data (_type_): _description_
            k (int, optional): _number of neighbors used_. Defaults to 5.
        Returns:
            _type_: _description_
            TODO: this knn imputer produces slightly different results (about 4%)
            from the one in R. Need to figure out why
        """
        name = abundance_data.columns[0]

        names = abundance_data[name]
        # x = abundance_data.select_dtypes(include=['float64', 'int64'])
        # imputer = KNNImputer(n_neighbors=k)
        # x_imputed = pd.DataFrame(imputer.fit_transform(x), columns=x.columns)


        x = abundance_data.select_dtypes(include=['float', 'int'])
        imputer = KNNImputer(n_neighbors=k)
        x_imputed = imputer.fit_transform(x)
        x_imputed = pd.DataFrame(x_imputed, columns=x.columns)



        abundance_data.loc[:, x.columns] = x_imputed.values
        abundance_data[name] = names
        return abundance_data


    def CalculatePCA(abundance_object, infotib,log2T = False):
        """_inpute PCA transformed and variance explained by each principal
        component_
        """
        name = abundance_object.columns[0]
        x = abundance_object
        
        sampleNames = x.columns[~x.columns.str.contains(
            name)].to_frame(index=False)

        if log2T: #apply log2 transformation
            x = np.log2(x.loc[:, ~x.columns.str.contains(name)].T.values)
        else:
            x = x.loc[:, ~x.columns.str.contains(name)].T.values
        # filter out columns with all zeros
        is_finite_col = np.isfinite(np.sum(x, axis=0))
        x_filtered = x[:, is_finite_col]

        
        # Instantiate PCA    
        pca = PCA()
        #
        # Determine transformed features
        #
        x_pca = pca.fit_transform(x_filtered)
        #
        # Determine explained variance using explained_variance_ration_ attribute
        #
        exp_var_pca = pca.explained_variance_ratio_
        #
        # Cumulative sum of eigenvalues; This will be used to create step plot
        # for visualizing the variance explained by each principal component.
        #
        cum_sum_eigenvalues = np.cumsum(exp_var_pca)
        #
        # convert numpy array to pandas dataframe for plotting
        
        pca_panda = pd.DataFrame(x_pca, columns=[
            'PC' + str(i+1) for i in range(x_pca.shape[1])])
        # add sample names to the dataframe
        pca_panda = pd.concat(
            [infotib, pca_panda], axis=1, join='inner')
        
        return pca_panda, exp_var_pca


    def filter_by_name(data_dict, runname_list):
        """_Filter the data_dict based on runname_list, only keep the columns
        of the data_dict that are in the runname_list_
        Args:

        Returns:
            _type_: _description_
        """

        # make dict for each runname, no Symbol/sequence
        nameDict = dict(zip(data_dict["run_metadata"]["Run Names"],data_dict["run_metadata"]["Run Identifier"]))
        
        identifier_list = []
        
        identifier_list_plus = []
        if "Annotated Sequence" in runname_list:
            runname_list.remove("Annotated Sequence")
        if "Symbol" in runname_list:
            runname_list.remove("Symbol")
        for eachName in runname_list:
            identifier_list.append(nameDict[eachName])

        for eachName in runname_list:
            identifier_list_plus.append(nameDict[eachName])


        filtered_data = {}
    # filtered_data["meta"] = data_dict["meta"]
        runname_list.extend(["Annotated Sequence","Symbol"])
        identifier_list_plus.extend(["Annotated Sequence","Symbol"])

        #filtered_data["run_metadata"] = [item for item in data_dict[
        #   "run_metadata"] if item in runname_list]
        
        filtered_data["run_metadata"] = data_dict["run_metadata"][
            data_dict["run_metadata"]["Run Names"].isin(
                runname_list)]  
        filtered_data["protein_abundance"] = data_dict["protein_abundance"][[
            col for col in data_dict["protein_abundance"].columns if any(
                word == col for word in identifier_list_plus)]]
        filtered_data["peptide_abundance"] = data_dict["peptide_abundance"][[
            col for col in data_dict["peptide_abundance"].columns if any(
                word == col for word in identifier_list_plus)]]
        filtered_data["protein_other_info"] = data_dict["protein_other_info"][[
            col for col in data_dict["protein_other_info"].columns if any(
                word == col for word in identifier_list_plus)]]
        filtered_data["peptide_other_info"] = data_dict["peptide_other_info"][[
            col for col in data_dict["peptide_other_info"].columns if any(
                word == col for word in identifier_list_plus)]]
        filtered_data["protein_ID_matrix"] = data_dict["protein_ID_matrix"][[
            col for col in data_dict["protein_ID_matrix"].columns if any(
                word == col for word in identifier_list_plus)]]
        filtered_data["peptide_ID_matrix"] = data_dict["peptide_ID_matrix"][[
            col for col in data_dict["peptide_ID_matrix"].columns if any(
                word == col for word in identifier_list_plus)]]
        filtered_data["protein_ID_Summary"] = data_dict["protein_ID_Summary"][
            data_dict["protein_ID_Summary"]["names"].isin(
                identifier_list)]
        filtered_data["peptide_ID_Summary"] = data_dict["peptide_ID_Summary"][
            data_dict["peptide_ID_Summary"]["names"].isin(
                identifier_list)]
        return filtered_data
