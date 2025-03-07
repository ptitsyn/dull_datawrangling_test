import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#don't forget to install openpyxl to read excel in python
dspecs=pd.read_excel('Technical Test - Data Wrangling - v20240923.xlsx', 'Data Specification')
dspecs=dspecs.iloc[:12]

df_clinical=pd.read_excel('Technical Test - Data Wrangling - v20240923.xlsx', 'Patient_clinical_data')
df_tsample=pd.read_excel('Technical Test - Data Wrangling - v20240923.xlsx', 'Tissue Sample Metadata')
df_sprot=pd.read_excel('Technical Test - Data Wrangling - v20240923.xlsx', 'Serum Protein data')
df_rpkm=pd.read_excel('Technical Test - Data Wrangling - v20240923.xlsx', 'RNA-seq (RPKM)')

df_example=pd.read_excel('Technical Test - Data Wrangling - v20240923.xlsx', 'Example report')

#transpose the rna-seq data
rpkm_cols=df_rpkm.pop('GeneID')
df_rpkm_t=df_rpkm.transpose()
df_rpkm_t.columns=rpkm_cols

#rename columns to conform the report template
df_tsample.columns=['Patient_ID','Sample_ID','Sample_General_Pathology','Material_Type','RIN','Total Reads(millions)']
df_tsample=df_tsample.drop(['RIN','Total Reads(millions)'],axis=1)

#join rna-seq results to rna-seq metadata co connect results to samples to patients
df_rpkm_sample=df_tsample.set_index('Sample_ID').join(df_rpkm_t,how='left')

#join proteomics and rna-seq tables on patient ID index
df_sprot_rpkm=df_rpkm_sample.set_index(df_rpkm_sample.columns[0]).join(df_sprot.set_index('Patient'), how='left')

#concatenate study and patient IDs element-wise in the clinical table to form unique ID as desired
sid=df_clinical[df_clinical.columns[0]].to_list()
pid=df_clinical[df_clinical.columns[1]].to_list()
uid=[x +'_'+ str(y) for x, y in zip(sid, pid)]
df_clinical['Unique_Patient_ID']=uid

#join clinical table to the protein-rnaseq table on patient IDs
df_w=df_sprot_rpkm.join(df_clinical.set_index(df_clinical.columns[1]),how='left')

block_df_list=[]
#transform the analysis cols to conform the report template
for patient in pid:
    p_table=df_w.loc[patient].fillna('NA')
    df_an=p_table[['ICAM1','IL6','IL6R','VCAM1','SELE','Serum IL-6 (g/L)','Serum IL-6 Receptor (mg/L)']]
    #form the column with gene names and
    gene_ids=df_an.columns.tolist()
    #form the measurement column
    a_measure=[]; a_geneid=[]; a_type=[];a_x=[];
    for name in gene_ids:
        a_x.append(patient)
        if('Serum' in name):
            m = name.split(' ')
            a_measure.append(m[-1])
            a_geneid.append(name[:(name.find('(')-1)])
            a_type.append('SERUM')
        else:
            a_measure.append('RPKM')
            a_geneid.append(name)
            a_type.append('RNA')

    #form geneid column by repeating the template
    a_gene=[a_geneid]*df_an.shape[0]
    a_gene = np.concatenate(a_gene).tolist()
    #form measure column by repeating the template
    a_metr=[a_measure]*df_an.shape[0]
    a_metr = np.concatenate(a_metr).tolist()
    #form meterial type column
    a_mt=[a_type]*df_an.shape[0]
    a_mt = np.concatenate(a_mt).tolist()
    #form index
    a_ind=[a_x]*df_an.shape[0]
    a_ind = np.concatenate(a_ind).tolist()
    #form the measured values into a single column
    a_values=df_an.values.flatten()
    a_block=pd.DataFrame(data=a_values,columns=['Result'],index=a_ind)
    a_block['Gene_Symbol']=a_gene
    a_block['Result_Units']=a_metr
    a_block['Material_Type']=a_mt

    df_cl=p_table[['Sample_General_Pathology', 'Sample', 'Study_ID', 'Sex', 'Age', 'Unique_Patient_ID']]
    df_cl=df_cl.rename(columns={'Sample' : 'Sample_ID'})

    block=a_block.join(df_cl,how='left')
    block_df_list.append(block)
#assemble the table from blocks
semifinal=pd.concat(block_df_list).rename_axis('Patient_ID')
semifinal=semifinal.reset_index()
#add status column
stat = ['NA'] * semifinal.shape[0]
semifinal['Status']=stat
#make sure the columns are in caps as specified
semifinal['Study_ID']=semifinal['Study_ID'].str.upper()
semifinal['Patient_ID']=semifinal['Study_ID'].str.upper()
semifinal['Sample_ID']=semifinal['Study_ID'].str.upper()
semifinal['Gene_Symbol']=semifinal['Study_ID'].str.upper()
semifinal['Unique_Patient_ID']=semifinal['Unique_Patient_ID'].str.upper()

final=semifinal[df_example.columns]
final.to_csv('curated_data.csv',index=False)
#print(final.columns)
#print(final)


exit(0)