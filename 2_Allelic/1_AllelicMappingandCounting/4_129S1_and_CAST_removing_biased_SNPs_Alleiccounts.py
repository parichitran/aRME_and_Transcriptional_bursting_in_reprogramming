#659 to be processed
import os
import pandas as pd
import sys


#mysub function to get count of only take count
def mysub(r):
    lst = [float(a) for a in r.split(',')]
    try:
        return sum(lst) / len(lst) 
    except ZeroDivisionError:
        return 0


#mysub_len funtion count ref and alt snp position 
def mysub_len(r):
    lst = [float(a) for a in r.split(',') ]
    try:
        return len(lst) # len(lst) 
    except ZeroDivisionError:
        return 0

file_path="/aRME_and_Transcriptional_bursting_in_reprogramming/2_Allelic/1_AllelicMappingandCounting"


#empty dataframe
df_129S1_SvImJ=pd.DataFrame(columns = ['ID'])
df_CAST=pd.DataFrame(columns = ['ID'])
df_both=pd.DataFrame(columns = ['ID'])

for filename in os.listdir(file_path):
    if filename.endswith("processed_gene.v3.1.txt"):#getting all gene file
       x=filename.split("_")
       name=str(x[0])
       data=pd.read_csv(file_path+"/"+filename, sep="\t")
       #data=data[data["contig"]=="X"]
       data=data[['ID', 'contig', 'aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']]
       #data=explode(data, ['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1'], fill_value='')
       df=data
       #https://stackoverflow.com/questions/63488123/ungroup-pandas-dataframe-column-values-separated-by-comma/63488203#63488203
       df[['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']]=df[['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']].apply(lambda x : x.str.split(','))
       df = df.join(pd.concat([df.pop(x).explode() for x in ['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']],axis=1))
       df["count"]= df['aggr_ref_rep1'].astype(int) + df["aggr_alt_rep1"].astype(int)
       df=df.loc[df["count"] >= 3]
       if df["count"].sum() >100:
           df.drop('count', axis=1, inplace=True)
           df =df.groupby(['ID','contig']).agg(', '.join).reset_index()
           df["count_snp"]=df["aggr_ref_rep1"].apply(mysub_len)
           df=df.loc[df["count_snp"] >=2]
           df["129S1_SvImJ"]=df["aggr_ref_rep1"].apply(mysub)
           df["CAST"]=df["aggr_alt_rep1"].apply(mysub)
           df=df.round(0)
           df["129S1_SvImJ_ratio"]=df["129S1_SvImJ"]/(df["129S1_SvImJ"]+df["CAST"])
           df["CAST_ratio"]=df["CAST"]/(df["129S1_SvImJ"]+df["CAST"])
           df_mo1=df.loc[df["129S1_SvImJ_ratio"]>=0.85] #monoallelic 129S1_SvImJ
           df_mo2=df.loc[df["CAST_ratio"]>=0.85] #monoallelic CAST
           df_mono=pd.concat([df_mo1, df_mo2])
           df_bi=df[~df["ID"].isin(list(df_mono["ID"]))]# biallelic genes cat
           df_bi[['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']]=df_bi[['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']].apply(lambda x : x.str.split(','))
           df_bi = df_bi.join(pd.concat([df_bi.pop(x).explode() for x in ['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']],axis=1))
           df_bi.drop(["129S1_SvImJ_ratio", "CAST_ratio","129S1_SvImJ", "CAST"], axis=1, inplace=True) #drop unwanted columns
           #df_bi=df_bi.astype(int)
           df_bi["129S1_SvImJ_snp_ratio"]=(df_bi["aggr_ref_rep1"].astype(int)/(df_bi["aggr_ref_rep1"].astype(int)+df_bi["aggr_alt_rep1"].astype(int)))*100 #percentage cal 129S1_SvImJ SNP
           df_bi["CAST_snp_ratio"]=(df_bi["aggr_alt_rep1"].astype(int)/(df_bi["aggr_ref_rep1"].astype(int)+df_bi["aggr_alt_rep1"].astype(int)))*100  #percentage cal CAST SNP
           # conflict SNPs identification in biallelic cat
           df_bi_conflict_snp_129S=df_bi[df_bi["129S1_SvImJ_snp_ratio"]<1]# SNP removing below 1% expressed allele 
           df_bi_conflict_snp_CAST=df_bi[df_bi["CAST_snp_ratio"]<1]# SNP removing below 1% expressed allele
           conflict_snp=list(df_bi_conflict_snp_129S["aggr_pos"])+list(df_bi_conflict_snp_CAST["aggr_pos"])#conflict SNPs position
           #remove above SNPs from the main dataframe 
           df[['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']]=df[['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']].apply(lambda x : x.str.split(','))
           df = df.join(pd.concat([df.pop(x).explode() for x in ['aggr_pos','aggr_ref_rep1', 'aggr_alt_rep1']],axis=1))
           df.drop(["129S1_SvImJ_ratio", "CAST_ratio","129S1_SvImJ", "CAST"], axis=1, inplace=True)
           df=df[~df["aggr_pos"].isin(conflict_snp)]
           #df.to_csv("Test1", sep="\t")
           df=df.astype(str)
           df =df.groupby(['ID','contig']).agg(','.join).reset_index()
           df["count_snp"]=df["aggr_ref_rep1"].apply(mysub_len)
           df=df[df["count_snp"] >=2]
           #df.to_csv("Test2", sep="\t")
           df["129S1_SvImJ"]=df["aggr_ref_rep1"].apply(mysub)
           df["CAST"]=df["aggr_alt_rep1"].apply(mysub)
           df=df.round(0)
           data=df
           data=data[["ID", "CAST", "129S1_SvImJ"]]
           data_129S1_SvImJ=data[["ID", "129S1_SvImJ"]]
           data_CAST=data[["ID","CAST"]]
           data.columns=data.columns.map(lambda x : x+"_"+name if x !='ID' else x)
           data_129S1_SvImJ.columns=data_129S1_SvImJ.columns.map(lambda x : x+"_"+name if x !='ID' else x)
           data_CAST.columns=data_CAST.columns.map(lambda x : x+"_"+name if x !='ID' else x)
           df_both = df_both.merge(data, on=['ID'], how='outer')#merge male exp
           df_both = df_both.fillna(0)
           df_both.to_csv("129S1_SvImJ_vs_CAST.txt", sep="\t", index=False)
           #129S1_SvImJ
           df_129S1_SvImJ = df_129S1_SvImJ.merge(data_129S1_SvImJ, on=['ID'], how='outer')#merge male exp
           df_129S1_SvImJ = df_129S1_SvImJ.fillna(0)
           df_129S1_SvImJ.columns = df_129S1_SvImJ.columns.str.replace(r'129S1_SvImJ_', '')
           df_129S1_SvImJ.to_csv("129S1_SvImJ.txt", sep="\t", index=False)
           #CAST
           df_CAST = df_CAST.merge(data_CAST, on=['ID'], how='outer')#merge male exp
           df_CAST = df_CAST.fillna(0)
           df_CAST.columns = df_CAST.columns.str.replace(r'CAST_', '')
           df_CAST.to_csv("CAST.txt", sep="\t", index=False)
       else:
           print(name)

