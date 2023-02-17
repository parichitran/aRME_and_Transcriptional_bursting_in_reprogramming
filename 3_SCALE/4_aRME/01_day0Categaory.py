import pandas as pd
import numpy as np

# Read the OTNine and CAST ratio file 
      
OTNine=pd.read_csv("AlleleAratio.txt", sep="\t")
CAST=pd.read_csv("AlleleBratio.txt", sep="\t")
OTNine=OTNine.rename_axis("ID").reset_index()
CAST=CAST.rename_axis("ID").reset_index()
print(OTNine)
print(CAST)
OTNine.columns = OTNine.columns.str.replace(r'_Male_OTNine', '')
CAST.columns = CAST.columns.str.replace(r'_Female_CAST', '')
print(OTNine)
print(CAST)


df = pd.DataFrame(columns=["ID"])
df2= pd.DataFrame(columns = ["ID"])

for col in OTNine.columns[1:]:
   count= OTNine[['ID',col]].merge(CAST[['ID',col]], on='ID', suffixes=('_OTNine','_CAST'))
   a=list(count.columns)
   #header=count.columns[1].split("_")
   #header=str(header[0])
   count["count"]=count[a[1]]+count[a[2]]
   count=count[count["count"]>0]
   count=count.drop(['count'], axis = 1)
   allele=count[["ID", count.columns[2], count.columns[1]]]
   count=count.fillna(0)
   df=df.merge(allele, on=["ID"], how="outer")
   header=count.columns[1].split("_")
   header=str(header[0])
   c1 = (5 < count[count.columns[1]]) & (count[count.columns[1]] < 95 )
   c2 = count[count.columns[1]] <= 5
   count[header] = np.select([c1, c2], ['Biallelic', 'CAST'], 'OTNine')
   count=count[["ID", header]]
   df2 = df2.merge(count, on=["ID"], how='outer')


df=df.fillna(0)
df.to_csv("day0_combine_10min_cell_gene_exp.txt", sep="\t",index=False)
df2=df2.fillna(0)
#print(df2)
dff = df2.copy()
dff['Biallelic']=(df2 == 'Biallelic').T.sum()
dff['OTNine']=(df2 == 'OTNine').T.sum()
dff['CAST']=(df2 == 'CAST').T.sum()
dff["Total"]=dff['Biallelic']+dff['OTNine']+dff['CAST']
dff["Biallelic_perc"]=dff['Biallelic']/dff["Total"]*100
dff["OTNine_perc"]=dff['OTNine']/dff["Total"]*100
dff["CAST_perc"]=dff['CAST']/dff["Total"]*100
#total allele kinetics file include all category
dff.to_csv("day0_allele_kinetics.txt", sep="\t", index=False)
#Cat 1 below pattern we should get 
####################################################################
#  Cat1                                                            #
#           OTNine OTNine OTNine OTNine OTNine OTNine OTNine OTNine OTNine OTNine                #
#                                    or                            #
#           CAST CAST CAST CAST CAST CAST CAST CAST                #
####################################################################
mono_OTNine=dff[dff["OTNine_perc"] == 100]#mono OTNine
mono_OTNine.to_csv("day0_mono_OTNine_cat1.txt", sep="\t", index=False)
mono_CAST=dff[dff["CAST_perc"] == 100]#mono CAST
mono_CAST.to_csv("day0_mono_CAST_cat1.txt", sep="\t", index=False)
Cat1=mono_OTNine.merge(mono_CAST, on=["ID"], how='outer')
Cat1=Cat1.fillna(0)
Cat1.to_csv("Cat1", sep="\t", index=False)
# for Pure Biallelic 
####################################################################################################
#  Cat1.5(>95%)                                                                                    #
#                                                                                                  #
#  Biallelic Biallelic Biallelic Biallelic Biallelic Biallelic Biallelic Biallelic Biallelic OTNine   #
####################################################################################################    
pure_Biallelic=dff[dff["Biallelic_perc"]>95]
pure_Biallelic.to_csv("Biallelic_Cat4_above_95_perc.txt", sep="\t", index=False)
#for cat2 remove the above all genes in the main dff which will give gines with cat2 and cat3
remove_genes=list(mono_OTNine["ID"])+ list(mono_CAST["ID"])+ list(pure_Biallelic["ID"])
print("Mono OTNine "+str(len(list(mono_OTNine["ID"]))))
print("Mono CAST "+ str(len(list(mono_CAST["ID"]))))
print("Mono Biallelic  "+ str(len(list(pure_Biallelic["ID"]))))
print("Removeable genes to get CAT2 and CAT3  "+ str(len(remove_genes)))
dff2=dff[~dff["ID"].isin(remove_genes)]
#Cat 2 get either zero perc from OTNine or CAST column which will have below patern for cat2
###############################################################
#  Cat 2                                                      #
#    OTNine OTNine Biallelic Biallelic OTNine OTNine OTNine Biallelic        #
#                          or                                 #
#      CAST CAST Biallelic CAST Biallelic CAST CAST           # 
###############################################################
dff2_OTNine=dff2[dff2["OTNine_perc"]==0]
dff2_OTNine.to_csv("day0_cat2_only_CAST.txt", sep="\t", index=False )
dff2_CAST=dff2[dff2["CAST_perc"]==0]
dff2_CAST.to_csv("day0_cat2_only_OTNine.txt", sep="\t", index=False )
Cat2=dff2_OTNine.merge(dff2_CAST, on=["ID"], how="outer")
Cat2=Cat2.fillna(0)
Cat2.to_csv("Cat2", sep="\t", index=False)
# Cat 3 remove the Cat2 genes from dff2 rest will be the CAT3 genes 
####################################################################
# Cat3                                                             #
#      CAST OTNine Biallelic CAST OTNine OTNine CAST OTNine                    #
#                                                                  #
####################################################################
Cat2list=list(Cat2["ID"])
print("Removeable genes to get CAT3 "+ str(len(list(Cat2["ID"]))))
Cat3=dff2[~dff2["ID"].isin(Cat2list)]
Cat3.to_csv("Cat3", sep="\t", index=False)




