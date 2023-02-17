import pandas as pd


def Gene_filter(Mat, Pat, min_cell_exp):
   df_cells=Mat.add(Pat.select_dtypes(include="number")).combine_first(Mat)
   df_cells=df_cells.fillna(0)
   cells=len(list(df_cells)) - min_cell_exp
   df_cells["Sum"]=(df_cells == 0).sum(1)
   df_cells=df_cells[df_cells["Sum"]<=cells]
   df_cells=df_cells.drop(["Sum"], axis=1)
   return(df_cells)


#129 allele
df_B6=pd.read_csv("Day0129rpkm.tsv",sep="\t")
df_B6=df_B6.reset_index()
df_B6.rename({"index":"ID"}, inplace=True, axis=1)
df_CAST=pd.read_csv("Day0castrpkm.tsv",sep="\t")
df_CAST=df_CAST.reset_index()
df_CAST.rename({"index":"ID"}, inplace=True, axis=1)

print(df_B6)
print(df_CAST)


df_add=df_B6.add(df_CAST.select_dtypes(include="number")).combine_first(df_B6)
df_add=df_add.fillna(0)
df_add["Average"]=df_add.mean(axis=1)
df_add=df_add[df_add["Average"]>10] #zero average genes 
df_B6=df_B6[df_B6["ID"].isin(list(df_add["ID"]))]
df_CAST=df_CAST[df_CAST["ID"].isin(list(df_add["ID"]))]

print(df_B6)
print(df_CAST)


Min_cell_exp=Gene_filter(df_B6, df_CAST, 26)#6 is return includes first column
df_B6=df_B6[df_B6["ID"].isin(list(Min_cell_exp["ID"]))]
df_CAST=df_CAST[df_CAST["ID"].isin(list(Min_cell_exp["ID"]))]

df_B6.to_csv("129S1_filtered_grt_10_f25cells_exp.tsv", sep="\t", index=False)
df_CAST.to_csv("CAST_filtered_grt_10_f25cells_exp.tsv", sep="\t", index=False)

print(df_B6)
print(df_CAST)

