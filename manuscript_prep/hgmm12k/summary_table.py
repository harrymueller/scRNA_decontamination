import pandas as pd
import os

INPUT="/data/Perkins/Decontamination Project/Analyses/hgmm12k"
FILES=["soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "decontx:paper", "fastcar", "cellbender"]
NEW_NAMES=["SoupX w autoEstCont", "SoupX w Backgorund Genes", "SoupX w Top Background Genes", "DecontX w No Cell Types", "DecontX w Cell Types", "DecontX w Paper Param", "FastCar", "CellBender"]
OUTPUT=os.path.join(INPUT, "origin_summary.xlsx")

# Excel writer
writer = pd.ExcelWriter(OUTPUT)

# template
main_df = pd.DataFrame(columns = ["celltype", "measure", "hg19", "mm10", "all", "method"])

# get data
for (i,f) in enumerate(FILES):
    df = pd.read_excel(os.path.join(INPUT, f, "transcript_origins.xlsx"), sheet_name = "summary", header = 1, index_col = None)
    df["method"] = f

    for a in ["exogenous percentage removed", "endogenous percentage removed"]:
        main_df = pd.concat([main_df, df[df.celltype.str.contains(a)]])

df = main_df

# create two dataframes, one for medians and one for std
a = df[df["measure"] == "median"].reset_index(drop = True)
b = df[df["measure"] == "std"].reset_index(drop = True)

main_df = a[["method", "celltype"]]
for n in ["hg19", "mm10", "all"]:
    main_df[n] = a[n]
    main_df[n + " sd"] = b[n]

names = ["method", "hg19", "hg19 sd", "mm10", "mm10 sd", "all", "all sd"]

# header row
for i in ["exogenous", "endogenous"]:
    df = pd.DataFrame(["Percentage (%) of " + i.capitalize() + " UMIs Removed"])
    df.to_excel(writer, sheet_name = i.capitalize(), index = False, header = False, startrow = 0)

    df = main_df[main_df.celltype.str.contains(i)][names]
    df.to_excel(writer, sheet_name = i.capitalize(), index = False, header = True, startrow = 1)


#print(df.celltype.str.contains("exogenous"))
writer.close()