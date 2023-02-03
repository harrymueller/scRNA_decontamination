import pandas as pd
#import openpyxl
import os

INPUT="/data/Perkins/Decontamination Project/Analyses/mouse kidney/reclus/tables"
FILES=["soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "fastcar", "cellbender"]
NEW_NAMES=["SoupX w autoEstCont", "SoupX w Backgorund Genes", "SoupX w Top Background Genes", "DecontX w No Cell Types", "DecontX w Cell Types", "FastCar", "CellBender"]
OUTPUT=os.path.join(INPUT, "all.xlsx")

#def get_df(num,)

# Excel writer
writer = pd.ExcelWriter(OUTPUT)

# HEADER
df = pd.read_excel(os.path.join(INPUT, FILES[0] + ".xlsx"), sheet_name = 1, nrows=5, header = None)
df.to_excel(writer, sheet_name = "Header", index = False, header = False)

# Contingency Tables
for (i,f) in enumerate(FILES):
    # name
    df = pd.DataFrame([NEW_NAMES[i],"Fresh"])
    df.to_excel(writer, sheet_name = NEW_NAMES[i], index = False, header = False, startrow = 0)

    df = pd.read_excel(os.path.join(INPUT, f + ".xlsx"), sheet_name = "Fresh", header = 6, index_col = 0)
    df.to_excel(writer, sheet_name = NEW_NAMES[i], startrow = 2)
    l = df.shape[0]+4
    df = pd.DataFrame(["MeOH"])
    df.to_excel(writer, sheet_name = NEW_NAMES[i], index = False, header = False, startrow = l)
    df = pd.read_excel(os.path.join(INPUT, f + ".xlsx"), sheet_name = "MeOH", header = 6, index_col = 0)
    df.to_excel(writer, sheet_name = NEW_NAMES[i], startrow = l+1)

writer.close()