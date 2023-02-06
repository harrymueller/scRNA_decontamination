import pandas as pd
import numpy as np
#import openpyxl
import os

INPUT="/data/Perkins/Decontamination Project/Analyses/mouse kidney/reclus"
FILES=["soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "fastcar", "cellbender"]
NEW_NAMES=["SoupX w autoEstCont", "SoupX w Background Genes", "SoupX w Top Background Genes", "DecontX w No Cell Types", "DecontX w Cell Types", "FastCar", "CellBender"]
OUTPUT=os.path.join(INPUT, "prior_post.xlsx")
NUM_CT = list(range(1, 14))

# Excel writer
writer = pd.ExcelWriter(OUTPUT)

# separately for fresh and methanol fixed samples
for fm in ["fresh", "MeOH"]:
    # create summary dataframe for number of genes DE in at least n cell types
    summary = pd.DataFrame(np.zeros([len(NUM_CT),len(NEW_NAMES)+1], dtype = int), columns = ["Number of CT"] + NEW_NAMES)
    summary["Number of CT"] = summary.index + 1

    # main df for genes
    main_df = pd.DataFrame()

    # for each method, read data and rename columns
    for (i,f) in enumerate(FILES):
        df = pd.read_excel(os.path.join(INPUT, f, fm + "_DEGs_between_prior_post.xlsx"), sheet_name = "summary", index_col = 0)
        df = df.rename(columns = {"greater_than_0.5": NEW_NAMES[i]})

        # if some DEGs, join with main dataframe
        if "No.DEGs...0.5......0.5" not in df.columns:
            main_df = pd.concat([main_df, df], axis = 1)
            
            # for each number of cell types, count the number of genes DE in at least that many CT
            for n in NUM_CT:
                summary[NEW_NAMES[i]][n-1] = np.sum(df[NEW_NAMES[i]]>=n)
        else:
            main_df[NEW_NAMES[i]] = np.NaN

    # SAVING
    # header
    df = pd.DataFrame([fm.capitalize() + " - Number of genes DE in at least n cell types."])
    df.to_excel(writer, sheet_name = fm.capitalize() + " Summary", index = False, header = False, startrow = 0)

    # summary
    print(summary)
    summary.to_excel(writer, sheet_name = fm.capitalize() + " Summary", startrow = 2, index = False)

    # header
    df = pd.DataFrame([fm.capitalize() + " - Number of cell types that a given genes is differentially expressed in."])
    df.to_excel(writer, sheet_name = fm.capitalize(), index = False, header = False, startrow = 0)

    # main data
    main_df = main_df[NEW_NAMES].fillna(0).astype(int).astype(str).replace("0","-")
    main_df.to_excel(writer, sheet_name = fm.capitalize(), startrow = 2)

writer.close()