import pandas as pd

dir = "/data/Perkins/Mouse_Kidney/results/reclus"
methods = ["no_decontamination", "soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "fastcar", "decontx:with_cell_types", "decontx:no_cell_types", "cellbender"]
cell_types = ["Podo", "PT", "MPH", "Endo", "B", "CD_IC", "aLOH", "MC", "CD_PC", "Fib", "DCT_CNT", "T", "NK"]

INCLUDE_NAMES = True

gt = {"cell_types": cell_types}
lt = {"cell_types": cell_types}

default_genes_gt = {}
default_genes_lt = {}

for m in methods:
    gt[m] = []
    lt[m] = []
    for ct in cell_types:
        x = pd.read_excel("{}/{}/DEGs.xlsx".format(dir, m), sheet_name=ct)
        x["avg_logFC"] = x["avg_logFC"].astype(float)
        
        gt_mask = x["avg_logFC" ] > 1
        lt_mask = x["avg_logFC" ] < 1
        if not INCLUDE_NAMES:
            gt[m].append(len(x[gt_mask]))
            lt[m].append(len(x[lt_mask]))
        elif INCLUDE_NAMES and m == "no_decontamination":
            gt[m].append(len(x[gt_mask]))
            lt[m].append(len(x[lt_mask]))

            default_genes_gt[ct] = list(x[gt_mask].iloc[:,0])
            default_genes_lt[ct] = list(x[lt_mask].iloc[:,0])
        else:
            gt_genes = x[gt_mask].iloc[:,0]
            gt_genes = [gene for gene in gt_genes if gene not in default_genes_gt[ct]]

            lt_genes = x[lt_mask].iloc[:,0]
            lt_genes = [gene for gene in lt_genes if gene not in default_genes_lt[ct]]

            if len(gt_genes) > 0: gt[m].append("{} (+ {})".format(len(x[gt_mask]), ", ".join(gt_genes)))
            else: gt[m].append(len(x[gt_mask]))

            if len(lt_genes) > 0: lt[m].append("{} (+ {})".format(len(x[lt_mask]), ", ".join(lt_genes)))
            else: lt[m].append(len(x[lt_mask]))

dfs = [pd.DataFrame(gt), pd.DataFrame(lt)]

with pd.ExcelWriter('output.xlsx') as writer:
    dfs[0].to_excel(writer, sheet_name = "higher expr in Fresh - logFC greater than 1", index = False)
    dfs[1].to_excel(writer, sheet_name = "higher expr in MeOH - logFC less than -1", index = False)