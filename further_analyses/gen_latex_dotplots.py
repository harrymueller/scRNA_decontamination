
FILES=["soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types", "fastcar", "cellbender"]
NEW_NAMES=["SoupX with autoEstCont", "SoupX with Three Marker Genes per Cell Type", "SoupX with the Top 25 Marker Genes", "DecontX without Cell Types", "DecontX with Cell Types", "FastCar", "CellBender"]

text = """
    % ##################################################
    % FULLNAME
    % ##################################################
    \section{FULLNAME}
    \subsection{\\nine}
    \\begin{figure}[h]
        \centering
            
        \caption{\dotplotA}
        \includegraphics[width=\imgwidth]{Nine/DEG_ABBR_Overexpr.png}
            
        \caption{\dotplotB}
        \includegraphics[width=\imgwidth]{Nine/DEG_ABBR_Underexpr.png}
    \end{figure}
    \celltypes

    \\newpage
    \subsection{\specific}
    \\begin{figure}[h]
        \centering
        \caption{\dotplotC}
        \includegraphics[width=\imgwidth]{Specific/DEG_ABBR_Overexpr.png}
            
        \caption{\dotplotD}
        \includegraphics[width=\imgwidth]{Specific/DEG_ABBR_Underexpr.png}
    \end{figure}
    \celltypes
"""

all = ""
for i in range(len(FILES)):
    all += text.replace("ABBR", FILES[i]).replace("FULLNAME", NEW_NAMES[i])

with open("file.txt", "w") as f:
    f.write(all)