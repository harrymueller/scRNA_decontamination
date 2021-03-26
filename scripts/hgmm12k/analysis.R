


### identify transcript origin
### tally transcripts by ct
### express as percentage
### CI?
identify_transcript_origin <- function (samples.combined) {
    # list of mm genes + list of hg genes
    # => matrix subsets containing only mm || hg genes

    # for each subset: for each cell: sum endo transcripts against non-endo transcripts -> calculate fraction of endo
    # DF (barcode, ct, endo_counts, non_endo_counts, endo_frac)

    # by ct: average fractions + sd
    
    # save to excel: first sheet is summary, then full DF (separated by CT)

    # return DF
}

### human gene counts in mouse cells w/ mouse genes x1000 on y && reverse
plot_transcripts <- function(transcripts) {
    # dotplot of transcripts

}


