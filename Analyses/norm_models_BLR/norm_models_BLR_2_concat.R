
library(data.table)
library(RMINC)

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM", "jacobians_abs", "jacobians_rel")
# names = c("FA")
# names = c("jacobians_abs", "jacobians_rel")


tissues = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')


for (n in 1:length(names)) {
    print(names[n])

    for (t in 1:length(tissues)) {
        print(tissues[t])

        # Find all nm results files for this micro and tissue type
        files = list.files(path="./results", pattern=paste0("nm_",names[n],"_",tissues[t],"_vox_*"), full.names=TRUE)
        end_vox = as.numeric(sub(".*_(\\d+)\\..*", "\\1", files))
        files = files[order(end_vox)]

        # load first file
        nm = fread(files[1])
        print(files[1])
        print(dim(nm))
        for (i in 2:length(files)) {
            print(files[i])
            print(dim(fread(files[i])))
            nm = rbind(nm, fread(files[i]))
        }
        print(dim(nm))

        # Write non-chunked nm for micro and tissue type
        fwrite(nm, paste0("./results/nm_",names[n],"_",tissues[t],".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

        # Delete chunk results
        file.remove(files)
    }
}


# Fix ROI columns
for (n in 1:length(names)) {
    print(names[n])

    for (t in 1:length(tissues)) {
        print(tissues[t])
        nm = as.data.frame(fread(paste0("./results/nm_",names[n],"_",tissues[t],".tsv")))
        nm$ROI = seq(1, nrow(nm))
        fwrite(nm, paste0("./results/nm_",names[n],"_",tissues[t],".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    }
}


