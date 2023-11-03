
library(data.table)

names = c("FA", "MD", "ICVF", "ISOVF", "OD", "T2star", "QSM")
# tissues = c('Cerebellum_GM', 'Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')
tissues = c('Cerebellum_WM', 'Brainstem', 'Subcortical_GM', 'Cortical_GM', 'Cerebral_NAWM')

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
        for (i in 2:length(files)) {
            nm = rbind(nm, fread(files[1]))
        }

        # Delete chunk results
        file.remove(files)

        # Write non-chunked nm for micro and tissue type
        fwrite(nm, paste0("./results/nm_",names[n],"_",tissues[t],".tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    }
}


# Check how many NAs...
files = list.files(path="./results", pattern=paste0("nm_",names[n],"*"), full.names=TRUE)

list_nm = list()

for (i in 1:length(files)) {
    list_nm[[i]] = as.data.frame(fread(files[i]))
    print(summary(list_nm[[i]]$N))
}
