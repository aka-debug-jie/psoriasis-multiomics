script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "../scripts/mr_core.R"))
source(file.path(dirname(script), "../scripts/mr_available.R"))
root <- tempfile("mr-fixture-")
dir.create(file.path(root, "research/data/processed/exposures"), recursive = TRUE)
dir.create(file.path(root, "research/data/processed/outcomes"), recursive = TRUE)
out <- file.path(root, "output"); dir.create(out)
x <- data.frame(cohort=c("Iceland","Iceland","UKB-PPP"), assay=c("a","a","b"), gene="GENE", uniprot="P1",
                variant=c("rs1","rs2","rs1"), chr="1", pos=c(100,200,100), EA=c("A","G","A"), OA=c("C","A","C"),
                eaf=.2, beta=c(.5,1,.25), neg_log10_p=c(20,40,20), build="GRCh38", se_reconstructed=c(.1,.1,.05))
y <- data.frame(rsid=c("rs1","rs2"), EA=c("A","A"), OA=c("C","G"), beta=c(.1,-.02),se=c(.05,.1),eaf=c(.2,.8),
                p_reported=c(.001,.5),chr="1",pos=c(100,200),build="GRCh38",se_source="reported",usable_effect=TRUE)
write.table(x,file.path(root,"research/data/processed/exposures/eldjarn_all_cis_unselected.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
for (name in c("mvp","finngen_r12")) write.table(y,file.path(root,"research/data/processed/outcomes",paste0(name,"_normalized.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
r <- mr_available(root,list(discovery=list(reported_p_threshold=1.89e-5)),out)
z <- read.delim(file.path(out,"cohort_wald.tsv"))
stopifnot(nrow(z)==2,z$rsid[z$cohort=="Iceland"]=="rs2",z$beta_mr[z$cohort=="Iceland"]==.02,
          !r$formal_candidate_list_produced,r$proteins==1,r$finngen_proteins==1)
# Exact duplicates must not double the precision of the replication result.
write.table(rbind(y,y[1,]),file.path(root,"research/data/processed/outcomes/finngen_r12_normalized.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
r2 <- mr_available(root,list(discovery=list(reported_p_threshold=1.89e-5)),out)
rep <- read.delim(file.path(out,"finngen_same_instrument_wald.tsv"))
stopifnot(nrow(rep)==2)
# Conflicting duplicates exclude that selected instrument, not average it.
bad <- y[1,]; bad$beta <- .3
write.table(rbind(y,bad),file.path(root,"research/data/processed/outcomes/finngen_r12_normalized.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
r3 <- mr_available(root,list(discovery=list(reported_p_threshold=1.89e-5)),out)
rep <- read.delim(file.path(out,"finngen_same_instrument_wald.tsv"))
stopifnot(nrow(rep)==1,rep$rsid=="rs2")
cat("PASS: exposure-only index selection, swapped outcome allele, fixed instruments for replication, partial-coverage status\n")
