# ============================================================
# fGSEA (Hallmark: MSigDB GMT / GO_BP-CC-MF: OrgDb 기반 직접 매핑)
# + KEGG(gseKEGG) + Reactome(gsePathway) + FerrDb
# + 모든 세트에 대해 Ferro-only CSV/3-panel PNG+SVG 저장
# ============================================================

suppressPackageStartupMessages({
  library(readxl); library(stringr)
  library(fgsea); library(ggplot2); library(enrichplot)
  library(clusterProfiler); library(ReactomePA)
  library(org.Mm.eg.db)
  library(AnnotationDbi); library(GO.db)
  library(dplyr)
})

# -----------------------------
# 1) 경로 설정
# -----------------------------
deg_path  <- "/Users/minsu/Downloads/ASGSR2025/poster/OSD-699_DEG/DEanalysis/CA_GCvsFLT-SAL.xlsx"
gmt_dir   <- "/Users/minsu/Downloads/MSigDB/Mouse"
out_root  <- "/Users/minsu/Downloads/ASGSR2025/poster/RRRM-2/OSD-562_DEG/fGSEA/FLT_OLD_ISS_vs_GC_OLD_ISS"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2) DEG 불러오기 및 랭킹 벡터 생성
# -----------------------------
df <- read_xlsx(deg_path)
col_symbol <- names(df)[grepl("symbol", names(df), ignore.case=TRUE)][1]
col_lfc    <- names(df)[grepl("log2",   names(df), ignore.case=TRUE)][1]
col_pval   <- names(df)[grepl("p.value|p_value|pvalue", names(df), ignore.case=TRUE)][1]
if (any(is.na(c(col_symbol, col_lfc, col_pval))))
  stop("SYMBOL/log2FC/p-value 열을 찾을 수 없습니다.")

df2 <- df %>%
  transmute(SYMBOL = as.character(.data[[col_symbol]]),
            log2FC = as.numeric(.data[[col_lfc]]),
            pval   = as.numeric(.data[[col_pval]])) %>%
  filter(!is.na(SYMBOL), !is.na(log2FC), !is.na(pval))

df2 <- df2 %>% mutate(rank_stat = sign(log2FC) * -log10(pval + 1e-300))
geneList <- tapply(df2$rank_stat, df2$SYMBOL, function(x) x[which.max(abs(x))])
geneList <- sort(geneList, decreasing=TRUE)
geneList <- setNames(as.numeric(geneList), names(geneList))
message(sprintf("✅ Ranked genes ready: %d (unique SYMBOLs)", length(geneList)))

entrez_map <- bitr(names(geneList), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
rank_entrez <- geneList[entrez_map$SYMBOL]
names(rank_entrez) <- entrez_map$ENTREZID
rank_entrez <- sort(rank_entrez[!is.na(rank_entrez)], decreasing=TRUE)

# -----------------------------
# 3) Ferroptosis 관련 키워드
# -----------------------------
ferro_keywords <- c(
  "ferroptosis","iron","heme","ferritin","transferrin",
  "oxidative stress","oxidative_stress","oxidative-stress","reactive oxygen","reactive_oxygen","reactive-oxygen",
  "nitric oxide","nitric_oxide","nitric-oxide","lipid","peroxidation",
  "glutathione","redox","hydrogen peroxide","hydrogen_peroxide","hydrogen-peroxide",
  "superoxide","mitochond","hypoxia","autophagy", "HIF"
)

# -----------------------------
# 4) gseaResult 저장 + Ferro-only 하위결과/그림
# -----------------------------
save_gsea_results <- function(gsea_obj, label, out_dir) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  df <- as.data.frame(gsea_obj@result)
  if (nrow(df) == 0) { message("⚠️ No terms for ", label); return(invisible(NULL)) }
  df <- df %>% arrange(desc(NES))
  if ("pval" %in% names(df)) names(df)[names(df)=="pval"] <- "pvalue"
  if ("padj" %in% names(df)) names(df)[names(df)=="padj"] <- "p.adjust"
  
  out_csv <- file.path(out_dir, paste0(label, "_fGSEA_results.csv"))
  write.csv(df, out_csv, row.names=FALSE)
  message("✅ Saved: ", out_csv)
  
  ferro_idx <- Reduce(`|`, lapply(ferro_keywords, function(k) grepl(k, tolower(df$Description))))
  if (any(ferro_idx)) {
    ferro_res <- df[ferro_idx, ]
    ferro_csv <- file.path(out_dir, paste0(label, "_FERRO_related.csv"))
    write.csv(ferro_res, ferro_csv, row.names=FALSE)
    message("🔎 Ferro subset saved → ", ferro_csv)
    
    ferro_dir <- file.path(out_dir, "Ferro_only")
    dir.create(ferro_dir, showWarnings=FALSE)
    for (pid in ferro_res$ID) {
      p <- tryCatch({
        enrichplot::gseaplot2(gsea_obj, geneSetID = pid,
                              subplots = 1:3, pvalue_table = TRUE,
                              title = ferro_res$Description[ferro_res$ID == pid])
      }, error=function(e) NULL)
      if (!is.null(p)) {
        fname <- paste0(label, "_", gsub("[^A-Za-z0-9]+","_", ferro_res$Description[ferro_res$ID==pid]), "_3panel")
        ggsave(file.path(ferro_dir, paste0(fname, ".png")), plot=p, width=6.5, height=6.5, dpi=300)
        ggsave(file.path(ferro_dir, paste0(fname, ".svg")), plot=p, width=6.5, height=6.5)
      }
    }
    message("🎨 Ferro-only plots → ", ferro_dir)
  }
}

# -----------------------------
# 5) OrgDb 기반 GO 세트
# -----------------------------
build_go_sets <- function(ont=c("BP","CC","MF")){
  ont<-match.arg(ont)
  all_keys <- keys(org.Mm.eg.db,keytype="ENTREZID")
  go_map <- AnnotationDbi::select(org.Mm.eg.db,keys=all_keys,
                                  columns=c("GO","ONTOLOGY"),keytype="ENTREZID")
  go_map <- go_map[!is.na(go_map$GO)&go_map$ONTOLOGY==ont,c("GO","ENTREZID")]
  colnames(go_map)<-c("term","gene")
  go_ids<-unique(go_map$term)
  go_name<-AnnotationDbi::select(GO.db,keys=go_ids,columns="TERM",keytype="GOID")
  colnames(go_name)<-c("term","name")
  list(TERM2GENE=go_map,TERM2NAME=go_name)
}

run_GO_fgsea <- function(geneList_entrez,ont=c("BP","CC","MF")){
  ont<-match.arg(ont)
  go_sets<-build_go_sets(ont)
  gsea_obj<-GSEA(geneList=geneList_entrez,TERM2GENE=go_sets$TERM2GENE,
                 TERM2NAME=go_sets$TERM2NAME,minGSSize=10,maxGSSize=500,pvalueCutoff=1)
  if ("pval"%in%names(gsea_obj@result)) names(gsea_obj@result)[names(gsea_obj@result)=="pval"]<-"pvalue"
  if ("padj"%in%names(gsea_obj@result)) names(gsea_obj@result)[names(gsea_obj@result)=="padj"]<-"p.adjust"
  gsea_obj
}

# -----------------------------
# 6) Hallmark (MSigDB GMT) — fgsea + fake_res(pvalue patch) + Ferro_related.csv
# -----------------------------
mh_sets <- fgsea::gmtPathways(file.path(gmt_dir, "mh.all.v2025.1.Mm.symbols.gmt"))
set.seed(42)
fg_mh <- fgseaMultilevel(pathways = mh_sets, stats = geneList, minSize = 10, maxSize = 500,nPermSimple = 1e6)
mh_dir <- file.path(out_root, "Hallmark")
dir.create(mh_dir, showWarnings=FALSE, recursive=TRUE)

fg_mh_df <- as.data.frame(fg_mh)
if ("leadingEdge"%in%names(fg_mh_df))
  fg_mh_df$leadingEdge<-sapply(fg_mh_df$leadingEdge,paste,collapse="; ")
if ("pval"%in%names(fg_mh_df)) names(fg_mh_df)[names(fg_mh_df)=="pval"]<-"pvalue"
if ("padj"%in%names(fg_mh_df)) names(fg_mh_df)[names(fg_mh_df)=="padj"]<-"p.adjust"
fg_mh_df <- fg_mh_df %>%
  mutate(Description=gsub("_"," ",pathway)) %>%
  dplyr::select(Description,everything(),-pathway)
write.csv(fg_mh_df,file.path(mh_dir,"Hallmark_fGSEA_results.csv"),row.names=FALSE)

ferro_idx_mh <- Reduce(`|`,lapply(ferro_keywords,function(k) grepl(k,tolower(fg_mh_df$Description))))
if (any(ferro_idx_mh)){
  ferro_res <- fg_mh_df[ferro_idx_mh,]
  ferro_csv <- file.path(mh_dir, "Hallmark_FERRO_related.csv")
  write.csv(ferro_res, ferro_csv, row.names=FALSE)
  message("🔎 Saved: Hallmark_FERRO_related.csv")
  
  ferro_dir<-file.path(mh_dir,"Ferro_only");dir.create(ferro_dir,showWarnings=FALSE)
  for(i in seq_len(nrow(ferro_res))){
    desc<-ferro_res$Description[i]; pname<-gsub(" ","_",desc)
    if (pname%in%names(mh_sets)){
      tmp_term<-data.frame(term=pname,gene=mh_sets[[pname]])
      fake_res<-suppressMessages(clusterProfiler::GSEA(
        geneList=geneList,TERM2GENE=tmp_term,minGSSize=1,maxGSSize=1e5,pvalueCutoff=1))
      # --- CSV 값 반영 ---
      fake_res@result$pvalue<-ferro_res$pvalue[i]
      fake_res@result$p.adjust<-ferro_res$p.adjust[i]
      fake_res@result$NES<-ferro_res$NES[i]
      p<-enrichplot::gseaplot2(fake_res,geneSetID=pname,subplots=1:3,
                               pvalue_table=TRUE,title=desc)
      fname<-paste0("Hallmark_",gsub("[^A-Za-z0-9]+","_",pname),"_3panel")
      ggsave(file.path(ferro_dir,paste0(fname,".png")),plot=p,width=6.5,height=6.5,dpi=300)
      ggsave(file.path(ferro_dir,paste0(fname,".svg")),plot=p,width=6.5,height=6.5)
    }
  }
  message("🎨 Hallmark Ferro-only plots → ", ferro_dir)
}

# -----------------------------
# 7) GO BP/CC/MF (OrgDb 기반)
# -----------------------------
fg_gobp<-run_GO_fgsea(rank_entrez,"BP"); save_gsea_results(fg_gobp,"GO_BP",file.path(out_root,"GO_BP"))
fg_gocc<-run_GO_fgsea(rank_entrez,"CC"); save_gsea_results(fg_gocc,"GO_CC",file.path(out_root,"GO_CC"))
fg_gomf<-run_GO_fgsea(rank_entrez,"MF"); save_gsea_results(fg_gomf,"GO_MF",file.path(out_root,"GO_MF"))

# -----------------------------
# 8) KEGG
# -----------------------------
ekegg<-gseKEGG(geneList=rank_entrez,organism="mmu",minGSSize=10,maxGSSize=500,pvalueCutoff=1)
ekegg<-setReadable(ekegg,OrgDb=org.Mm.eg.db,keyType="ENTREZID")
save_gsea_results(ekegg,"KEGG",file.path(out_root,"KEGG"))

# -----------------------------
# 9) Reactome
# -----------------------------
eReact<-gsePathway(geneList=rank_entrez,organism="mouse",minGSSize=10,maxGSSize=500,pvalueCutoff=1)
eReact<-setReadable(eReact,OrgDb=org.Mm.eg.db,keyType="ENTREZID")
save_gsea_results(eReact,"Reactome",file.path(out_root,"Reactome"))

# -----------------------------
# 10) FerrDb (fake_res + pvalue patch + dummy gene 보정)
# -----------------------------
ferrdb_dir<-file.path(out_root,"FerrDb_GSEA");dir.create(ferrdb_dir,showWarnings=FALSE)
ferrdb_files<-list(
  Marker = "/Users/minsu/Downloads/ASGSR2025/poster/FerrDb/FerrDb_Marker.csv",
  Driver = "/Users/minsu/Downloads/ASGSR2025/poster/FerrDb/FerrDb_Driver.csv",
  Suppressor = "/Users/minsu/Downloads/ASGSR2025/poster/FerrDb/FerrDb_Suppressor.csv"
)
normalize_symbol<-function(x){sapply(x,function(g){g<-tolower(g);paste0(toupper(substr(g,1,1)),substr(g,2,nchar(g)))})}
ferrdb_sets<-lapply(ferrdb_files,function(f) unique(normalize_symbol(read.csv(f)[,2])))
cat("FerrDb gene counts after deduplication:\n")
sapply(ferrdb_sets, length)

ferrdb_pathways<-list(FerrDb_Marker=ferrdb_sets$Marker,FerrDb_Driver=ferrdb_sets$Driver,FerrDb_Suppressor=ferrdb_sets$Suppressor)

fg_ferrdb<-fgseaMultilevel(pathways=ferrdb_pathways,stats=geneList,minSize=5,maxSize=1000,nPermSimple = 1e6)
fg_ferrdb_df<-as.data.frame(fg_ferrdb)
if ("leadingEdge"%in%names(fg_ferrdb_df)) fg_ferrdb_df$leadingEdge<-sapply(fg_ferrdb_df$leadingEdge,paste,collapse="; ")
if ("pval"%in%names(fg_ferrdb_df)) names(fg_ferrdb_df)[names(fg_ferrdb_df)=="pval"]<-"pvalue"
if ("padj"%in%names(fg_ferrdb_df)) names(fg_ferrdb_df)[names(fg_ferrdb_df)=="padj"]<-"p.adjust"
write.csv(fg_ferrdb_df,file.path(ferrdb_dir,"FerrDb_fGSEA_results.csv"),row.names=FALSE)

# ES plot 복원 + CSV값 반영
for(i in seq_len(nrow(fg_ferrdb_df))){
  pname<-fg_ferrdb_df$pathway[i]
  genes<-ferrdb_pathways[[pname]]
  if (length(genes)<2) genes<-c(genes,"DUMMYGENE")  # dummy gene 보정
  tmp_term<-data.frame(term=pname,gene=genes)
  fake_res<-suppressMessages(clusterProfiler::GSEA(
    geneList=geneList,TERM2GENE=tmp_term,minGSSize=1,maxGSSize=1e5,pvalueCutoff=1))
  fake_res@result$pvalue<-fg_ferrdb_df$pvalue[i]
  fake_res@result$p.adjust<-fg_ferrdb_df$p.adjust[i]
  fake_res@result$NES<-fg_ferrdb_df$NES[i]
  p<-enrichplot::gseaplot2(fake_res,geneSetID=pname,subplots=1:3,
                           pvalue_table=TRUE,title=pname)
  fname<-paste0("FerrDb_",gsub("[^A-Za-z0-9]+","_",pname),"_3panel")
  ggsave(file.path(ferrdb_dir,paste0(fname,".png")),plot=p,width=6.5,height=6.5,dpi=300)
  ggsave(file.path(ferrdb_dir,paste0(fname,".svg")),plot=p,width=6.5,height=6.5)
}
message("🎨 FerrDb plots saved → ", ferrdb_dir)

message("\n🎯 fGSEA 완료 (Hallmark/FerrDb plot + Ferro_related.csv + pvalue box 일치)")
message("All outputs → ", out_root)