# =============================================================================
# GSE1009 — Nefropatia Diabética vs. Controle (glomérulos isolados)
# Plataforma: Affymetrix HG-U95Av2 (GPL8300)
# Versão CORRIGIDA e comentada do pipeline
# =============================================================================
# Mudanças principais em relação à versão original (justificativa completa no
# documento de auditoria):
#   (1) Bloco de debug gráfico removido.
#   (2) Extração de metadados robusta + checagens que FALHAM em vez de seguir
#       silenciosamente com NA.
#   (3) Tratamento explícito das réplicas 1a/1b (não-independência das amostras).
#   (4) Cortes de significância definidos UMA vez e reutilizados em todo lugar
#       (limma, volcano, heatmap, venn, GO/KEGG).
#   (5) unique() nas listas de genes e no universo do enriquecimento.
#   (6) Guardas de NULL/vazio em todos os passos que podem não retornar nada.
#   (7) Remoção de comandos redundantes.
# =============================================================================


# =============================================================================
# 0. PARÂMETROS DA ANÁLISE (fonte única de verdade)
# =============================================================================
# Definir os cortes como constantes evita a classe de bug mais comum deste tipo
# de script: a tabela de DEGs usar um corte e o gráfico usar outro.

FDR_CUT   <- 0.05   # corte de p-valor ajustado (BH)
LFC_CUT   <- 1      # corte de |log2 fold change| (1 = 2x)
TOP_HEAT  <- 50     # nº máximo de genes no heatmap

# Como tratar as réplicas dentro do mesmo indivíduo (Control 1a/1b, Diabetes 1a/1b):
#   "duplicateCorrelation" -> mantém os 6 arrays e modela a correlação intra-indivíduo
#   "average"              -> colapsa para 1 array por indivíduo (4 arrays, 2 vs 2)
#   "ignore"               -> comportamento do script original (NÃO recomendado)
ESTRATEGIA_REPLICA <- "duplicateCorrelation"

DIR_DATA <- "data"
DIR_OUT  <- "results"


# =============================================================================
# 1. PACOTES
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs_cran <- c("tidyverse", "pheatmap", "ggrepel")
pkgs_bioc <- c("affy", "GEOquery", "AnnotationDbi", "hgu95av2.db", "limma",
               "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db",
               "ggVennDiagram", "affyPLM")

falta_cran <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
if (length(falta_cran)) install.packages(falta_cran)

falta_bioc <- pkgs_bioc[!sapply(pkgs_bioc, requireNamespace, quietly = TRUE)]
if (length(falta_bioc)) BiocManager::install(falta_bioc)

suppressPackageStartupMessages({
  library(tidyverse)
  library(affy)
  library(GEOquery)
  library(AnnotationDbi)
  library(hgu95av2.db)
  library(limma)
  library(EnhancedVolcano)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggVennDiagram)
})

options(timeout = 600)
dir.create(DIR_DATA, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_OUT,  showWarnings = FALSE, recursive = TRUE)
data_execucao <- Sys.Date()


# =============================================================================
# 2. DOWNLOAD DOS DADOS BRUTOS (.CEL) E DOS METADADOS
# =============================================================================
# Dois downloads distintos e complementares:
#   - getGEOSuppFiles  -> os arquivos CEL (sinal bruto, por sonda)
#   - getGEO           -> os metadados da série (quem é controle, quem é doente)

tar_file <- file.path(DIR_DATA, "GSE1009", "GSE1009_RAW.tar")
if (!file.exists(tar_file)) {
  getGEOSuppFiles("GSE1009", baseDir = DIR_DATA)
}
stopifnot(file.exists(tar_file))

gse_meta <- getGEO("GSE1009", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pdata    <- pData(gse_meta)

# --- CHECAGEM: a plataforma anunciada bate com o pacote de anotação usado? ----
plataforma <- unique(as.character(pdata$platform_id))
message("Plataforma declarada no GEO: ", paste(plataforma, collapse = ", "))
stopifnot(identical(plataforma, "GPL8300"))   # GPL8300 == HG-U95Av2 == hgu95av2.db

# --- Onde está a informação de grupo nesta série? ----------------------------
# Séries de 2003 não têm characteristics_ch1 padronizado. Aqui o grupo está no
# campo de texto livre "description". IMPORTANTE: o SOFT do GSE1009 tem VÁRIAS
# linhas !Sample_description por amostra ("Keywords = Kidney" etc.), e o GEOquery
# transforma isso em description, description.1, description.2. Só a primeira
# carrega o grupo — por isso selecionamos a coluna explicitamente e conferimos.
stopifnot("description" %in% colnames(pdata))
print(table(pdata$description, useNA = "ifany"))
print(pdata[, c("title", "description")])


# =============================================================================
# 3. DESCOMPACTAR E LOCALIZAR OS CELs
# =============================================================================

dir_cel <- file.path(DIR_DATA, "GSE1009", "CEL")
dir.create(dir_cel, showWarnings = FALSE, recursive = TRUE)
untar(tar_file, exdir = dir_cel)

cel_files <- sort(list.files(dir_cel, pattern = "\\.CEL(\\.gz)?$",
                             full.names = TRUE, recursive = TRUE,
                             ignore.case = TRUE))
if (length(cel_files) == 0) stop("Nenhum arquivo .CEL encontrado em ", dir_cel)
message("CELs encontrados: ", length(cel_files))


# =============================================================================
# 4. CASAR ARQUIVOS CEL <-> GSM <-> GRUPO
# =============================================================================
# NUNCA assumir que a ordem alfabética dos arquivos corresponde à ordem do pdata.
# O casamento é feito pelo nome do arquivo declarado em supplementary_file.

col_supp <- grep("^supplementary_file", colnames(pdata), value = TRUE)[1]
if (is.na(col_supp)) stop("Coluna supplementary_file não encontrada no pData.")

limpar_nome <- function(x) sub("\\.CEL$", "", sub("\\.gz$", "", x, ignore.case = TRUE),
                               ignore.case = TRUE)

mapa <- data.frame(
  gsm   = rownames(pdata),
  chave = limpar_nome(basename(as.character(pdata[[col_supp]]))),
  stringsAsFactors = FALSE
)

idx <- match(limpar_nome(basename(cel_files)), mapa$chave)
if (any(is.na(idx))) {
  print(data.frame(cel = basename(cel_files), casou = !is.na(idx)))
  stop("Falha ao casar arquivos CEL com os GSM do pData.")
}

gsm_ids  <- mapa$gsm[idx]
pdata_ord <- pdata[gsm_ids, ]

# --- Vetor de grupo, com falha explícita se algum rótulo não for reconhecido --
niveis_geo <- c("mRNA from isolated glomeruli from normal kidney",
                "mRNA from isolated glomeruli from kidneys with diabetic nephropathy")

texto_grupo <- trimws(as.character(pdata_ord$description))
nao_reconhecido <- setdiff(unique(texto_grupo), niveis_geo)
if (length(nao_reconhecido)) {
  stop("Rótulo de grupo não reconhecido no campo description:\n  ",
       paste(nao_reconhecido, collapse = "\n  "))
}

grupo <- factor(texto_grupo, levels = niveis_geo,
                labels = c("Controle", "NefropatiaDiabetica"))
stopifnot(!any(is.na(grupo)))

# --- Indivíduo (bloco): "Control 1a"/"Control 1b" são o MESMO indivíduo -------
# Os títulos do GEO são: Control 1a, Control 1b, Control 2, Diabetes 1a,
# Diabetes 1b, Diabetes 2. O sufixo a/b indica réplicas do mesmo indivíduo.
# Ignorar isso = pseudo-replicação (ver auditoria, Problema #1).
individuo <- factor(sub("[ab]$", "", trimws(as.character(pdata_ord$title))))

metadados <- data.frame(
  arquivo   = basename(cel_files),
  gsm       = gsm_ids,
  titulo    = as.character(pdata_ord$title),
  individuo = individuo,
  grupo     = grupo,
  stringsAsFactors = FALSE
)
print(metadados)
print(table(metadados$grupo))
print(table(metadados$individuo, metadados$grupo))
write.csv(metadados, file.path(DIR_OUT, "GSE1009_metadados_amostras.csv"),
          row.names = FALSE)


# =============================================================================
# 5. LEITURA DOS ARQUIVOS CEL
# =============================================================================
# Não passamos compress=: o pacote affy detecta .gz sozinho, e forçar um único
# valor quebra se o diretório tiver arquivos comprimidos e não comprimidos.

raw_data <- ReadAffy(filenames = cel_files)
sampleNames(raw_data) <- metadados$gsm      # nomes estáveis e rastreáveis
print(raw_data)

# Data de escaneamento: proxy barato para efeito de lote (batch).
scan_date <- tryCatch(as.character(protocolData(raw_data)$ScanDate),
                      error = function(e) rep(NA_character_, ncol(raw_data)))
print(data.frame(gsm = metadados$gsm, grupo = metadados$grupo, scan = scan_date))


# =============================================================================
# 6. CONTROLE DE QUALIDADE (dados brutos)
# =============================================================================

cores_grupo <- c(Controle = "#2C7FB8", NefropatiaDiabetica = "#D95F02")

pdf(file.path(DIR_OUT, "01_QC_bruto.pdf"), width = 11, height = 7)
boxplot(raw_data, las = 2, col = cores_grupo[metadados$grupo],
        main = "GSE1009 — log2(intensidade) bruta por array")
hist(raw_data, main = "GSE1009 — densidade dos sinais brutos")
deg <- AffyRNAdeg(raw_data)
plotAffyRNAdeg(deg, cols = cores_grupo[metadados$grupo])
title(sub = "Inclinações muito diferentes = degradação desigual de RNA")
dev.off()

# QC quantitativo (melhor que o boxplot): RLE e NUSE.
# Regra prática: NUSE mediano > 1.05 ou RLE mediano longe de 0 => array suspeito.
if (requireNamespace("affyPLM", quietly = TRUE)) {
  plm <- affyPLM::fitPLM(raw_data, background = TRUE, normalize = TRUE)
  pdf(file.path(DIR_OUT, "02_QC_RLE_NUSE.pdf"), width = 10, height = 7)
  affyPLM::RLE(plm,  main = "RLE — Relative Log Expression", las = 2)
  abline(h = 0, lty = 2)
  affyPLM::NUSE(plm, main = "NUSE — Normalized Unscaled Standard Error", las = 2)
  abline(h = 1.05, lty = 2, col = "red")
  dev.off()
}


# =============================================================================
# 7. NORMALIZAÇÃO RMA
# =============================================================================
# RMA = correção de background (modelo normal+exponencial) + normalização por
# quantis (força a MESMA distribuição em todos os arrays) + sumarização por
# median polish em escala log2 (várias sondas -> 1 valor por probe set).

eset <- rma(raw_data)

# Sondas AFFX- são controles do fabricante (spike-ins, housekeeping de controle),
# não genes de interesse. Removidas DEPOIS do RMA de propósito: elas participam
# legitimamente da estimativa de background/quantis, mas não devem entrar no
# teste estatístico (inflam a correção de múltiplos testes).
expr <- exprs(eset)
n_affx <- sum(grepl("^AFFX", rownames(expr)))
expr <- expr[!grepl("^AFFX", rownames(expr)), , drop = FALSE]
message("Sondas AFFX removidas: ", n_affx, " | probe sets restantes: ", nrow(expr))

pdf(file.path(DIR_OUT, "03_boxplot_RMA.pdf"), width = 11, height = 7)
boxplot(as.data.frame(expr), las = 2, col = cores_grupo[metadados$grupo],
        main = "GSE1009 — após RMA (as caixas devem ficar alinhadas)")
dev.off()


# =============================================================================
# 8. PCA
# =============================================================================
# t() porque prcomp espera amostras nas LINHAS e variáveis nas COLUNAS; a matriz
# de expressão vem no formato oposto (genes nas linhas).
# scale. = FALSE porque os dados já estão em log2 e na mesma escala — escalar
# aqui daria peso igual a sondas de ruído e a sondas de sinal forte.

pca <- prcomp(t(expr), scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(pca$x[, 1:2], metadados)

g_pca <- ggplot(pca_df, aes(PC1, PC2, color = grupo, shape = individuo)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = titulo), size = 3, show.legend = FALSE) +
  scale_color_manual(values = cores_grupo) +
  labs(title = "PCA — GSE1009 (RMA, todas as sondas)",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  theme_classic()

ggsave(file.path(DIR_OUT, "04_PCA.pdf"), g_pca, width = 9, height = 7)
write.csv(pca_df, file.path(DIR_OUT, "GSE1009_PCA_coordenadas.csv"), row.names = FALSE)


# =============================================================================
# 9. ANOTAÇÃO DAS SONDAS
# =============================================================================
# Objeto NÃO chamado "annotation": esse nome colide com Biobase::annotation().

anot_probes <- AnnotationDbi::select(
  hgu95av2.db,
  keys     = rownames(expr),
  columns  = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype  = "PROBEID"
)

# Quantas sondas são realmente ambíguas (1 sonda -> vários genes)?
ambiguas <- anot_probes %>% count(PROBEID) %>% filter(n > 1)
message("Sondas com mapeamento 1:muitos: ", nrow(ambiguas),
        " de ", length(unique(anot_probes$PROBEID)))
write.csv(anot_probes %>% semi_join(ambiguas, by = "PROBEID"),
          file.path(DIR_OUT, "GSE1009_probes_ambiguas.csv"), row.names = FALSE)

# Estratégia adotada: manter 1 linha por PROBEID (primeira ocorrência).
# Isso é aceitável AQUI porque o teste estatístico roda no nível de PROBE SET —
# a escolha só afeta o RÓTULO, não o p-valor. Mas as sondas ambíguas ficam
# registradas no CSV acima para conferência dos genes que forem reportados.
anot_1p1 <- anot_probes %>% distinct(PROBEID, .keep_all = TRUE)

taxa_anot <- 100 * mean(!is.na(anot_1p1$SYMBOL))
message(sprintf("Taxa de anotação: %.1f%%", taxa_anot))


# =============================================================================
# 10. EXPRESSÃO DIFERENCIAL (limma)
# =============================================================================
# ~0 + grupo => matriz de médias por grupo (sem intercepto). Cada coluna é a
# média do grupo; o contraste subtrai uma coluna da outra.

design <- model.matrix(~ 0 + grupo)
colnames(design) <- levels(grupo)

contrastes <- makeContrasts(NefropatiaDiabetica - Controle, levels = design)

if (ESTRATEGIA_REPLICA == "average") {
  # Colapsa réplicas do mesmo indivíduo -> 1 array por indivíduo (2 vs 2).
  expr_fit   <- avearrays(expr, ID = as.character(individuo))
  grupo_fit  <- factor(tapply(as.character(grupo), individuo,
                              function(x) x[1])[colnames(expr_fit)],
                       levels = levels(grupo))
  design     <- model.matrix(~ 0 + grupo_fit); colnames(design) <- levels(grupo_fit)
  contrastes <- makeContrasts(NefropatiaDiabetica - Controle, levels = design)
  fit <- lmFit(expr_fit, design)

} else if (ESTRATEGIA_REPLICA == "duplicateCorrelation") {
  # Mantém os 6 arrays mas informa ao limma que arrays do mesmo indivíduo são
  # correlacionados -> os graus de liberdade deixam de ser inflados.
  corfit <- duplicateCorrelation(expr, design, block = individuo)
  message("Correlação intra-indivíduo (consensus): ",
          round(corfit$consensus.correlation, 3))
  fit <- lmFit(expr, design, block = individuo,
               correlation = corfit$consensus.correlation)

} else {
  warning("ESTRATEGIA_REPLICA = 'ignore': arrays do mesmo indivíduo tratados ",
          "como independentes (pseudo-replicação).")
  fit <- lmFit(expr, design)
}

fit2 <- eBayes(contrasts.fit(fit, contrastes))

print(summary(decideTests(fit2, adjust.method = "BH",
                          p.value = FDR_CUT, lfc = LFC_CUT)))

res <- topTable(fit2, adjust.method = "BH", number = Inf, sort.by = "P") %>%
  rownames_to_column("PROBEID") %>%
  left_join(anot_1p1, by = "PROBEID") %>%
  mutate(rotulo = ifelse(is.na(SYMBOL), PROBEID, SYMBOL)) %>%
  relocate(PROBEID, SYMBOL, GENENAME, ENTREZID)

degs <- res %>% filter(adj.P.Val < FDR_CUT, abs(logFC) > LFC_CUT)
message("DEGs (adj.P < ", FDR_CUT, " e |logFC| > ", LFC_CUT, "): ", nrow(degs))

write.csv(as.data.frame(expr), file.path(DIR_OUT, "GSE1009_expressao_RMA.csv"))
write.csv(res,  file.path(DIR_OUT, "GSE1009_DE_completo.csv"), row.names = FALSE)
write.csv(degs, file.path(DIR_OUT, "GSE1009_DE_significativos.csv"), row.names = FALSE)


# =============================================================================
# 11. VOLCANO PLOT
# =============================================================================
# Usa EXATAMENTE FDR_CUT e LFC_CUT — os mesmos do passo 10.
# y = adj.P.Val é uma escolha deliberada, então o eixo é rotulado como tal.

g_volc <- EnhancedVolcano(
  res,
  lab = res$rotulo, x = "logFC", y = "adj.P.Val",
  pCutoff = FDR_CUT, FCcutoff = LFC_CUT,
  ylab = bquote(~-Log[10] ~ "p-valor ajustado (BH)"),
  title = "GSE1009 — Nefropatia Diabética vs. Controle",
  subtitle = "Affymetrix HG-U95Av2 (GPL8300)",
  caption = sprintf("Corte: adj.P.Val < %.2f e |logFC| > %.2f", FDR_CUT, LFC_CUT),
  labSize = 3, pointSize = 2
)
ggsave(file.path(DIR_OUT, "05_volcano.pdf"), g_volc, width = 10, height = 9)


# =============================================================================
# 12. HEATMAP
# =============================================================================
# ATENÇÃO CONCEITUAL: os genes plotados foram ESCOLHIDOS por diferirem entre os
# grupos. O heatmap vai separar os grupos por construção. Ele ilustra o padrão;
# ele NÃO é evidência independente de separação. (Ver auditoria, Problema #7.)

if (nrow(degs) > 0) {
  probes_heat <- degs %>% arrange(adj.P.Val) %>% slice_head(n = TOP_HEAT)
  titulo_heat <- sprintf("Top %d DEGs (adj.P < %.2f, |logFC| > %.1f)",
                         nrow(probes_heat), FDR_CUT, LFC_CUT)
} else {
  message("Nenhum DEG passou nos dois cortes — heatmap exploratório com os ",
          TOP_HEAT, " menores p-valores BRUTOS (NÃO são DEGs).")
  probes_heat <- res %>% arrange(P.Value) %>% slice_head(n = TOP_HEAT)
  titulo_heat <- sprintf("Top %d sondas por p-valor BRUTO — exploratório, ",
                         nrow(probes_heat))
  titulo_heat <- paste0(titulo_heat, "sem significância após FDR")
}

mat_heat <- expr[probes_heat$PROBEID, , drop = FALSE]
rownames(mat_heat) <- make.unique(probes_heat$rotulo)
anot_col <- data.frame(Grupo = grupo, Individuo = individuo,
                       row.names = colnames(mat_heat))

pdf(file.path(DIR_OUT, "06_heatmap.pdf"), width = 9, height = 11)
pheatmap(mat_heat, scale = "row", annotation_col = anot_col,
         annotation_colors = list(Grupo = cores_grupo),
         main = titulo_heat, fontsize_row = 7)
dev.off()


# =============================================================================
# 13. DIAGRAMA DE VENN — sobreposição dos dois critérios
# =============================================================================
# Os conjuntos são de PROBE SETS (não de genes) — o rótulo reflete isso.

set_p  <- res$PROBEID[res$adj.P.Val < FDR_CUT]
set_fc <- res$PROBEID[abs(res$logFC) > LFC_CUT]

if (length(set_p) > 0 && length(set_fc) > 0) {
  lista_venn <- setNames(list(set_p, set_fc),
                         c(sprintf("adj.P.Val < %.2f", FDR_CUT),
                           sprintf("|logFC| > %.1f", LFC_CUT)))
  g_venn <- ggVennDiagram(lista_venn) +
    labs(title = "GSE1009 — probe sets por critério de corte") +
    theme(legend.position = "none")
  ggsave(file.path(DIR_OUT, "07_venn.pdf"), g_venn, width = 7, height = 7)
} else {
  message("Venn não gerado: pelo menos um dos conjuntos está vazio (",
          length(set_p), " por FDR, ", length(set_fc), " por logFC).")
}


# =============================================================================
# 14. ENRIQUECIMENTO FUNCIONAL (GO e KEGG)
# =============================================================================
# unique() é obrigatório: um gene com 5 sondas contaria 5 vezes e distorceria
# tanto a lista de interesse quanto o universo do teste hipergeométrico.
# O universo é o conjunto de genes MEDIDOS na plataforma, não o genoma inteiro —
# senão todo termo enriquecido em "genes que o HG-U95Av2 mede" aparece como hit.

genes_sig <- unique(na.omit(degs$ENTREZID))
universo  <- unique(na.omit(res$ENTREZID))
message("Genes únicos na lista: ", length(genes_sig), " | universo: ", length(universo))

salvar_enriquecimento <- function(obj, nome, titulo) {
  if (is.null(obj)) { message("Nenhum resultado para ", nome); return(invisible(NULL)) }
  df <- as.data.frame(obj)
  write.csv(df, file.path(DIR_OUT, paste0("GSE1009_", nome, ".csv")), row.names = FALSE)
  if (nrow(df) == 0) { message("Nenhum termo significativo em ", nome); return(invisible(NULL)) }
  ggsave(file.path(DIR_OUT, paste0("08_", nome, "_dotplot.pdf")),
         enrichplot::dotplot(obj, showCategory = 15) + labs(title = titulo),
         width = 10, height = 9)
}

if (length(genes_sig) >= 5) {
  ego <- enrichGO(gene = genes_sig, universe = universo, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  salvar_enriquecimento(ego, "GO_BP", "GSE1009 — GO Biological Process")

  ekegg <- enrichKEGG(gene = genes_sig, universe = universo, organism = "hsa",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  salvar_enriquecimento(ekegg, "KEGG", "GSE1009 — vias KEGG")
} else {
  message("Lista de genes muito curta (", length(genes_sig),
          ") para ORA confiável. Alternativa recomendada: GSEA sobre a estatística t ",
          "de TODOS os genes (clusterProfiler::gseGO), que não exige corte prévio.")
}


# =============================================================================
# 15. REGISTRO DE REPRODUTIBILIDADE
# =============================================================================

writeLines(c(
  paste("Data de execução:", data_execucao),
  paste("Estratégia de réplicas:", ESTRATEGIA_REPLICA),
  paste("Cortes: adj.P.Val <", FDR_CUT, "| |logFC| >", LFC_CUT),
  paste("Arrays:", ncol(expr), "| probe sets analisados:", nrow(expr)),
  "", "Delineamento:", capture.output(print(metadados)),
  "", capture.output(sessionInfo())
), file.path(DIR_OUT, "sessionInfo.txt"))

message("Pipeline concluído. Saídas em: ", normalizePath(DIR_OUT))
# =============================================================================
# FIM
# =============================================================================
