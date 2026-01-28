# 플라스미드 한개로 테스트: t-test 분석
# ==========================================================
# 플라스미드 전사체 분석 (data source: Pistek et al., 2023)
# ==========================================================
setwd("D:/rnaseq_raw_SRP424429/RNA-seq-test-nextflow-results")
list.files()
list.files("results_27samples/salmon")
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(writexl) 

# ========================================
# RepCap DEG 분석: 전체 워크플로우
# ========================================

# ========================================
# 1. Salmon 결과 불러오기
# ========================================

cat("=== 1단계: Salmon 결과 불러오기 ===\n\n")

# Salmon 결과가 있는 샘플 폴더 확인
salmon_dir <- "results_27samples/salmon"
actual_samples <- list.dirs(salmon_dir, recursive = FALSE, full.names = FALSE)

cat("Salmon 결과 샘플:", length(actual_samples), "개\n")
print(actual_samples)

# samplesheet에서 전체 샘플 확인
samples <- read.csv("samplesheet_27.csv", stringsAsFactors = FALSE)
all_sample_names <- samples$sample

# Mock 샘플 (Salmon 결과 없음)
mock_samples <- setdiff(all_sample_names, actual_samples)

cat("\nMock 샘플:", length(mock_samples), "개\n")
print(mock_samples)

# ========================================
# 2. tximport로 Salmon 데이터 불러오기
# ========================================

cat("\n=== 2단계: tximport 실행 ===\n\n")

# Salmon quant.sf 파일 경로
files <- file.path(salmon_dir, actual_samples, "quant.sf")
names(files) <- actual_samples

# 파일 존재 확인
cat("모든 파일 존재:", all(file.exists(files)), "\n")

if(!all(file.exists(files))) {
  cat("존재하지 않는 파일:\n")
  print(files[!file.exists(files)])
  stop("일부 파일이 존재하지 않습니다!")
}

# tximport 실행
txi <- tximport(
  files,
  type = "salmon",
  txOut = TRUE  # 전사체 수준
)

cat("tximport 완료\n")
cat("데이터 크기:", dim(txi$counts), "\n")

# ========================================
# 3. Mock 샘플 추가 (count = 0)
# ========================================

cat("\n=== 3단계: Mock 샘플 추가 ===\n\n")

# Mock 샘플에 대한 0 count 행렬 생성
n_transcripts <- nrow(txi$counts)
n_mock <- length(mock_samples)

mock_counts <- matrix(0, nrow = n_transcripts, ncol = n_mock)
colnames(mock_counts) <- mock_samples
rownames(mock_counts) <- rownames(txi$counts)

mock_abundance <- matrix(0, nrow = n_transcripts, ncol = n_mock)
colnames(mock_abundance) <- mock_samples
rownames(mock_abundance) <- rownames(txi$abundance)

mock_length <- matrix(mean(txi$length), nrow = n_transcripts, ncol = n_mock)
colnames(mock_length) <- mock_samples
rownames(mock_length) <- rownames(txi$length)

# AAV와 Mock 통합
txi_complete <- list(
  counts = cbind(txi$counts, mock_counts),
  abundance = cbind(txi$abundance, mock_abundance),
  length = cbind(txi$length, mock_length),
  countsFromAbundance = "no"
)

cat("통합 데이터 크기:", dim(txi_complete$counts), "\n")
cat("총 샘플 수:", ncol(txi_complete$counts), "\n")

# ========================================
# 4. 메타데이터 생성
# ========================================

cat("\n=== 4단계: 메타데이터 생성 ===\n\n")

all_samples <- c(actual_samples, mock_samples)

metadata <- data.frame(
  sample = all_samples,
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Mock 여부 판단
    is_mock = sample %in% mock_samples,
    
    # 샘플명 파싱
    # 예: cell_line_1_0h_rep1 -> cell_line_1, 0h, rep1
    cell_line = str_extract(sample, "cell_line_\\d+"),
    timepoint = str_extract(sample, "\\d+h"),
    replicate = str_extract(sample, "rep\\d+"),
    
    # Condition
    condition = ifelse(is_mock, "Mock", "AAV"),
    
    # 시간 숫자로 변환
    time_numeric = as.numeric(gsub("h", "", timepoint))
  ) %>%
  # txi_complete의 컬럼 순서와 맞추기
  arrange(match(sample, colnames(txi_complete$counts)))

cat("메타데이터 생성 완료:", nrow(metadata), "개 샘플\n\n")

# 확인
cat("Cell lines:", paste(unique(metadata$cell_line), collapse = ", "), "\n")
cat("Timepoints:", paste(unique(metadata$timepoint), collapse = ", "), "\n")
cat("Conditions:", paste(unique(metadata$condition), collapse = ", "), "\n\n")

# 샘플 순서 일치 확인
if(!all(metadata$sample == colnames(txi_complete$counts))) {
  stop("메타데이터와 count 데이터의 샘플 순서가 일치하지 않습니다!")
}

cat("✓ 메타데이터와 count 데이터 순서 일치 확인\n")

# 메타데이터 저장
write.csv(metadata, "metadata_complete.csv", row.names = FALSE)

# ========================================
# 5. Count 데이터프레임 생성
# ========================================

cat("\n=== 5단계: Count 데이터프레임 생성 ===\n\n")

counts_df <- as.data.frame(t(txi_complete$counts)) %>%
  rownames_to_column("sample") %>%
  left_join(metadata, by = "sample")

cat("counts_df 생성 완료\n")
cat("크기:", dim(counts_df), "\n\n")

# 데이터 확인
cat("=== 데이터 미리보기 ===\n")
print(head(counts_df, 10))

cat("\n=== 조건별 샘플 수 ===\n")
print(table(counts_df$condition))
print(table(counts_df$cell_line, counts_df$timepoint))

# 저장
write.csv(counts_df, "RepCap_counts_complete.csv", row.names = FALSE)


# ========================================
# 올바른 DEG 분석: 0h를 baseline으로 24h, 72h 비교
# ========================================

# counts_df 확인
print("=== 데이터 확인 ===")
print(table(counts_df$cell_line, counts_df$timepoint))

# ========================================
# DEG 분석: 각 cell line에서 0h vs 24h, 0h vs 72h
# ========================================

deg_table <- data.frame()

for(cl in unique(na.omit(counts_df$cell_line))) {
  
  cat(sprintf("\n=== %s 분석 ===\n", cl))
  
  # 0h 샘플 (baseline)
  baseline_0h <- counts_df %>% 
    filter(cell_line == cl, timepoint == "0h")
  
  cat(sprintf("0h 샘플 수: %d\n", nrow(baseline_0h)))
  
  # 24h와 72h 각각 비교
  for(tp in c("24h", "72h")) {
    
    cat(sprintf("분석 중: %s, 0h vs %s\n", cl, tp))
    
    test_samples <- counts_df %>% 
      filter(cell_line == cl, timepoint == tp)
    
    cat(sprintf("%s 샘플 수: %d\n", tp, nrow(test_samples)))
    
    if(nrow(baseline_0h) >= 2 && nrow(test_samples) >= 2) {
      
      # 평균, 표준편차
      baseline_mean <- mean(baseline_0h$RepCap, na.rm = TRUE)
      baseline_sd <- sd(baseline_0h$RepCap, na.rm = TRUE)
      if(is.na(baseline_sd)) baseline_sd <- 0
      
      test_mean <- mean(test_samples$RepCap, na.rm = TRUE)
      test_sd <- sd(test_samples$RepCap, na.rm = TRUE)
      if(is.na(test_sd)) test_sd <- 0
      
      # Fold Change
      if(baseline_mean == 0) {
        fold_change <- Inf
        log2fc <- Inf
      } else {
        fold_change <- test_mean / baseline_mean
        log2fc <- log2(test_mean / baseline_mean)
      }
      
      # t-test
      test_result <- t.test(test_samples$RepCap, baseline_0h$RepCap)
      
      # TPM
      baseline_sample_names <- baseline_0h$sample
      test_sample_names <- test_samples$sample
      
      baseline_tpm <- mean(txi_complete$abundance["RepCap", baseline_sample_names], na.rm = TRUE)
      test_tpm <- mean(txi_complete$abundance["RepCap", test_sample_names], na.rm = TRUE)
      
      # 결과 저장
      result <- data.frame(
        Gene = "RepCap",
        Cell_Line = cl,
        Comparison = paste0("0h_vs_", tp),
        Baseline = "0h",
        Test_Timepoint = tp,
        
        Baseline_N = nrow(baseline_0h),
        Baseline_Mean_Count = baseline_mean,
        Baseline_SD_Count = baseline_sd,
        Baseline_TPM = baseline_tpm,
        
        Test_N = nrow(test_samples),
        Test_Mean_Count = test_mean,
        Test_SD_Count = test_sd,
        Test_TPM = test_tpm,
        
        Fold_Change = fold_change,
        Log2_Fold_Change = log2fc,
        
        P_Value = test_result$p.value,
        T_Statistic = as.numeric(test_result$statistic),
        DF = as.numeric(test_result$parameter),
        CI_Lower = test_result$conf.int[1],
        CI_Upper = test_result$conf.int[2],
        
        stringsAsFactors = FALSE
      )
      
      deg_table <- rbind(deg_table, result)
      cat(sprintf("  ✓ Fold Change: %.2f, P-value: %.4f\n", fold_change, test_result$p.value))
      
    } else {
      cat(sprintf("  경고: 샘플 수 부족 (0h: %d, %s: %d)\n", 
                  nrow(baseline_0h), tp, nrow(test_samples)))
    }
  }
}

# ========================================
# FDR (Benjamini-Hochberg) 보정 추가
# ========================================

deg_table <- deg_table %>%
  mutate(
    # FDR 보정 (Benjamini-Hochberg method)
    FDR = p.adjust(P_Value, method = "BH"),
    
    # 보정된 p-value 기준 유의성
    Significant_FDR_0.05 = FDR < 0.05,
    Significant_FDR_0.01 = FDR < 0.01,
    Significant_FDR_0.001 = FDR < 0.001,
    
    # 원래 p-value 기준 유의성 (참고용)
    Significant_P_0.05 = P_Value < 0.05,
    Significant_P_0.01 = P_Value < 0.01
  )

cat("\n✓ DEG 분석 완료 (FDR 보정 포함):", nrow(deg_table), "개 비교\n")
print(deg_table %>% select(Cell_Line, Comparison, Fold_Change, P_Value, FDR))

# ========================================
# Publication Table (FDR 포함)
# ========================================

publication_table <- deg_table %>%
  mutate(
    # 0h (Baseline)
    Baseline_Count = sprintf("%.1f ± %.1f", Baseline_Mean_Count, Baseline_SD_Count),
    
    # Test timepoint
    Test_Count = sprintf("%.1f ± %.1f", Test_Mean_Count, Test_SD_Count),
    
    # TPM
    Baseline_TPM_Display = sprintf("%.2f", Baseline_TPM),
    Test_TPM_Display = sprintf("%.2f", Test_TPM),
    
    # Fold Change
    FC_Display = ifelse(is.infinite(Fold_Change), 
                        "Inf", 
                        sprintf("%.2f", Fold_Change)),
    
    # Log2 Fold Change
    Log2FC_Display = ifelse(is.infinite(Log2_Fold_Change),
                            "Inf",
                            sprintf("%.2f", Log2_Fold_Change)),
    
    # P-value
    P_Value_Display = case_when(
      P_Value < 0.0001 ~ "< 0.0001",
      TRUE ~ sprintf("%.4f", P_Value)
    ),
    
    # FDR
    FDR_Display = case_when(
      FDR < 0.0001 ~ "< 0.0001",
      TRUE ~ sprintf("%.4f", FDR)
    ),
    
    # Significance (FDR 기준)
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(Cell_Line, Comparison,
         Baseline_Count, Test_Count,
         Baseline_TPM_Display, Test_TPM_Display,
         FC_Display, Log2FC_Display,
         P_Value_Display, FDR_Display, Significance) %>%
  rename(
    `Cell Line` = Cell_Line,
    `Comparison` = Comparison,
    `0h Count (Mean ± SD)` = Baseline_Count,
    `Test Count (Mean ± SD)` = Test_Count,
    `0h TPM` = Baseline_TPM_Display,
    `Test TPM` = Test_TPM_Display,
    `Fold Change` = FC_Display,
    `Log2 FC` = Log2FC_Display,
    `P-value` = P_Value_Display,
    `FDR (adj. P-value)` = FDR_Display,
    `Sig.` = Significance
  )

# ========================================
# Excel 저장
# ========================================

excel_output <- list(
  "DEG_Results_Detailed" = deg_table,
  "Publication_Table" = publication_table,
  "Individual_Samples" = counts_df %>% 
    select(sample, cell_line, timepoint, replicate, RepCap) %>%
    arrange(cell_line, timepoint, replicate)
)

# 바탕화면에 저장
desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "RepCap_DEG_with_FDR.xlsx")
write_xlsx(excel_output, desktop_path)

cat("\n")
cat(rep("=", 70), "\n", sep="")
cat("✓ DEG 분석 완료 (FDR 보정 포함)!\n")
cat(rep("=", 70), "\n", sep="")
cat("\n파일 저장:", desktop_path, "\n\n")

cat("=== 결과 요약 ===\n")
cat(sprintf("총 비교: %d개\n", nrow(deg_table)))
cat(sprintf("유의한 차등 발현 (FDR < 0.05): %d개\n", 
            sum(deg_table$Significant_FDR_0.05, na.rm = TRUE)))
cat(sprintf("유의한 차등 발현 (P < 0.05, 보정 전): %d개\n", 
            sum(deg_table$Significant_P_0.05, na.rm = TRUE)))

cat("\n=== Publication Table ===\n")
print(publication_table)

cat("\n=== P-value vs FDR 비교 ===\n")
print(deg_table %>% 
        select(Cell_Line, Comparison, P_Value, FDR, Significant_P_0.05, Significant_FDR_0.05))

# 파일 자동 열기
shell.exec(desktop_path)
