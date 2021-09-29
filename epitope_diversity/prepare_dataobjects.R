library("tidyverse")
all_tsr <- read_csv("epitope_diversity/data/All_TSR_mTECs.csv")

hi <- all_tsr %>%
    select(Chr, Start, Stop, contains("hi_TPM"))
hi$sum_TPM <- rowSums(select(hi, contains("hi_TPM")))

lo <- all_tsr %>%
    select(Chr, Start, Stop, contains("lo_TPM"))
lo$sum_TPM <- rowSums(select(lo, contains("lo_TPM")))

write_delim(select(hi, Chr, Start, Stop, sum_TPM),
          "epitope_diversity/data/tsr_tpm_hi_sum.bedGraph",
          col_names = FALSE,
          delim="\t")

write_delim(select(lo, Chr, Start, Stop, sum_TPM),
          "epitope_diversity/data/tsr_tpm_lo_sum.bedGraph",
          col_names = FALSE,
          delim="\t")
