# SLU, GWAS Lepidium
# Author: Felix Seifert, cropSeq bioinformatics

# phenotype data adjustment
# paths are hardcoded for analysis

library(lme4)

filter_outliers <- function(data, trait) {
    quartiles <- quantile(unlist(data[trait]), probs=c(.25, .75), na.rm = TRUE)
    inter_quartile_range <- IQR(unlist(data[trait]))
    lower_border <- quartiles[1] - 1.5 * inter_quartile_range
    upper_border <- quartiles[2] + 1.5 * inter_quartile_range
    data <- subset(data, data[trait] > lower_border & data[trait] < upper_border)
    return(data)
}

# trait: seed weight
trait <- "SEED_WEIGHT"
data <- read.csv("analysis/phenotype/seed_weight_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
SEED_WEIGHT <- as.numeric(SEED_WEIGHT)
LINE <- as.factor(LINE)
ENVIRONMENT <- as.factor(ENVIRONMENT)
LATE_SOWING <- as.factor(LATE_SOWING)

model <- lmer(SEED_WEIGHT ~ (1|LINE) + (1|LINE:ENVIRONMENT) + (1|LATE_SOWING:ENVIRONMENT))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "seed_weight")
write.table(adjusted_output, "analysis/phenotype/adjusted/seed_weight.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait: seed weight (only first batch without late sowing)
trait <- "SEED_WEIGHT"
data <- read.csv("analysis/phenotype/seed_weight_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
SEED_WEIGHT <- as.numeric(SEED_WEIGHT)
LINE <- as.factor(LINE)
ENVIRONMENT <- as.factor(ENVIRONMENT)

model <- lmer(SEED_WEIGHT ~ (1|LINE) + (1|LINE:ENVIRONMENT))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "seed_weight")
write.table(adjusted_output, "analysis/phenotype/adjusted/seed_weight_adjusted_first_batch_subset.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait: flowering time
trait <- "FLOWERING_DAP"
data <- read.csv("analysis/phenotype/flowering_time_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
LINE <- as.factor(LINE)
FLOWERING_DAP <- as.numeric(FLOWERING_DAP)
ENVIRONMENT <- as.factor(ENVIRONMENT)
LATE_SOWING <- as.factor(LATE_SOWING)

model <- lmer(FLOWERING_DAP ~ (1|LINE) + (1|ENVIRONMENT:LATE_SOWING))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "flowering_time")
write.table(adjusted_output, "analysis/phenotype/adjusted/flowering_time.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait group: glucosinolate content traits
# trait: glucoallyssin
trait <- "GLUCOALLYSSIN"
data <- read.csv("analysis/phenotype/glucoallyssin_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
LINE <- as.factor(LINE)
GLUCOALLYSSIN_CONCENTRATION <- as.numeric(GLUCOALLYSSIN)

model <- lmer(GLUCOALLYSSIN_CONCENTRATION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "glucoallyssin_concentration")
write.table(adjusted_output, "analysis/phenotype/adjusted/glucoallyssin_concentration.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait: sinalbin
trait <- "SINALBIN"
data <- read.csv("analysis/phenotype/sinalbin_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
LINE <- as.factor(LINE)
SINALBIN_CONCENTRATION <- as.numeric(SINALBIN)

model <- lmer(SINALBIN_CONCENTRATION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "sinalbin_concentration")
write.table(adjusted_output, "analysis/phenotype/adjusted/sinalbin_concentration.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait: glucoallyssin + sinalbin (total glucosinolate)
trait <- "GLUCOSINOLATE"
data <- read.csv("analysis/phenotype/glucosinolate_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
LINE <- as.factor(LINE)
GLUCOSINOLATE <- as.numeric(GLUCOSINOLATE)

model <- lmer(GLUCOSINOLATE ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "glucosinolate_concentration")
write.table(adjusted_output, "analysis/phenotype/adjusted/glucosinolate_concentration.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait group: oil content/composition traits
# trait: oil content
trait <- "OIL_CONTENT"
data <- read.csv("analysis/phenotype/oil_content_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
LINE <- as.factor(LINE)
OIL_CONTENT <- as.numeric(OIL_CONTENT)
ENVIRONMENT <- as.factor(ENVIRONMENT)
REPLICATE <- as.factor(REPLICATE)

model <- lmer(OIL_CONTENT ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "oil_content")
write.table(adjusted_output, "analysis/phenotype/adjusted/oil_content.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait: oleic acid fraction
trait <- "OLEIC_ACID"
data <- read.csv("analysis/phenotype/oleic_acid_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
LINE <- as.factor(LINE)
OLEIC_ACID_FRACTION <- as.numeric(OLEIC_ACID)
ENVIRONMENT <- as.factor(ENVIRONMENT)
REPLICATE <- as.factor(REPLICATE)

model <- lmer(OLEIC_ACID_FRACTION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "oleic_acid_fraction")
write.table(adjusted_output, "analysis/phenotype/adjusted/oleic_acid_fraction.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))

# trait: erucic acid fraction
trait <- "ERUCIC_ACID"
data <- read.csv("analysis/phenotype/erucic_acid_processed_renamed.csv", sep="\t")
print(trait)
print(paste("unfiltered: ", nrow(data)))
data <- filter_outliers(data, trait)
print(paste("outlier-filtered: ", nrow(data)))
print("")
attach(data)
LINE <- as.factor(LINE)
ERUCIC_ACID_FRACTION <- as.numeric(ERUCIC_ACID)
ENVIRONMENT <- as.factor(ENVIRONMENT)
REPLICATE <- as.factor(REPLICATE)

model <- lmer(ERUCIC_ACID_FRACTION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted <- coef(model)$LINE
adjusted_output <- data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) <- c("FID", "IID", "erucic_acid_fraction")
write.table(adjusted_output, "analysis/phenotype/adjusted/erucic_acid_fraction.adjusted_filtered.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(data)
rm(list=setdiff(ls(), "filter_outliers"))
