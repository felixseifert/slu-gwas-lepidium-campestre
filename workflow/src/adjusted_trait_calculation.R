library(lme4)

# trait: seed weight
seed_weight_data = read.csv("analysis/phenotype/seed_weight_processed_renamed.csv", sep="\t")
attach(seed_weight_data)
SEED_WEIGHT=as.numeric(SEED_WEIGHT)
LINE=as.factor(LINE)
ENVIRONMENT=as.factor(ENVIRONMENT)
LATE_SOWING=as.factor(LATE_SOWING)

model = lmer(SEED_WEIGHT ~ (1|LINE) + (1|LINE:ENVIRONMENT) + (1|LATE_SOWING:ENVIRONMENT))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "seed_weight")
write.table(adjusted_output, "analysis/phenotype/adjusted/seed_weight.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(seed_weight_data)
rm(list=ls())

# trait: seed weight (only first batch without late sowing)
seed_weight_data = read.csv("analysis/phenotype/seed_weight_processed_renamed.csv", sep="\t")
seed_weight_data = seed_weight_data[seed_weight_data$LATE_SOWING==0,]
attach(seed_weight_data)
SEED_WEIGHT=as.numeric(SEED_WEIGHT)
LINE=as.factor(LINE)
ENVIRONMENT=as.factor(ENVIRONMENT)

model = lmer(SEED_WEIGHT ~ (1|LINE) + (1|LINE:ENVIRONMENT))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "seed_weight")
write.table(adjusted_output, "analysis/phenotype/adjusted/seed_weight_adjusted_first_batch_subset.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(seed_weight_data)
rm(list=ls())

# trait: flowering time
flowering_time_data = read.csv("analysis/phenotype/flowering_time_processed_renamed.csv", sep="\t")
attach(flowering_time_data)
LINE=as.factor(LINE)
FLOWERING_DAP=as.numeric(FLOWERING_DAP)
ENVIRONMENT=as.factor(ENVIRONMENT)
LATE_SOWING=as.factor(LATE_SOWING)

model = lmer(FLOWERING_DAP ~ (1|LINE) + (1|ENVIRONMENT:LATE_SOWING))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "flowering_time")
write.table(adjusted_output, "analysis/phenotype/adjusted/flowering_time.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(flowering_time_data)
rm(list=ls())

# trait group: glucosinolate content traits
# trait: glucoallyssin
glucoallyssin_content_data = read.csv("analysis/phenotype/glucoallyssin_processed_renamed.csv", sep="\t")
attach(glucoallyssin_content_data)
LINE=as.factor(LINE)
GLUCOALLYSSIN_CONCENTRATION=as.numeric(GLUCOALLYSSIN)
model = lmer(GLUCOALLYSSIN_CONCENTRATION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "glucoallyssin_concentration")
write.table(adjusted_output, "analysis/phenotype/adjusted/glucoallyssin_concentration.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(glucoallyssin_content_data)
rm(list=ls())

# trait: sinalbin
sinalbin_content_data = read.csv("analysis/phenotype/sinalbin_processed_renamed.csv", sep="\t")
attach(sinalbin_content_data)
LINE=as.factor(LINE)
SINALBIN_CONCENTRATION=as.numeric(SINALBIN)
model = lmer(SINALBIN_CONCENTRATION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "sinalbin_concentration")
write.table(adjusted_output, "analysis/phenotype/adjusted/sinalbin_concentration.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(sinalbin_content_data)
rm(list=ls())

# trait: glucoallyssin + sinalbin (total glucosinolate)
glucosinolate_content_data = read.csv("analysis/phenotype/glucosinolate_processed_renamed.csv", sep="\t")
attach(glucosinolate_content_data)
LINE=as.factor(LINE)
GLUCOSINOLATE=as.numeric(GLUCOSINOLATE)
model = lmer(GLUCOSINOLATE ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "glucosinolate_concentration")
write.table(adjusted_output, "analysis/phenotype/adjusted/glucosinolate_concentration.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(glucosinolate_content_data)
rm(list=ls())

# trait group: oil content/composition traits
# trait: oil content
oilcontent_oilcomposition_data = read.csv("analysis/phenotype/oil_content_processed_renamed.csv", sep="\t")
attach(oilcontent_oilcomposition_data)
LINE=as.factor(LINE)
OIL_CONTENT=as.numeric(OIL_CONTENT)
ENVIRONMENT=as.factor(ENVIRONMENT)
REPLICATE=as.factor(REPLICATE)
model = lmer(OIL_CONTENT ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "oil_content")
write.table(adjusted_output, "analysis/phenotype/adjusted/oil_content.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(oilcontent_oilcomposition_data)
rm(list=ls())

# trait: oleic acid fraction
oleic_acid_data = read.csv("analysis/phenotype/oleic_acid_processed_renamed.csv", sep="\t")
attach(oleic_acid_data)
LINE=as.factor(LINE)
OLEIC_ACID_FRACTION=as.numeric(OLEIC_ACID)
ENVIRONMENT=as.factor(ENVIRONMENT)
REPLICATE=as.factor(REPLICATE)
model = lmer(OLEIC_ACID_FRACTION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "oleic_acid_fraction")
write.table(adjusted_output, "analysis/phenotype/adjusted/oleic_acid_fraction.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(oleic_acid_data)
rm(list=ls())

# trait: erucic acid fraction
erucic_acid_data = read.csv("analysis/phenotype/erucic_acid_processed_renamed.csv", sep="\t")
attach(erucic_acid_data)
LINE=as.factor(LINE)
ERUCIC_ACID_FRACTION=as.numeric(ERUCIC_ACID)
ENVIRONMENT=as.factor(ENVIRONMENT)
REPLICATE=as.factor(REPLICATE)
model = lmer(ERUCIC_ACID_FRACTION ~ (1|LINE) + (1|ENVIRONMENT:REPLICATE))
adjusted = coef(model)$LINE
adjusted_output = data.frame(rep(0, length(adjusted[,1])), rownames(adjusted), adjusted[,1])
names(adjusted_output) = c("FID", "IID", "erucic_acid_fraction")
write.table(adjusted_output, "analysis/phenotype/adjusted/erucic_acid_fraction.adjusted.tsv", col.names=T, quote=F, row.names=F, sep="\t")

detach(erucic_acid_data)
rm(list=ls())
