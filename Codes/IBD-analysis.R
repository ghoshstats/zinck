library(ggplot2)
library(rstan)
library(knockoff)
library(randomForest)
library(phyloseq)
library(dplyr)
library(zinck)
library(reshape2)
library(topicmodels)
library(VennDiagram)
library(grid)
library(RColorBrewer)

################################################################################
###################### IBD data genus level ####################################

load("genera.RData") ## Loading the meta IBD studies
combined_studies <- as.data.frame(t(physeq_genera@otu_table))
study_names <- physeq_genera@sam_data[["dataset_name"]]

## RISK ##
risk_indices <- which(study_names == "RISK")
risk_otu <- combined_studies[risk_indices, ]

## CS-PRISM ##
prism_indices <- which(study_names == "CS-PRISM")
prism_otu <- combined_studies[prism_indices, ]

## HMP2 ##
HMP_indices <- which(study_names == "HMP2")
hmp_otu <- combined_studies[HMP_indices, ]

## Pouchitis ##
pouchitis_indices <- which(study_names == "Pouchitis")
pouchitis_otu <- combined_studies[pouchitis_indices, ]

## MucosalIBD ##
mi_indices <- which(study_names == "MucosalIBD")
mi_otu <- combined_studies[mi_indices, ]

metadata <- physeq_genera@sam_data

# Get the indices of the last occurrence of each unique subject_accession
latest_indices <- sapply(unique(metadata$subject_accession), function(x) {
  max(which(metadata$subject_accession == x))
})

# Subset the metadata to keep only the latest samples
latest_metadata <- metadata[latest_indices, ]

study_names <- c("CS-PRISM","HMP2","MucosalIBD","Pouchitis","RISK")
meta_IBD <- latest_metadata[latest_metadata$dataset_name %in% study_names,]

# Extract row names for each study
risk_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "RISK", ])
prism_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "CS-PRISM", ])
hmp_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "HMP2", ])
mi_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "MucosalIBD", ])
pouchitis_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "Pouchitis", ])

# Subset the OTU matrices and Y vectors based on these row names
risk_otu_latest <- risk_otu[row.names(risk_otu) %in% risk_row_names, ]
prism_otu_latest <- prism_otu[row.names(prism_otu) %in% prism_row_names, ]
hmp_otu_latest <- hmp_otu[row.names(hmp_otu) %in% hmp_row_names, ]
mi_otu_latest <- mi_otu[row.names(mi_otu) %in% mi_row_names, ]
pouchitis_otu_latest <- pouchitis_otu[row.names(pouchitis_otu) %in% pouchitis_row_names, ]

combined_otu <- rbind(prism_otu_latest,hmp_otu_latest,mi_otu_latest,
                      pouchitis_otu_latest,risk_otu_latest)
combined_otu <- combined_otu[ match(rownames(meta_IBD), rownames(combined_otu)), ] 
## To make sure samples in OTU correspond to the samples in meta

IBD_resp <- meta_IBD$disease
Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

X <- combined_otu

dlt <- rep(0,ncol(X)) ## Initializing the deltas with the sparsity of each column

for(t in (1:ncol(X)))
{
  dlt[t] <- 1-mean(X[,t]>0)
  if(dlt[t]==0)
  {
    dlt[t] = dlt[t]+0.01
  }
  if (dlt[t]==1)
  {
    dlt[t] = dlt[t]-0.01
  }
  
}

zinLDA_stan_data <- list(
  K = 27,
  V = ncol(X),
  D = nrow(X),
  n = X,
  delta = dlt
)

## Fitting the zinck model ##
set.seed(1)
fitIBD <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", importance_resampling=TRUE, iter=10000,tol_rel_obj=0.01,elbo_samples=500)

theta <- fitIBD@sim[["est"]][["theta"]]
beta <- fitIBD@sim[["est"]][["beta"]]
X_tilde <- zinck::generateKnockoff(X,theta,beta,seed=1) ## Generating the kncokoff copy

X_aug <- cbind(X,X_tilde) ## Creating the augmented matrix

######### Tuning the Random Forest model ####################
bestmtry <- tuneRF(X_aug,as.factor(Y),stepFactor = 1.5,improve=1e-5,ntree=1000) ## Getting the best mtry hyperparameter
m <- bestmtry[as.numeric(which.min(bestmtry[,"OOBError"])),1]
# m = 109
df_X <- as.data.frame(W_aug)
colnames(df_X) <- NULL
rownames(df_X) <- NULL
df_X$Y <- Y
set.seed(1)
model_rf <- randomForest(formula=as.factor(Y)~.,ntree=1000,mtry=m,importance=TRUE,data=as.matrix(df_X)) ## Fitting the tuned Random Forest model
cf <- as.data.frame(importance(model_rf))[,3] ## Extracting the Mean Decrease in Impurities for each variable
W <- abs(cf[1:249])-abs(cf[250:498])
T <- knockoff.threshold(W,fdr = 0.1, offset = 0) ## This is the knockoff threshold
T <- zinck.threshold(W,fdr=0.1)
names_zinck <- colnames(X[,which(W>=T)])

# 2, 30,35, 36, 39, 55,58,69,77,79,90, 94 selected

############## Feature Importance Plot #############

names_zinck <- c("Bilophila genus (f. Desulfovibrionaceae)","Unknown genus (f. Ruminococcaceae)","Holdemania genus (f. Erysipelotrichaceae)","Parabacteroides genus (f. Porphyromonadaceae)","Unknown genus (f. Rikenellaceae)", "Bacteroides genus (f. Bacteroidaceae)","Oscillospira genus (f. Ruminococcaceae)","Phascolarctobacterium genus (f. Veillonellaceae)","Coprococcus genus (f. Lachnospiraceae)","Blautia genus (f. Lachnospiraceae)", "Lachnospira genus (f. Lachnospiraceae)","Faecalibacterium genus (f. Ruminococcaceae)")
Xnorm <- X / rowSums(X)

# Calculate the column sums for cases and controls
case_sums <- colMeans(Xnorm[Y == 1, which(W>=T)])
control_sums <- colMeans(Xnorm[Y == 0, which(W>=T)])

colors <- ifelse(case_sums > control_sums, "red", "blue")


# Create a vector to indicate which labels should be bold
bold_indices <- c(3,10)
bold_labels <- ifelse(seq_along(names_zinck) %in% bold_indices, "bold", "plain")

# Create a data frame for plotting
data.species <- data.frame(
  impscores = sort(W[which(W >= T)], decreasing = FALSE),
  name = factor(names_zinck, levels = names_zinck),
  y = seq(length(names_zinck)) * 0.9,
  color = colors,
  fontface = bold_labels
)

tiff("Zinck_IBD_feature_scores_notlognormalized.tiff", units="in", width=10, height=6, res=600)

# Create the plot
plt.species <- ggplot(data.species) +
  geom_col(aes(impscores, name, fill = color), width = 0.6) +
  scale_fill_identity() +
  theme_bw() +
  ylab("Genera") +
  xlab("Feature Statistic") +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, face = data.species$fontface) # Adjust fontface based on the condition
  )
print(plt.species)
dev.off()

############################ Plotting the Venn Diagram ################################
#######################################################################################
                      
zinCK <- c(2, 30,35, 36, 39, 55,58,69,77,79,90, 94)
DL <- c(12,17,18,30,36,44,49,52,70,88,90,94,96)
KF <- c(2, 10, 12, 39, 49, 59, 64, 72, 88, 90)
vanilla_LDA <- c(2, 30,36, 39,55,58,69,77, 90,94)
venn.plot <- venn.diagram(
  x = list(zinCK = zinCK, KF = KF, DeepLINK = DL, vanilla_LDA = vanilla_LDA),
  category.names = c("Zinck" , "MX-KF" , "DeepLINK", "LDA-KF"),
  filename = 'IBD_venn_diagram_not_lognormalized.png',
  output = TRUE,  # Set to FALSE so that we can draw it with grid.draw
  imagetype="png",
  height = 600,
  width = 600,
  resolution = 400,
  lwd = 1,
  col=c("red", 'blue', 'green', "orange"),
  fill = c(alpha("red", 0.3), alpha('blue', 0.3), alpha('green', 0.3), alpha('orange', 0.3)),
  cex = 0.7,
  fontfamily = "sans",
  cat.cex = 0,  # Modified for readability
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0, 0),  # Modified for 4 categories
  cat.dist = c(0.055, 0.055, 0.055, 0.055),  # Modified for 4 categories
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', 'green', 'orange')  # Modified for 4 categories
)




