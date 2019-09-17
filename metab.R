### analysis for metabolomics data, Nat Com paper

# This script analyses the data from the metabolomics data set, 
# from experiments with different conditions (Bactopepton, antibiotics,
# UV, Heat...) that affect lifespan in C. elegans


# libraries

library(tidyverse)
library(readxl)
library(ggrepel)
library(here)
library(caret)
library(Rtsne)
library(FactoMineR) 
library(factoextra) 
library(openxlsx)
library(broom)
library(ComplexHeatmap)
library(circlize)

# my library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')

# session options
options(width = 220)

metadata = read_xlsx(here::here('data', 'Metabolites.xlsx'), sheet = 'Metadata')
met = read_xlsx(here::here('data', 'Metabolites.xlsx'), sheet = 'Putative_metabolites')
qmet = read_xlsx(here::here('data', 'Metabolites.xlsx'), sheet = 'Quantitative_estimation')

# transform data
met = met %>%
	gather(Sample, Score, CT1:P4) %>%
	left_join(metadata) %>%
	mutate(ID = as.factor(ID),
		   Group = as.factor(Group),
		   Replicate = str_match_all(Sample,'[:digit:]{1}'), # gets numbers from Sample name
		   Replicate = as.numeric(Replicate),
		   Sample = as.factor(Sample)) %>%
	select(ID, Metabolite, Pathway_label, PubChem, HMDB, Sample, Group, Replicate, m_z, MT_RT, Score)

# some loops to clean the rubish 
for (i in 1:dim(met)[1]) { met$Metabolite[i] = strsplit(met$Metabolite[i], "\\r")[[1]][1] }
for (i in 1:dim(met)[1]) { met$Pathway_label[i] = strsplit(met$Pathway_label[i], "\\r")[[1]][1] }
for (i in 1:dim(met)[1]) { met$PubChem[i] = strsplit(met$PubChem[i], "\\r")[[1]][1] }
for (i in 1:dim(met)[1]) { met$HMDB[i] = strsplit(met$HMDB[i], "\\r")[[1]][1] }

# summarise data
met.sum = met %>%
	group_by(ID, Group, Metabolite, PubChem, HMDB) %>%
	summarise(Mean = mean(Score, na.rm = TRUE),
			  SD = sd(Score, na.rm = TRUE))

#########
### dataset for the paper
#########

final_met = met %>%
	mutate(Score = replace_na(Score, 2E-52)) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F) 

rownames(final_met) = final_met[,1]; final_met = final_met[,-1]

final_met = t(final_met)

write.xlsx(final_met, here('data', 'Raw_metabolomics.xlsx') , sheetName = c('Raw_metabolomics'), colNames = T, rowNames = T) 

cosa = t(scale(t(final_met), center = TRUE, scale = TRUE))




### Exploration
# check that some metabolites does not have reads/scores 
met %>% filter(Group == 'Carbenicillin', Metabolite == 'Pyruvic acid')


# lets filter our dataset
# this will count the non-NA values, so we can filter those that are == 1
removals = met %>%
	group_by(Metabolite) %>%
	summarise(N = sum(!is.na(Score))) %>%
	filter(N == 1) %>% select(Metabolite) %>% t %>% as.vector



# specify colors
cols = colsel(7, palette = 'sat1')

met.sum %>% 
	filter(Metabolite == 'Choline') %>%
	ggplot(aes(x = Group, y = Mean, colour = Group, fill = Group)) +
	geom_bar(stat = 'identity', color = "black", size = 0.5) +
	geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), colour = 'black', width = 0.2) +
	scale_fill_manual(values = cols) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))



# thousands of barplots
met.sum %>% 
ggplot(aes(x = Group, y = Mean, colour = Group, fill = Group)) +
	geom_bar(stat = 'identity', color = "black", size = 0.5) +
	geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), colour = 'black', width = 0.2) +
	scale_fill_manual(values = cols) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	facet_wrap(~ Metabolite, scales = "free")


ggsave(file = here('exploration', "Barplot_complete_metabolites.pdf"), 
	width = 400, height = 600, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")


### how many complete cases do we have in general

# title aesthetics
title = theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"))


# number of NAs in mean values
met.sum %>% 
	group_by(Group) %>% 
	summarise(sum(is.na(Mean))) %>%
	rename(NAs = `sum(is.na(Mean))`) %>%
	mutate(NAs = round((NAs / 228) * 100)) %>%
	ggplot(aes(x = Group, y = NAs, fill = Group)) +
	geom_bar(stat = 'identity', colour = 'black') +
	# geom_text(aes(label = NAs), vjust=1.6, color="black", size=5.5) +
	scale_fill_manual(values = colsel(7, palette = 'sat1')) +
	labs(title = 'Number of NAs in Mean values',
         x = 'Condition',
         y = 'Missing values (in %)') +
	theme_classic() + title +
	scale_y_continuous(limits = c(0, 80)) +
	theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12, angle = 45, hjust = 1),
  	    axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12))

quartz.save(type = 'pdf', 
	file = here('exploration', 'NAs_mean.pdf'), 
	width = 8, height = 8, family = 'Arial')

# number of NAs in sd values
met.sum %>% 
	group_by(Group) %>% 
	summarise(sum(is.na(SD))) %>%
	rename(NAs = `sum(is.na(SD))`) %>%
	mutate(NAs = (NAs / 228) * 100) %>%
	ggplot(aes(x = Group, y = NAs, fill = Group)) +
	geom_bar(stat = 'identity', colour = 'black') +
	scale_fill_manual(values = colsel(7, palette = 'sat1')) +
	labs(title = 'Number of NAs in SD values',
         x = 'Condition',
         y = 'Proportion of NAs') +
	theme_classic() + title +
	theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12, angle = 45, hjust = 1),
  	    axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12))

quartz.save(type = 'pdf', 
	file = here('exploration', 'NAs_sd.pdf'), 
	width = 8, height = 8, family = 'Arial')




### complete cases
# rearrange a bit the data
comp = met.sum %>%
	select(-SD) %>%
	spread(Group, Mean) %>%
	ungroup %>%
	select(Metabolite, Carbenicillin:Uvirradiation)

# different combinations of case withdrawal
Complete = sum(complete.cases(comp))
No_Heat = sum(complete.cases(comp %>% select(-HeatKill)))
No_Trime = sum(complete.cases(comp %>% select(-Trimethoprim)))
No_UV = sum(complete.cases(comp %>% select(-Uvirradiation)))
No_Heat_Trime = sum(complete.cases(comp %>% select(-HeatKill, -Trimethoprim)))
No_Heat_UV = sum(complete.cases(comp %>% select(-HeatKill, -Uvirradiation)))
No_ALL = sum(complete.cases(comp %>% select(-HeatKill, -Uvirradiation, -Trimethoprim)))

nums = c(Complete, No_Heat, No_Trime, No_UV, No_Heat_Trime, No_Heat_UV, No_ALL)
Conditions = c('Complete', 'No_Heat', 'No_Trime', 'No_UV', 'No_Heat_Trime', 'No_Heat_UV', 'No_ALL')

df = data.frame(nums, Conditions)

ggplot(df, aes(x = Conditions, y = nums, fill = Conditions)) +
	geom_bar(stat = 'identity', colour = 'black') +
	labs(title = 'Number of comlete cases (without any NA)',
         x = 'Group',
         y = 'Complete cases') +
	scale_fill_manual(values = colsel(7, palette = 'sat1')) +
	theme_classic() + title +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here('exploration', "Complete_cases.pdf"), 
	width = 100, height = 90, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")



## histogram

met %>%
	group_by(Metabolite) %>%
	summarise(N = sum(!is.na(Score))) %>%
	ggplot(aes(N)) +
	geom_histogram() +
	labs(title = 'Number of non-NA values per metabolite',
         x = 'Number of non-NA values',
         y = 'Count') +
	theme_classic()

ggsave(file = here('exploration', "Hist_non_NA.pdf"), 
	width = 100, height = 90, units = 'mm', scale = 2, device = cairo_pdf, family = "Arial")






# correlation plots

no_heat = comp[complete.cases(comp %>% select(-HeatKill)),] %>% 
	select(-HeatKill) %>% data.frame

rownames(no_heat) = no_heat[,1]
no_heat[,1] = NULL

# correlation plot
corrplot::corrplot.mixed(cor(no_heat), upper = "square", lower = 'number')

dev.copy2pdf(device = cairo_pdf,
             file = here('exploration', 'corrplot_no_heat.pdf'),
             width = 8, height = 8, useDingbats = FALSE)

corrplot::corrplot.mixed(cor(comp[,2:8], use = "pairwise.complete.obs"), upper = "square", lower = 'number')

dev.copy2pdf(device = cairo_pdf,
             file = here('exploration', 'corrplot.pdf'),
             width = 8, height = 8, useDingbats = FALSE)



####################
### PCA analysis ###
####################




pca_b_data = met %>%
	# mutate(Score = replace_na(Score, 1E-20)) %>%
	mutate(Score = replace_na(Score, 4.23e-05)) %>%  # half of the minimum, a method that is used
	# filter(!Metabolite %in% removals) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F)
rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]

# lets compute the PCA
# warning, there are missing values of this

res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 6, graph = F)

# extract info about the individuals
ind = get_pca_ind(res.pca)


ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], ind$coord[,4], ind$coord[,5], ind$coord[,6])
ind_df['Sample'] = rownames(ind_df) 
ind_df = left_join(ind_df, metadata)
colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Dim4', 'Dim5', 'Dim6' ,'Sample', 'Group')

# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                  scale = c(sd(x) * sc, sd(y) * sc),
                                  centre = c(mean(x), mean(y))))
}

# make a data frame from ellipses
ell = ind_df %>% group_by(Group) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Group)) + 
  	geom_point(size = 4, show.legend = NA, alpha = 0.7) + 
  	geom_path(data = ell, aes(x = x, y = y, group = interaction(Group)), size = 1) +
  	geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Group), fill = Group), size = 1, alpha = 0.3) +
  	xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
  	ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
  	theme_classic() +
  	scale_colour_manual(values = colsel(7, palette = 'sat1')) +
  	scale_fill_manual(values = colsel(7, palette = 'sat1')) +
  	theme(plot.title = element_text(hjust = 0.5),
  	  panel.grid.major = element_blank(),
  	  panel.grid.minor = element_blank(),
  	  axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12))

quartz.save(type = 'pdf', 
	file = here('analysis', 'PCA_main.pdf'), 
	width = 6, height = 6, family = 'Arial')


# save PCA info
list_of_datasets = list('PCA score' = ind_df[,-6], 'Variance' = res.pca$eig[,2:3], 'Variable contribution' = res.pca$var$contrib)

write.xlsx(list_of_datasets, here('analysis', 'PCA_coord.xlsx'), colNames = T, rowNames = T) 





### secondary analyses from PCA data


# variables information
var = get_pca_var(res.pca)

### secondary plots
# plot variables
fviz_pca_var(res.pca, col.var = "black")

quartz.save(type = 'pdf', 
	file = here('analysis', 'PCA_variables.pdf'), 
	width = 9, height = 9)

# quality of representation
corrplot(var$cos2, is.corr = FALSE,  tl.col = "black", tl.cex = 0.3, tl.srt = 70)

quartz.save(type = 'pdf', 
	file = here('analysis', 'PCA_quality.pdf'), 
	width = 10, height = 49)

p1 = fviz_cos2(res.pca, choice = "var", axes = 1,  top = 50)
p2 = fviz_cos2(res.pca, choice = "var", axes = 2,  top = 50)
p3 = fviz_cos2(res.pca, choice = "var", axes = 3,  top = 50)
p6 = fviz_cos2(res.pca, choice = "var", axes = 4,  top = 50)
p7 = fviz_cos2(res.pca, choice = "var", axes = 5,  top = 50)
p4 = fviz_cos2(res.pca, choice = "var", axes = 1:2,top = 50)
p5 = fviz_cos2(res.pca, choice = "var", axes = 2:3,top = 50)


pdf(file = here('analysis', 'PCA_quality_dimensions.pdf'))
p1 
p2
p3
p6
p7
p4
p5
dev.off()


# Contributions of variables to PC1
p1 = fviz_contrib(res.pca, choice = "var", axes = 1, top = 50)

# Contributions of variables to PC2
p2 = fviz_contrib(res.pca, choice = "var", axes = 2, top = 50)

# Contributions of variables to PC3
p3 = fviz_contrib(res.pca, choice = "var", axes = 3, top = 50)

# Contributions of variables to PC3
p4 = fviz_contrib(res.pca, choice = "var", axes = 4, top = 50)

# Contributions of variables to PC3
p5 = fviz_contrib(res.pca, choice = "var", axes = 5, top = 50)


pdf(file = here( 'analysis' ,"/PCA_contrib.pdf"))
p1 
p2
p3
p4
p5
dev.off()




# contribution of individuals to PCA
p1 = fviz_contrib(res.pca, choice = "ind", axes = 1)
p2 = fviz_contrib(res.pca, choice = "ind", axes = 2)
p3 = fviz_contrib(res.pca, choice = "ind", axes = 3)
p4 = fviz_contrib(res.pca, choice = "ind", axes = 4)
p5 = fviz_contrib(res.pca, choice = "ind", axes = 5)

pdf(file = here('analysis', 'PCA_ind_contr.pdf'))
p1 
p2
p3
p4
p5
dev.off()


p1 = fviz_cos2(res.pca, choice = "ind", axes = 1,  top = 50)
p2 = fviz_cos2(res.pca, choice = "ind", axes = 2,  top = 50)
p3 = fviz_cos2(res.pca, choice = "ind", axes = 3,  top = 50)
p6 = fviz_cos2(res.pca, choice = "ind", axes = 4,  top = 50)
p7 = fviz_cos2(res.pca, choice = "ind", axes = 5,  top = 50)
p4 = fviz_cos2(res.pca, choice = "ind", axes = 1:2,top = 50)
p5 = fviz_cos2(res.pca, choice = "ind", axes = 2:3,top = 50)


pdf(file = here('analysis', 'PCA_ind_quality_dimensions.pdf'))
p1 
p2
p3
p6
p7
p4
p5
dev.off()


# corrplots are giving some problems with Illustrator, so here it is a new approach to plot these data

gradcolours = c('white', 'black')
data.frame(ind$cos2[,1:3]) %>%
  mutate(Sample = rownames(.)) %>%
  gather(Dimension, Values, Dim.1:Dim.3) %>%
  ggplot(aes(x = Dimension, y = Sample)) + 
  geom_tile(aes(fill = Values, width = Values, height = Values)) +
  xlab('Dimensions') + 
  ylab('Individuals') +
  ggtitle(paste("Quality of representation \n(Individuals)") )+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = NA),
        axis.line.x = element_line(colour = NA),
        axis.line.y = element_line(colour = NA),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_equal() +
  scale_fill_gradientn(colours = gradcolours, guide = "legend")

ggsave(file = here('analysis', 'PCA_ind_quality.pdf'), width = 42, height = 90, units = 'mm', scale = 2, device = cairo_pdf)


gradcolours = c('white', 'black')
data.frame(var$cos2[,1:3]) %>%
  mutate(Sample = rownames(.)) %>%
  gather(Dimension, Values, Dim.1:Dim.3) %>%
  ggplot(aes(x = Dimension, y = Sample)) + 
  geom_tile(aes(fill = Values, width = Values, height = Values)) +
  xlab('Dimensions') + 
  ylab('Individuals') +
  ggtitle(paste("Quality of representation \n(Variables)"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = NA),
        axis.line.x = element_line(colour = NA),
        axis.line.y = element_line(colour = NA),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_equal() +
  scale_fill_gradientn(colours = gradcolours, guide = "legend")

ggsave(file = here('analysis', 'PCA_var_quality.pdf'), width = 50, height = 200, units = 'mm', scale = 3, device = cairo_pdf)






#############################################
### Summary statistics, statistical tests ###
#############################################



# save all the summary statistics in one file
write.xlsx(filter(met.sum, !Metabolite %in% removals), here('analysis', 'Summary_stats.xlsx') , sheetName = c('results'), colNames = T, rowNames = F) 


# let's start with the statistical analysis
met.clean = filter(met, !Metabolite %in% removals)

# cleaning a bit the data set 
met.stat = met.clean %>% 
	select(Metabolite, Group, Sample, Replicate, Score) %>%
	mutate(Sample = str_replace_all(Sample,'[:digit:]', ''))


# test the general differences between groups and metabolites
model = aov(Score ~ Group * Metabolite, data = met.stat)
summary(model)



# test with subset of data
sub = met.stat %>% filter(Metabolite == 'Adenine')
model = lm('Score ~ Group', sub)
summary(model)
multi_model = multcomp::glht(model, linfct = multcomp::mcp(Group = 'Tukey'))
results = summary(multi_model)
results
# this is to have a tidy version of the statistical test
tidy(results)



# lets go for the complete dataset
# BEFORE THAT, WE HAVE REMOVE THOSE METABOLITES WITH ONLY ONE GROUP
# this is because, of course, we cannot 

stat.removals = met.sum %>% group_by(Metabolite) %>% summarise(N = sum(!is.na(Mean))) %>% 
	filter(N < 2) %>% select(Metabolite) %>%  t %>% as.vector

stat.removals = c(stat.removals, 'N-Acetylhistidine', 'Pyruvic acid', 'Tetrahydrouridine') # removing these ones because the are note good enough

nested = met.stat %>%
	filter(!Metabolite %in% stat.removals) %>%
    group_by(Metabolite) %>%
    nest()


# a couple of ad-hoc functions
apply_lm = function(df, formula){
      lm(data = df, as.formula(formula))      
}
apply_multilm = function(model){
	multcomp::glht(model, linfct = multcomp::mcp(Group = 'Tukey'))
}

# calculate linear models per metabolite
nested = nested %>%
      mutate(models = map(.f = apply_lm, .x = data, formula = 'Score ~ Group'))


# generate results and multiple hypothesis testing with glth
res = nested %>%
	mutate(multcomp = map(.f = apply_multilm, .x = models)) %>%
	mutate(results = map(.f = summary, .x = multcomp)) %>%
	mutate(results = map(.f = tidy, .x = results)) %>%
	dplyr::select(Metabolite, results) %>%
	unnest


# this code is to see which samples are not running well with the analysis
# for (i in 177:211) {
# 	print(i)
# 	summary(res$multcomp[i][[1]])
# }

# lets 
results = res %>% rename(Contrast = lhs) %>% 
	mutate(FDR = p.adjust(p.value, method = 'fdr'),
		   p_stars = gtools::stars.pval(p.value),
		   FDR_stars = gtools::stars.pval(FDR),
		   G1 = str_split(Contrast, ' - '),
		   G2 = str_split(Contrast, ' - '))

for (i in 1:2512) {
	results$G1[i] = unlist(results$G1[i][[1]][1])
	results$G2[i] = unlist(results$G2[i][[1]][2])
}

results$G1 = unlist(results$G1)
results$G2 = unlist(results$G2)

results = select(results, Metabolite, Contrast, G1, G2, everything())

# save all the summary statistics in one file
write.xlsx(results, here('analysis', 'statistical_comparisons.xlsx') , sheetName = c('stat'), colNames = T, rowNames = F) 




#####################
### Modulus stuff ###
#####################

mod.old = read_xlsx(here::here('data', 'modulus.xlsx'), sheet = 'old_worms')


# quick test to see wether there are differences between the three controls
w.cnt = mod.old %>% select(N2_day18_1, N2_day18_2, N2_day18_3) %>% gather(Sample, Score, N2_day18_1:N2_day18_3) %>% mutate(Sample = as.factor(Sample))
model = aov(Score ~ Sample, data = w.cnt)
tidy(TukeyHSD(model))

mod.stat = mod.old %>% gather(Sample, Score, everything()) %>% mutate(Sample = as.factor(Sample))

# create contrast
levels(mod.stat$Sample)
# 1:"Carb_day24"      2:"gltA_day23"   3:"HK_day24"    4:"N2_day18_1"   5:"N2_day18_2"   6:"N2_day18_3"   7:"N2_trypto_day18"   8:"TMP_day33"    9:"UV_day28"

M = rbind( c( 0 , 1 , 0 ,-1 , 0 , 0 , 0 , 0 , 0 ),     # gltA vs Control
		   c( 0 , 0 , 0 ,-1 , 0 , 0 , 0 , 1 , 0 ),     # TMP vs Control
		   c( 0 , 0 , 0 , 0 ,-1 , 0 , 1 , 0 , 0 ),     # Trypto vs Control
		   c( 0 , 0 , 1 , 0 , 0 ,-1 , 0 , 0 , 0 ),     # HK vs Control
		   c( 1 , 0 , 0 , 0 , 0 ,-1 , 0 , 0 , 0 ),     # Carb vs Control
		   c( 0 , 0 , 0 , 0 , 0 ,-1 , 0 , 0 , 1 ),     # UV vs Control
		   c( 0 , 1 , 0 , 0 , 0 , 0 , 0 ,-1 , 0 ),     # gltA vs TMP
		   c( 0 , 1 , 0 , 0 , 0 , 0 ,-1 , 0 , 0 ),     # gltA vs Trypto
		   c( 0 , 1 ,-1 , 0 , 0 , 0 , 0 , 0 , 0 ),     # gltA vs HK
		   c(-1 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ),     # gltA vs Carb
		   c( 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 ,-1 ),     # gltA vs UV
		   c( 0 , 0 , 0 , 0 , 0 , 0 ,-1 , 1 , 0 ),     # TMP vs Trypto
		   c( 0 , 0 ,-1 , 0 , 0 , 0 , 0 , 1 , 0 ),     # TMP vs HK
		   c(-1 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ),     # TMP vs Carb
		   c( 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,-1 ),     # TMP vs UV
		   c( 0 , 0 ,-1 , 0 , 0 , 0 , 1 , 0 , 0 ),     # Trypto vs HK
		   c(-1 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 ),     # Trypto vs Carb
		   c( 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ,-1 ),     # Trypto vs UV
		   c(-1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 ),     # HK vs Carb
		   c( 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,-1 ),     # HK vs UV
		   c( 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,-1 )     # Carb vs UV
		   )

# Optional: we can name each column
colnames(M) = c("Carb_day24", "gltA_day23", "HK_day24", "N2_day18_1", "N2_day18_2" , "N2_day18_3" , "N2_trypto_day18" , "TMP_day33" , "UV_day28")

# Not-so-optional: name what is each comparison
rownames(M) = c('gltA vs Control',
				'TMP vs Control',
				'Trypto vs Control',
				'HK vs Control',
				'Carb vs Control',
				'UV vs Control',
				'gltA vs TMP',
				'gltA vs Trypto',
				'gltA vs HK',
				'gltA vs Carb',
				'gltA vs UV',
				'TMP vs Trypto',
				'TMP vs HK',
				'TMP vs Carb',
				'TMP vs UV',
				'Trypto vs HK',
				'Trypto vs Carb',
				'Trypto vs UV',
				'HK vs Carb',
				'HK vs UV',
				'Carb vs UV')


model = lm('Score ~ Sample', mod.stat)
summary(model)
multi_model = multcomp::glht(model, linfct = multcomp::mcp(Sample = M))
results = summary(multi_model)
tidy(results)








##############################
### gltA, TMP, UV analysis ###
##############################

# Here's a short description of what this piece of code does, just for clarity and 
# to refresh my bad memory in the future.
# First of all, I look for all the specific metabolites that are shared between the 
# control (Bacto) and each one of the other groups. Then, I perform an intersect to 
# get the shared metabolites between some of the groups to create groups of groups. 
# I mean, the UV, TM and gltA behave in similar ways for the worm phenotype, so I take
# all the metabolites that are shared between them. After each group is created, I check
# that they all have means and values, and are complete cases

# Aside from this, I have calculated the log2 of each value, summarised them by their means 
# and standard deviations, and then, filtered only the control values (the means) to store 
# them as a new column in the general data set (met). By doing this, we can calculate 
# the fold change of each metabolite just substrating the control value. 
# Remember that {log(x) - log(y) = log(x/y)}. 

# Also, I have calculated the multi univariate comparisons between all metabolites per
# condition using the log transformed data, and corrected the multiple comparisons 
# with FDR values. 

# With all this, we can now proceed. First, we filter only the metabolites present
# in the group of interest, and then, we filter only those that does not have any 
# difference between them (as the group is formed of conditions that behave similarly).
# Then, I also filter the p-values of this subset of metabolites in respect to the 
# control, and I make a Heatmap with all these information. 
# Summarising the heatmap information, each square represents the FC(log2) between 
# each condition in respect to the control, per metabolite which seems to be 
# equal between conditions. The stars represent the p-value of the condition-control 
# comparison.



temp = met.sum %>% filter(Group %in% c('ControlBacto', 'Uvirradiation')) %>%
			ungroup %>%
			select(Group, Metabolite, Mean) %>%
			spread(Group, Mean)

temp = temp[complete.cases(temp),]
UV.met = temp$Metabolite

temp = met.sum %>% filter(Group %in% c('ControlBacto', 'gltAdel')) %>%
			ungroup %>%
			select(Group, Metabolite, Mean) %>%
			spread(Group, Mean)

temp = temp[complete.cases(temp),]
gltA.met = temp$Metabolite

temp = met.sum %>% filter(Group %in% c('ControlBacto', 'Trimethoprim')) %>%
			ungroup %>%
			select(Group, Metabolite, Mean) %>%
			spread(Group, Mean)
temp = temp[complete.cases(temp),]
TM.met = temp$Metabolite


temp = met.sum %>% filter(Group %in% c('ControlBacto', 'Carbenicillin')) %>%
			ungroup %>%
			select(Group, Metabolite, Mean) %>%
			spread(Group, Mean)

temp = temp[complete.cases(temp),]
carb.met = temp$Metabolite

temp = met.sum %>% filter(Group %in% c('ControlBacto', 'ControlTrypt')) %>%
			ungroup %>%
			select(Group, Metabolite, Mean) %>%
			spread(Group, Mean)

temp = temp[complete.cases(temp),]
tryp.met = temp$Metabolite

temp = met.sum %>% filter(Group %in% c('ControlBacto', 'HeatKill')) %>%
			ungroup %>%
			select(Group, Metabolite, Mean) %>%
			spread(Group, Mean)

temp = temp[complete.cases(temp),]
HK.met = temp$Metabolite

# metabolites common for the three conditions (plus Bacto)
# we will use it later
UAT.met = intersect(intersect(UV.met, gltA.met), TM.met)

# metabolites in Carb and Trypto
CT.met = intersect(tryp.met, carb.met)

# metabolites in Carb, Trypto and HK
CTH.met = intersect(intersect(tryp.met, HK.met), carb.met)

# check that this works

met.sum %>% filter(Metabolite %in% UV.met, Group == 'Uvirradiation') %>% data.frame

met.sum %>% filter(Metabolite %in% gltA.met, Group == 'gltAdel') %>% data.frame

met.sum %>% filter(Metabolite %in% TM.met, Group == 'Trimethoprim') %>% data.frame

met.sum %>% filter(Metabolite %in% carb.met, Group == 'Carbenicillin') %>% data.frame

met.sum %>% filter(Metabolite %in% tryp.met, Group == 'ControlTrypt') %>% data.frame

met.sum %>% filter(Metabolite %in% HK.met, Group == 'HeatKill') %>% data.frame

met.sum %>% filter(Metabolite %in% UAT.met, Group %in% c('Uvirradiation', 'gltAdel', 'Trimethoprim')) %>% data.frame

met.sum %>% filter(Metabolite %in% UAT.met, Group %in% c('ControlBacto')) %>% data.frame


# create 
met = met %>% mutate(logScore = log(Score, base = 2))


# summarise data
met.sum = met %>%
	group_by(ID, Group, Metabolite, PubChem, HMDB) %>%
	summarise(Mean = mean(Score, na.rm = TRUE),
			  SD = sd(Score, na.rm = TRUE),
			  logMean = mean(logScore, na.rm = TRUE),
			  logSD = sd(logScore, na.rm = TRUE))


control = met.sum %>% ungroup %>% 
	filter(Group == 'ControlBacto') %>% 
	select(Metabolite, Mean, logMean) %>%
	rename(Bacto_logMean = logMean,
			Bacto_Mean = Mean)


met = met %>% left_join(control) %>% ungroup
# after all these little transformations, lets rebuild again everything

UAT.sum = met %>%
	ungroup %>%
	# filter(Group %in% c('gltAdel', 'Trimethoprim', 'Uvirradiation')) %>%
	mutate(logFC = logScore - Bacto_logMean) %>%
	group_by(ID, Group, Metabolite, PubChem, HMDB) %>%
	summarise(Mean_FC = mean(logFC, na.rm = TRUE),
			  SD_FC = sd(logFC, na.rm = TRUE)) %>%
	ungroup







########
## statistical comparison for the log-transformed data

# cleaning a bit the data set 
met.stat = met %>% 
	select(Metabolite, Group, Sample, Replicate, logScore) %>%
	mutate(Sample = str_replace_all(Sample,'[:digit:]', ''))


# test the general differences between groups and metabolites
model = aov(logScore ~ Group * Metabolite, data = met.stat)
summary(model)


stat.removals = c(stat.removals, 'N-Acetylhistidine', 'Pyruvic acid', 'Tetrahydrouridine') # removing these ones because the are note good enough

nested = met.stat %>%
	filter(!Metabolite %in% stat.removals) %>%
    group_by(Metabolite) %>%
    nest()

# calculate linear models per metabolite
nested = nested %>%
      mutate(models = map(.f = apply_lm, .x = data, formula = 'logScore ~ Group'))


# generate results and multiple hypothesis testing with glth
res = nested %>%
	mutate(multcomp = map(.f = apply_multilm, .x = models)) %>%
	mutate(results = map(.f = summary, .x = multcomp)) %>%
	mutate(results = map(.f = tidy, .x = results)) %>%
	dplyr::select(Metabolite, results) %>%
	unnest


# this code is to see which samples are not running well with the analysis
# for (i in 177:211) {
# 	print(i)
# 	summary(res$multcomp[i][[1]])
# }

# lets 
log.res = res %>% rename(Contrast = lhs) %>% 
	mutate(FDR = p.adjust(p.value, method = 'fdr'),
		   p_stars = gtools::stars.pval(p.value),
		   FDR_stars = gtools::stars.pval(FDR),
		   G1 = str_split(Contrast, ' - '),
		   G2 = str_split(Contrast, ' - '))

for (i in 1:2512) {
	log.res$G1[i] = unlist(log.res$G1[i][[1]][1])
	log.res$G2[i] = unlist(log.res$G2[i][[1]][2])
}

log.res$G1 = unlist(log.res$G1)
log.res$G2 = unlist(log.res$G2)

log.res = select(log.res, Metabolite, Contrast, G1, G2, everything())

# correct for values = 0
log.res[log.res$FDR == 0,]$FDR = 1e-20

# save all the summary statistics in one file
write.xlsx(log.res, here('analysis', 'stat_log_comparisons.xlsx') , sheetName = c('stat'), colNames = T, rowNames = F) 




####################
######

# once we have the statistical comparisons between groups, we can proceed
# to plot the data in heatmaps, scatterplots and volcano plots

#########################
### gltA, TM and UV group

met.sum %>% filter(Metabolite %in% UAT.met, Group %in% c('Uvirradiation', 'gltAdel', 'Trimethoprim')) %>% data.frame



UAT.sum	%>% 
	filter(Metabolite %in% UAT.met, Group %in% c('Uvirradiation', 'gltAdel', 'Trimethoprim'))



## lets select only those metabolites that are not different between the three combinations of groups

contr = c('Uvirradiation - gltAdel', 'Trimethoprim - gltAdel', 'Uvirradiation - Trimethoprim')

temp = log.res %>%
	filter(Metabolite %in% UAT.met, Contrast %in% contr) %>%
	select(Metabolite, Contrast, FDR) %>%
	spread(Contrast, FDR) %>% data.frame

rownames(temp) = temp[,1]; temp[,1] = NULL

# this complicated and almost unnecessary command is 
# to extract the not significant metabolites of the three group combination

not.sig.UAT = rownames(temp[rowSums(temp > 0.05) == 3,])


# with this subset, we will  extract those contrasts against the bacto

stat.description = log.res %>% 
	filter(Metabolite %in% not.sig.UAT, Contrast %in% c('gltAdel - ControlBacto', 'Uvirradiation - ControlBacto', 'Trimethoprim - ControlBacto')) %>%
	select(Metabolite, G1, FDR, FDR_stars) %>%
	rename(Group = G1)





minv <- -2
maxv <- 2
nstep <- maxv-minv
nstep <-8
clrbrks = seq(-2, 2 ,by = 0.5)
brks = seq(minv, maxv , by = (maxv-minv)/(nstep))
bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)

UAT.sum	%>% 
	filter(Metabolite %in% not.sig.UAT, Group %in% c('Uvirradiation', 'gltAdel', 'Trimethoprim')) %>% 
	filter(!Metabolite %in% c('XA0027', 'Terephthalic acid', 'Rhein')) %>%
	left_join(stat.description) %>%
	mutate(Metabolite = as.factor(Metabolite)) %>%
	mutate(Metabolite = factor(Metabolite, levels = rev(levels(Metabolite)))) %>%
	ggplot(aes(x = Group, y = Metabolite)) +
	geom_tile(aes(fill = Mean_FC)) +
	theme_minimal() +
	scale_fill_gradientn(colours = clrscale, breaks = clrbrks, limits = c(-2,2)) +
	geom_text(aes(label = as.character(FDR_stars)), size = 6) +
	theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12, angle = 45, hjust = 1),
  	    axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
	file = here('analysis', '/Comparison_heatmap_gltA_TM_UV.pdf'), 
	width = 6, height = 7, family = 'Arial')





##############
### Carb, tryp
# check the data
met.sum %>% filter(Metabolite %in% CT.met, Group %in% c('Carbenicillin', 'ControlTrypt')) %>% data.frame



UAT.sum	%>% 
	filter(Metabolite %in% CT.met, Group %in% c('Carbenicillin', 'ControlTrypt'))



## lets select only those metabolites that are not different between the three combinations of groups

contr = c('ControlTrypt - Carbenicillin')

not.sig.CT = log.res %>%
	filter(Metabolite %in% CT.met, Contrast %in% contr) %>%
	select(Metabolite, Contrast, FDR) %>%
	spread(Contrast, FDR) %>% 
	filter(`ControlTrypt - Carbenicillin` > 0.05) %>%
	select(Metabolite) %>% t %>% as.vector

# this complicated and almost unnecessary command is 
# to extract the not significant metabolites of the three group combination


# with this subset, we will  extract those contrasts against the bacto

stat.description = log.res %>% 
	filter(Metabolite %in% not.sig.CT, Contrast %in% c('ControlTrypt - ControlBacto', 'ControlBacto - Carbenicillin')) %>%
	select(Metabolite, G1, G2, FDR, FDR_stars) %>%
	mutate(Group = ifelse(G2 == 'Carbenicillin', 'Carbenicillin', 'ControlTrypt'))





minv <- -2.1
maxv <- 2
nstep <- maxv-minv
nstep <-8
clrbrks = seq(-2.1, 2 ,by = 0.5)
brks = seq(minv, maxv , by = (maxv-minv)/(nstep))
bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)

UAT.sum	%>% 
	filter(Metabolite %in% not.sig.CT, Group %in% c('Carbenicillin', 'ControlTrypt')) %>% 
	filter(!Metabolite %in% c('XA0002', 'XA0033', 'XC0017', 'XC0137', 'Terephthalic acid')) %>%
	left_join(stat.description) %>%
	mutate(Metabolite = as.factor(Metabolite)) %>%
	mutate(Metabolite = factor(Metabolite, levels = rev(levels(Metabolite)))) %>%
	ggplot(aes(x = Group, y = Metabolite)) +
	geom_tile(aes(fill = Mean_FC)) +
	theme_minimal() +
	scale_fill_gradientn(colours = clrscale, breaks = clrbrks, limits = c(-2.1,2)) +
	geom_text(aes(label = as.character(FDR_stars))) +
	theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12, angle = 45, hjust = 1),
  	    axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 


quartz.save(type = 'pdf', 
	file = here('analysis', '/Comparison_heatmap_Carb_Tryp.pdf'), 
	width = 7, height = 11, family = 'Arial')




###############
### Carb, tryp, HK
# check the data
met.sum %>% filter(Metabolite %in% CTH.met, Group %in% c('Carbenicillin', 'ControlTrypt', 'HeatKill')) %>% data.frame



UAT.sum	%>% 
	filter(Metabolite %in% CTH.met, Group %in% c('Carbenicillin', 'ControlTrypt', 'HeatKill'))



## lets select only those metabolites that are not different between the three combinations of groups

contr = c('ControlTrypt - Carbenicillin', 'HeatKill - ControlTrypt', 'HeatKill - Carbenicillin')


temp = log.res %>%
	filter(Metabolite %in% CTH.met, Contrast %in% contr) %>%
	select(Metabolite, Contrast, FDR) %>%
	spread(Contrast, FDR) %>% data.frame

rownames(temp) = temp[,1]; temp[,1] = NULL

# this complicated and almost unnecessary command is 
# to extract the not significant metabolites of the three group combination

not.sig.CTH = rownames(temp[rowSums(temp > 0.05) == 3,])

# with this subset, we will  extract those contrasts against the bacto

stat.description = log.res %>% 
	filter(Metabolite %in% not.sig.CTH, Contrast %in% c('ControlTrypt - ControlBacto', 'ControlBacto - Carbenicillin', 'HeatKill - ControlBacto')) %>%
	select(Metabolite, G1, FDR, FDR_stars) %>%
	rename(Group = G1)

stat.description$Group = ifelse(stat.description$Group == 'ControlBacto', 'Carbenicillin', stat.description$Group)




minv <- -0.65
maxv <- 1.1
nstep <- maxv-minv
nstep <-8
clrbrks = seq(-0.65, 1.1 ,by = 0.4)
brks = seq(minv, maxv , by = (maxv-minv)/(nstep))
bgg <- colorRampPalette(c("blue", "gray90", "red"))(n = nstep)
clrscale <- colorRampPalette(c("blue4","blue", "gray90", "red","red4"))(n = nstep)

UAT.sum	%>% 
	filter(Metabolite %in% not.sig.CTH, Group %in% c('Carbenicillin', 'ControlTrypt', 'HeatKill')) %>% 
	left_join(stat.description) %>%
	filter(!Metabolite == 'Terephthalic acid') %>%
	mutate(Metabolite = as.factor(Metabolite)) %>%
	mutate(Metabolite = factor(Metabolite, levels = rev(levels(Metabolite)))) %>%
	ggplot(aes(x = Group, y = Metabolite)) +
	geom_tile(aes(fill = Mean_FC)) +
	theme_minimal() +
	scale_fill_gradientn(colours = clrscale, breaks = clrbrks, limits = c(-0.65, 1.1)) +
	geom_text(aes(label = as.character(FDR_stars))) +
		theme(axis.ticks=element_blank(),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12, angle = 45, hjust = 1),
  	    axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	    axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 


quartz.save(type = 'pdf', 
	file = here('analysis', '/Comparison_heatmap_Carb_Tryp_HK.pdf'), 
	width = 7, height = 3.5, family = 'Arial')












##########################################
### scatterplots from different groups ###
##########################################

scat = function(data, mets, groups) {
    cosa = data %>% filter(Metabolite %in% mets, Group %in% groups) %>%
		ungroup %>%
		select(Group, Metabolite, Mean_FC, SD_FC) %>%
		gather(Stat, Value, -(Group:Metabolite)) %>%
		unite(Var, Group, Stat) %>%
		spread(Var, Value) 
	colnames(cosa)[2:5] = c('x', 'xsd', 'y', 'ysd')
	
	cosa %>% 
		ggplot(aes_string(x = x, y = y)) +
		geom_errorbar(aes_string(ymin = y - ysd, ymax = y + ysd), colour = "grey50", alpha = 0.5) +
		geom_errorbarh(aes_string(xmin = x - xsd, xmax = x + xsd), colour = "grey50", alpha = 0.5) +
		geom_point(alpha = 0.8, size = 1.5) +
		theme_classic() +
		labs(
    		title = 'gltA vs TM',
    		x = 'gltA vs Control, FC(log2)',
    		y = 'Trimethoprim vs Control, FC(log2)') +
    	geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) + 
    	geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    	geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", alpha = 0.8)
}



scat(UAT.sum, mets = UAT.met, group = c('gltAdel', 'Trimethoprim'))









UAT.sum %>% filter(Metabolite %in% not.sig, Group %in% c('gltAdel', 'Trimethoprim')) %>%
	ungroup %>%
	select(Group, Metabolite, Mean_FC, SD_FC) %>%
	gather(Stat, Value, -(Group:Metabolite)) %>%
	unite(Var, Group, Stat) %>%
	spread(Var, Value) %>%
	ggplot(aes(x = gltAdel_Mean_FC, y = Trimethoprim_Mean_FC)) +
	geom_errorbar(aes(ymin = Trimethoprim_Mean_FC - Trimethoprim_SD_FC, ymax = Trimethoprim_Mean_FC + Trimethoprim_SD_FC), colour = "grey50", alpha = 0.5) +
	geom_errorbarh(aes(xmin = gltAdel_Mean_FC - gltAdel_SD_FC, xmax = gltAdel_Mean_FC + gltAdel_SD_FC), colour = "grey50", alpha = 0.5) +
	geom_point(alpha = 0.8, size = 1.5) +
	theme_classic() +
	labs(
    	title = 'gltA vs TM',
    	x = 'gltA vs Control, FC(log2)',
    	y = 'Trimethoprim vs Control, FC(log2)') +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", alpha = 0.8)




########################
### Complete heatmap ###
########################


library(ComplexHeatmap)
library(circlize)

heatdata = read_xlsx(here::here('data', 'Metabolites.xlsx'), sheet = 'Heatmap')

heatdata = heatdata %>%
	rename(Metabolite = `Compound name`) %>%
	select(Metabolite, CT1:P4) %>% data.frame


# some loops to clean the rubish 
for (i in 1:dim(heatdata)[1]) { heatdata$Metabolite[i] = strsplit(heatdata$Metabolite[i], "\\r")[[1]][1] }


rownames(heatdata) = heatdata[,1]; heatdata[,1] = NULL


Heatmap(heatdata, col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")),
	name = 'Values',
	column_title = "Samples", 
	column_title_side = "bottom",
	row_title = "Metabolites",
	show_row_names = FALSE)


quartz.save(type = 'pdf', 
	file = here('analysis', 'Samples_heatmap.pdf'), 
	width = 11, height = 11, family = 'Arial', dpi = 500)


cairo_ps(filename = here('analysis', 'Samples_heatmap.eps'),
         width = 9, height = 9, pointsize = 12,
         fallback_resolution = 500)
Heatmap(heatdata, col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")),
	name = 'Values',
	column_title = "Samples", 
	column_title_side = "bottom",
	row_title = "Metabolites",
	show_row_names = FALSE)
dev.off()






#####################
### Volcano plots ###
#####################


# Volcano plots with density contour

### gltA

volcan = log.res %>% 
	filter(Metabolite %in% UV.met, Contrast == 'Uvirradiation - ControlBacto') %>% 
	select(Metabolite, estimate:FDR_stars) %>% 
	left_join(UAT.sum %>% 
		filter(Metabolite %in% UV.met, Group == 'Uvirradiation') %>% 
		select(Metabolite, Mean_FC)) %>%
	rename(logFC = Mean_FC) %>%
	mutate(log10P = -log10(FDR))

# simple volcano plot

volcan %>%
	ggplot(aes(x = logFC, y = log10P)) +
	expand_limits(y = c(0, 1)) +
	geom_vline(xintercept = 0) +
	geom_hline(yintercept = 0) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	geom_hline(yintercept = -log10(0.01), linetype = "longdash") +
	geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
	geom_point(data = subset(volcan, abs(logFC) < 2), alpha = 0.8) +
	geom_point(data = subset(volcan, abs(logFC) >= 2 & log10P > -log10(0.01)), colour = 'red', size = 2.1) +
	geom_text_repel(data = subset(volcan, abs(logFC) >= 2 & log10P > -log10(0.01)), aes(label = Metabolite)) +
	theme_bw() + 
	theme(panel.grid = element_blank()) +
	xlab("Fold change (log2)") +
	ylab("-log10(FDR)")




# second version


ggplot(volcan, aes(x = logFC, y = log10P)) +
	scale_fill_gradient(low = "lightgray", high = "navy") +
	scale_color_gradient(low = "lightgray", high = "navy") +
	# scale_y_continuous(trans = revlog_trans(), expand = c(0.005, 0.005)) +
	stat_density_2d(aes(fill = ..level..), geom = "polygon", show.legend = FALSE) +
	geom_point(data = subset(volcan, abs(logFC) >= 2 & log10P > -log10(0.01)), 
		colour = 'red', size = 2.1) +
	geom_vline(xintercept = 0) +
	geom_hline(yintercept = 0) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	geom_hline(yintercept = -log10(0.01), linetype = "longdash") +
	geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
	theme_linedraw() + 
	theme(panel.grid = element_blank()) +
	geom_text_repel(data = subset(volcan, abs(logFC) >= 2 & log10P > -log10(0.01)), aes(label = Metabolite)) +
	xlab("Fold change (log2)") +
	ylab("-log10(FDR)")



# function to plot every condition

vol2d = function(met = gltA.met, group = 'gltAdel', contrast = 'gltAdel - ControlBacto', lim = 2, pval = 0.05) {
	volcan = log.res %>% 
	filter(Metabolite %in% met, Contrast == 'gltAdel - ControlBacto') %>% 
	select(Metabolite, estimate:FDR_stars) %>% 
	left_join(UAT.sum %>% 
		filter(Metabolite %in% met, Group == group) %>% 
		select(Metabolite, Mean_FC)) %>%
	rename(logFC = Mean_FC) %>%
	mutate(log10P = -log10(FDR))
	# return(volcan)

	p = ggplot(volcan, aes_string(x = 'logFC', y = 'log10P')) +
		scale_fill_gradient(low = "lightgray", high = "navy") +
		scale_color_gradient(low = "lightgray", high = "navy") +
		# scale_y_continuous(trans = revlog_trans(), expand = c(0.005, 0.005)) +
		stat_density_2d(aes_string(fill = '..level..'), geom = "polygon", show.legend = FALSE) +
		geom_point(data = subset(volcan, abs(logFC) >= lim & log10P > -log10(pval)), 
			colour = 'red', size = 2.1) +
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = 0) +
		geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
		geom_hline(yintercept = -log10(0.01), linetype = "longdash") +
		geom_vline(xintercept = c(-lim, lim), linetype = "dashed") +
		theme_linedraw() + 
		theme(panel.grid = element_blank()) +
		geom_text_repel(data = subset(volcan, abs(logFC) >= lim & log10P > -log10(pval)), aes_string(label = 'Metabolite')) +
		xlab("Fold change (log2)") +
		ylab("-log10(FDR)")
	return(p)
}

# generate and save plots

# gltA
vol2d(met = gltA.met, group = 'gltAdel', contrast = 'gltAdel - ControlBacto', lim = 0, pval = 0.01)

quartz.save(type = 'pdf', 
	file = here('analysis', 'Volcano_2d_gltA.pdf'), 
	width = 9, height = 9, family = 'Arial')

# UV
vol2d(met = UV.met, group = 'Uvirradiation', contrast = 'Uvirradiation - ControlBacto',lim = 0, pval = 0.01)
quartz.save(type = 'pdf', 
	file = here('analysis', 'Volcano_2d_UV.pdf'), 
	width = 9, height = 9, family = 'Arial')

# TM
vol2d(met = TM.met, group = 'Trimethoprim', contrast = 'Trimethoprim - ControlBacto',lim = 0, pval = 0.01)
quartz.save(type = 'pdf', 
	file = here('analysis', 'Volcano_2d_TM.pdf'), 
	width = 9, height = 9, family = 'Arial')


# Carb
vol2d(met = carb.met, group = 'Carbenicillin', contrast = 'ControlBacto - Carbenicillin', lim = 0, pval = 0.01)
quartz.save(type = 'pdf', 
	file = here('analysis', 'Volcano_2d_Carb.pdf'), 
	width = 9, height = 9, family = 'Arial')


# Tryp
vol2d(met = tryp.met, group = 'ControlTrypt', contrast = 'ControlTrypt - ControlBacto', lim = 0, pval = 0.01)
quartz.save(type = 'pdf', 
	file = here('analysis', 'Volcano_2d_Tryp.pdf'), 
	width = 9, height = 9, family = 'Arial')


# HK
vol2d(met = HK.met, group = 'HeatKill', contrast = 'HeatKill - ControlBacto', lim = 0, pval = 0.01)
quartz.save(type = 'pdf', 
	file = here('analysis', 'Volcano_2d_HK.pdf'), 
	width = 9, height = 9, family = 'Arial')




######




###############################
###############################

## EXTRA: PCA and heatmap without epsilon

## PCA
na_met = met %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F) 



pca_b_data = met %>%
	# mutate(Score = replace_na(Score, 1E-20)) %>%
	# filter(!Metabolite %in% removals) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F)
rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]

# lets compute the PCA
# warning, there are missing values of this

res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 6, graph = F)

# extract info about the individuals
ind = get_pca_ind(res.pca)


ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], ind$coord[,4], ind$coord[,5], ind$coord[,6])
ind_df['Sample'] = rownames(ind_df) 
ind_df = left_join(ind_df, metadata)
colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Dim4', 'Dim5', 'Dim6' ,'Sample', 'Group')

# get ellipses based on the correlation
getellipse = function(x, y, sc = 1) {
  as.data.frame(ellipse::ellipse(cor(x, y),
                                  scale = c(sd(x) * sc, sd(y) * sc),
                                  centre = c(mean(x), mean(y))))
}

# make a data frame from ellipses
ell = ind_df %>% group_by(Group) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Group)) + 
  	geom_point(size = 4, show.legend = NA, alpha = 0.7) + 
  	geom_path(data = ell, aes(x = x, y = y, group = interaction(Group)), size = 1) +
  	geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Group), fill = Group), size = 1, alpha = 0.3) +
  	xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
  	ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
  	theme_classic() +
  	scale_colour_manual(values = colsel(7, palette = 'sat1')) +
  	scale_fill_manual(values = colsel(7, palette = 'sat1')) +
  	theme(plot.title = element_text(hjust = 0.5),
  	  panel.grid.major = element_blank(),
  	  panel.grid.minor = element_blank(),
  	  axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12))

quartz.save(type = 'pdf', 
	file = here('exploration', 'PCA_main_withNAs.pdf'), 
	width = 6, height = 6, family = 'Arial')

## Missing values are imputed by the mean of the variable




####


## previous version
final_met = met %>%
	mutate(Score = replace_na(Score, 2E-52)) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F) 
met.scale = t(scale(t(final_met), center = TRUE, scale = TRUE))


Heatmap(met.scale, col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")),
	name = 'Values',
	column_title = "Samples", 
	column_title_side = "bottom",
	row_title = "Metabolites",
	show_row_names = FALSE)





na_met = met %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>% as.matrix
	data.frame(check.names = F)

na_met2 = matrix(as.numeric(na_met), nrow = dim(na_met)[1], ncol = dim(na_met)[2])
colnames(na_met2) = colnames(na_met)
rownames(na_met2) = na_met[,1]



cosa = t(scale((na_met2), center = TRUE, scale = TRUE))
cosa = cosa[-1,]

# remove every row with everything as NA values, produced because 
# there was only one sample with a value in the raw data
# cosa = cosa[rowSums(is.na(cosa)) != ncol(cosa), ]

# Im having a lot of problems trying to order the rows with clustering due
# to the high quantity of NAs
# So, let's  try to cluster the data with epsilons, and then, pass that clustering 
# to the matrix without epsilons

mat = hclust(amap::Dist(met.scale, method = 'euclidean'), method = 'complete')
mat2 = hclust(amap::Dist(t(met.scale), method = 'euclidean'), method = 'complete')

# mat = hclust(dist(met.scale), method = 'complete')

columns = c('HK1', 'HK4', 'HK2', 'HK3',
			'TP3', 'TP2', 'TP4', 'TP1',
			'UV1', 'UV3', 'UV4', 'UV2',
			'P2', 'P1', 'P4', 'P3',
			'GL1', 'GL2', 'GL4', 'GL3',
			'CT4', 'CT3', 'CT1', 'CT2',
			'CB3', 'CB2', 'CB4', 'CB1')

Heatmap(cosa, col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")),
	name = 'Values',
	column_title = "Samples", 
	column_title_side = "bottom",
	row_title = "Metabolites",
	show_row_names = FALSE,
	cluster_rows = mat,
	cluster_columns = F,
	column_order = columns)

quartz.save(type = 'pdf', 
	file = here('exploration', 'Heatmap_withNAs.pdf'), 
	width = 11, height = 11, family = 'Arial')




# ----------------------------------------------------------#


############################################
### Different strategies for NAs in data ###
############################################


### half of the value of minimum

pca_b_data = met %>%
	# mutate(Score = replace_na(Score, 1E-20)) %>%
	mutate(Score = replace_na(Score, 4.23e-05)) %>%  # half of the minimum, a method that is used
	# filter(!Metabolite %in% removals) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F)
rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]



### half of the value of minimum

pca_b_data = met %>%
	# filter(!Metabolite %in% removals) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F)
rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]

# loop to substitute every NA by the minimum value per metabolite
for (i in 1:228) {
	x = pca_b_data[,i]
	x[is.na(x)] = min(na.omit(x))/2
	pca_b_data[,i] = x
}


### Random Forest imputation

library(missForest)

pca_b_data = met %>%
	# mutate(Score = replace_na(Score, 1E-20)) %>%
	# mutate(Score = replace_na(Score, 4.23e-05)) %>%  # half of the minimum, a method that is used
	# filter(!Metabolite %in% removals) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F)


rf = missForest(pca_b_data)
pca_b_data = rf$ximp


rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]



### SVD (singular value decomposition) method

library(pcaMethods)
pca_b_data = met %>%
	# filter(!Metabolite %in% removals) %>%
	select(Sample, Metabolite, Score) %>%
	spread(Metabolite, Score) %>%
	data.frame(check.names = F)
rownames(pca_b_data) = pca_b_data[,1]; pca_b_data = pca_b_data[,-1]


result <- pca(pca_b_data, method = "svdImpute" , nPcs = 3, center = TRUE) 
pca_b_data <- completeObs(result)







## quick heatmap
met.scale = t(scale((pca_b_data), center = TRUE, scale = TRUE))
met.scale = met.scale[-2,]
Heatmap(met.scale, col = colorRamp2(c(-3, 0, 3), c("green", "white", "red")),
	name = 'Values',
	column_title = "Samples", 
	column_title_side = "bottom",
	row_title = "Metabolites",
	show_row_names = FALSE)





# MICE method
colnames(pca_b_data) = 1:228
imputed_Data = mice(t(pca_b_data), m=5, maxit = 50, method = 'pmm', seed = 500)

complete_data = complete(imputed_Data,2)

pca_b_data = t(complete_data)












# lets compute the PCA
# warning, there are missing values of this

res.pca = PCA(pca_b_data, scale.unit = TRUE, ncp = 6, graph = F)

# extract info about the individuals
ind = get_pca_ind(res.pca)


ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3], ind$coord[,4], ind$coord[,5], ind$coord[,6])
ind_df['Sample'] = rownames(ind_df) 
ind_df = left_join(ind_df, metadata)
colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3', 'Dim4', 'Dim5', 'Dim6' ,'Sample', 'Group')

# make a data frame from ellipses
ell = ind_df %>% group_by(Group) %>% do(getellipse(.$Dim1, .$Dim2, 1)) %>% data.frame

# plot!
ggplot(ind_df, aes(x = Dim1, y = Dim2, color = Group)) + 
  	geom_point(size = 4, show.legend = NA, alpha = 0.7) + 
  	geom_path(data = ell, aes(x = x, y = y, group = interaction(Group)), size = 1) +
  	geom_polygon(data = ell, aes(x = x, y = y, group = interaction(Group), fill = Group), size = 1, alpha = 0.3) +
  	xlab(paste("PC1 - ", round(res.pca$eig[1,2], 1), " % of variance", sep = "")) + 
  	ylab(paste("PC2 - ", round(res.pca$eig[2,2], 1), " % of variance", sep = "")) +
  	theme_classic() +
  	scale_colour_manual(values = colsel(7, palette = 'sat1')) +
  	scale_fill_manual(values = colsel(7, palette = 'sat1')) +
  	theme(plot.title = element_text(hjust = 0.5),
  	  panel.grid.major = element_blank(),
  	  panel.grid.minor = element_blank(),
  	  axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
  	  axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12))

quartz.save(type = 'pdf', 
	file = here('exploration/NA_methods', 'PCA_half_min_permetabolite.pdf'), 
	width = 6, height = 6, family = 'Arial')

























