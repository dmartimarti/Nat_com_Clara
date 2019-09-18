# BGA analysis

library(tidyverse)
library(readxl)
library(ggrepel)
library(here)
library(openxlsx)

options(width = 220)

# my own library of functions
source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')


odir = 'Plots'
dir.create(odir, showWarnings = TRUE, recursive = FALSE, mode = "0777")

# directory 
# "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Filipe/Metab_natcom/data/BGA"


# Get timeseries data
# It will take a while
time.data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
	filter(Data %in% c('595nm_f', '750nm_f')) %>%
	select(-c(Variable1, Replicate_x)) %>%
	gather(Time_s, OD, matches('\\d')) %>%
	filter(!is.na(Media)) %>% # Remove empty values if there are missmatches
	rename(Str = Media) %>%
	mutate(Str = as.character(Str), # Change strain namings
           Time_s = as.numeric(Time_s),
           Time_h = Time_s/3600,
           Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
           Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
           Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
           Col = factor(Col, levels = LETTERS[1:8])) %>%
 	select(c(-File, -Pattern, -Strain)) %>%
 	rename(Strain = Str, Replicate = Replicate_y) %>%
 	mutate_at(c('Well', 'Strain'), as.factor)



tsum = time.data %>%
	group_by(Data, Strain, Time_h) %>%
	summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
	ungroup 


# write all time series data for paper revision

write.xlsx(tsum, 'FigureS5_growth_curves.xlsx', colNames = T, rowNames = F) 




# generate labels for plots

od595 = tsum %>% filter(Data == '595nm_f', Time_h == 18) %>% select(Strain, Time_h, Mean) %>% mutate(Time_h = Time_h + 0.85)
od750 = tsum %>% filter(Data == '750nm_f', Time_h == 18) %>% select(Strain, Time_h, Mean) %>% mutate(Time_h = Time_h + 0.85)

od595[od595$Strain == 'UV',]$Time_h = 18.3
od595[od595$Strain == 'Carb',]$Time_h = 19.5
od595[od595$Strain == 'HK',]$Time_h = 18.8

# complete plot
tsum %>%
  filter(Data == '595nm_f', Time_h > 0.2) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  # facet_wrap(vars(Strain)) +
  scale_colour_manual(values = c('#B31714', '#FFFF12', '#FF8E8C', '#1881CC', '#6893B3', '#B314A5', '#FF8CF5', '#18CCB9', '#3D8880', '#FF8012', '#3814B3', '#FF8012', '#06CC10', '#3D8741')) +
  scale_fill_manual(  values = c('#B31714', '#FFFF12', '#FF8E8C', '#1881CC', '#6893B3', '#B314A5', '#FF8CF5', '#18CCB9', '#3D8880', '#FF8012', '#3814B3', '#FF8012', '#06CC10', '#3D8741')) +
  geom_text(data = od595, aes(label = Strain, x = Time_h, y = Mean)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'BGA_595_complete.pdf'), 
  width = 10, height = 8, family = 'Arial')




group1 = c('Control', 'UV', 'TMP', 'TMPR', 'HK', 'TMPR_TMP', 'Carb')
group2 = c('CtrSoy', 'CtrTryp', 'Control', 'sucA', 'gltA')
group3 = c('Control', 'Met', 'MR', 'MR_Met')

# partial plots
# group 1
tsum %>%
  # mutate(Strain = as.character(Strain)) %>%
  filter(Data == '595nm_f', Time_h > 0.2) %>%
  filter(Strain %in% group1) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  # ylim(0,0.7) +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  # facet_wrap(vars(Strain)) +
  scale_colour_manual(values = c('blue','black','green','red','gray50','orange','brown')) +
  scale_fill_manual(values =   c('blue','black','green','red','gray50','orange','brown')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'BGA_595_HK_Carb_UV.pdf'), 
  width = 10, height = 8, family = 'Arial')



# partial plots
# group 2
tsum %>%
  # mutate(Strain = as.character(Strain)) %>%
  filter(Data == '595nm_f', Time_h > 0.2) %>%
  filter(Strain %in% group2) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  # facet_wrap(vars(Strain)) +
  scale_colour_manual(values = c('black', '#2253B3','green','red','gray50','orange')) +
  scale_fill_manual(values =   c('black', '#2253B3','green','red','gray50','orange')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'BGA_595_sucA_gltA_Controls.pdf'), 
  width = 10, height = 8, family = 'Arial')






# partial plots
# group 3
tsum %>%
  # mutate(Strain = as.character(Strain)) %>%
  filter(Data == '595nm_f', Time_h > 0.2) %>%
  filter(Strain %in% group3) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  # facet_wrap(vars(Strain)) +
  scale_colour_manual(values = c('black', '#2253B3','gray50','red','gray50','orange')) +
  scale_fill_manual(values =   c('black', '#2253B3','gray50','red','gray50','orange')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'BGA_595_Metformin.pdf'), 
  width = 10, height = 8, family = 'Arial')



### same for 750nm


# partial plots
# group 1
tsum %>%
  # mutate(Strain = as.character(Strain)) %>%
  filter(Data == '750nm_f', Time_h > 0.2) %>%
  filter(Strain %in% group1) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  # ylim(0,0.7) +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  # facet_wrap(vars(Strain)) +
  scale_colour_manual(values = c('blue','black','green','red','gray50','orange','brown')) +
  scale_fill_manual(values =   c('blue','black','green','red','gray50','orange','brown')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'BGA_750_HK_Carb_UV.pdf'), 
  width = 10, height = 8, family = 'Arial')



# partial plots
# group 2
tsum %>%
  # mutate(Strain = as.character(Strain)) %>%
  filter(Data == '750nm_f', Time_h > 0.2) %>%
  filter(Strain %in% group2) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  # facet_wrap(vars(Strain)) +
  scale_colour_manual(values = c('black', '#2253B3','green','red','gray50','orange')) +
  scale_fill_manual(values =   c('black', '#2253B3','green','red','gray50','orange')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'BGA_750_sucA_gltA_Controls.pdf'), 
  width = 10, height = 8, family = 'Arial')






# partial plots
# group 3
tsum %>%
  # mutate(Strain = as.character(Strain)) %>%
  filter(Data == '750nm_f', Time_h > 0.2) %>%
  filter(Strain %in% group3) %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  # facet_wrap(vars(Strain)) +
  scale_colour_manual(values = c('black', '#2253B3','gray50','red','gray50','orange')) +
  scale_fill_manual(values =   c('black', '#2253B3','gray50','red','gray50','orange')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'BGA_750_Metformin.pdf'), 
  width = 10, height = 8, family = 'Arial')







######################
######################
######################


data = read_csv('Output/Summary.csv', quote = "\"") %>%
    filter(!is.na(Media)) %>%
    select(-c(Variable1, Replicate_x, Strain, Metformin_mM, File, Pattern, Data, Reader)) %>%
    rename(logAUC750_raw = `750nm_f_logAUC`,
           AUC750_raw = `750nm_f_AUC`, 
           logAUC595_raw = `595nm_f_logAUC`,
           AUC595_raw = `595nm_f_AUC`,
           Strain = Media,
           Replicate = Replicate_y) %>% 
    mutate(
        AUC750_scale =  AUC750_raw * 2, 
        AUC595_scale =  AUC595_raw * 2,
        logAUC750_raw = log2(AUC750_scale),
        logAUC595_raw = log2(AUC595_scale),
        Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
        Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
        Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
        Col = factor(Col, levels = LETTERS[1:8])) %>% #Change Type column coding
    select(Strain, Replicate, Well, Col, Row, logAUC595_raw, logAUC750_raw, AUC595_scale, AUC750_scale)  


# dataset to include also the growth rate of bacteria
data.gr = read_csv('Output/Summary.csv', quote = "\"") %>%
    filter(!is.na(Media)) %>%
    select(-c(Variable1, Replicate_x, Strain, Metformin_mM, File, Pattern, Data, Reader)) %>%
    rename(logAUC750_raw = `750nm_f_logAUC`,
           AUC750_raw = `750nm_f_AUC`, 
           logAUC595_raw = `595nm_f_logAUC`,
           AUC595_raw = `595nm_f_AUC`,
           GR595_raw = `595nm_dt_Max`,
           Strain = Media,
           Replicate = Replicate_y) %>% 
    mutate(
        AUC750_scale =  AUC750_raw * 2, 
        AUC595_scale =  AUC595_raw * 2,
        logAUC750_raw = log2(AUC750_scale),
        logAUC595_raw = log2(AUC595_scale),
        Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
        Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
        Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
        Col = factor(Col, levels = LETTERS[1:8])) %>% #Change Type column coding
    select(Strain, Replicate, Well, Col, Row, logAUC595_raw, logAUC750_raw, AUC595_scale, GR595_raw)  



data750 = data[1:42,]  %>% select(c(-logAUC595_raw, -AUC595_scale)) %>% mutate(AUC750_scale = AUC750_scale/2)

data595 = data[43:84,] %>% select(c(-logAUC750_raw, -AUC750_scale)) %>% mutate(AUC595_scale = AUC595_scale/2)


d595.sum = data595 %>% 
    group_by(Strain) %>% 
    summarise(Mean = mean(AUC595_scale), SD = sd(AUC595_scale))

d750.sum = data750 %>% 
    group_by(Strain) %>% 
    summarise(Mean = mean(AUC750_scale), SD = sd(AUC750_scale))





# summary of the growth rate of bacteria
gr595.sum = data.gr %>%
  group_by(Strain) %>%
  summarise(Mean = mean(GR595_raw, na.rm = TRUE), SD = sd(GR595_raw, na.rm = TRUE))




# write data from growth rate in excel

list_of_datasets = list('Growth Rate Summary' = gr595.sum, 'Growth Rate raw' = data.frame(data.gr)[43:84,])

write.xlsx(list_of_datasets,'Growth_rate.xlsx', colNames = T, rowNames = F) 



########################
### statistical test ###
########################

library(broom)


model = aov((AUC595_scale/2) ~ Strain, data = data595)

results = tidy(TukeyHSD(model))

contrasts = c('Met-Control', 'MR-Met', 'MR-Control', 'MR_Met-Control', 'MR_Met-Met', 
    'TMP-Control', 'TMPR-Control', 'TMPR_TMP-Control', 'TMPR_TMP-TMPR', 'TMPR_TMP-TMP', 
    'Control-Carb', 'CtrSoy-Control', 'CtrTryp-Control', 'gltA-Control', 'HK-Control', 
    'Met-Control', 'sucA-Control', 'UV-Control')


res.595 = results %>% 
    filter(comparison %in% c(contrasts)) %>%
    mutate(p_stars = gtools::stars.pval(adj.p.value))
    



model = aov((AUC750_scale/2) ~ Strain, data = data750)

results = tidy(TukeyHSD(model))

contrasts = c('Met-Control', 'MR-Met', 'MR-Control', 'MR_Met-Control', 'MR_Met-Met', 
    'TMP-Control', 'TMPR-Control', 'TMPR_TMP-Control', 'TMPR_TMP-TMPR', 'TMPR_TMP-TMP', 
    'Control-Carb', 'CtrSoy-Control', 'CtrTryp-Control', 'gltA-Control', 'HK-Control', 
    'Met-Control', 'sucA-Control', 'UV-Control')


res.750 = results %>% 
    filter(comparison %in% c(contrasts)) %>%
    mutate(p_stars = gtools::stars.pval(adj.p.value))


# save all the summary statistics in one file

list_of_datasets <- list('raw_data_595' = data595, 'summary_stats_595' = d595.sum, "anova_595" = res.595, 
                         'raw_data_750' = data750, 'summary_stats_750' = d750.sum, "anova_750" = res.750)


write.xlsx(list_of_datasets, 'stats.xlsx', colNames = T, rowNames = F) 









##############################################################################################################
##############################################################################################################


######################################################
### statistical test with endpoints of time series ###
######################################################



od595 = time.data %>% filter(Data == '595nm_f', Time_h == 18) %>% select(Strain, Replicate, OD)
od750 = time.data %>% filter(Data == '750nm_f', Time_h == 18) %>% select(Strain, Replicate, OD)


model = aov(OD ~ Strain, data = od595)

results = tidy(TukeyHSD(model))

contrasts = c('Met-Control', 'MR-Met', 'MR-Control', 'MR_Met-Control', 'MR_Met-Met', 
    'TMP-Control', 'TMPR-Control', 'TMPR_TMP-Control', 'TMPR_TMP-TMPR', 'TMPR_TMP-TMP', 
    'Control-Carb', 'CtrSoy-Control', 'CtrTryp-Control', 'gltA-Control', 'HK-Control', 
    'Met-Control', 'sucA-Control', 'UV-Control')


res.595 = results %>% 
    filter(comparison %in% c(contrasts)) %>%
    mutate(p_stars = gtools::stars.pval(adj.p.value))
    


# model = aov(OD ~ Strain, data = od750)

# results = tidy(TukeyHSD(model))

# contrasts = c('Met-Control', 'MR-Met', 'MR-Control', 'MR_Met-Control', 'MR_Met-Met', 
#     'TMP-Control', 'TMPR-Control', 'TMPR_TMP-Control', 'TMPR_TMP-TMPR', 'TMPR_TMP-TMP', 
#     'Control-Carb', 'CtrSoy-Control', 'CtrTryp-Control', 'gltA-Control', 'HK-Control', 
#     'Met-Control', 'sucA-Control', 'UV-Control')


# res.750 = results %>% 
#     filter(comparison %in% c(contrasts)) %>%
#     mutate(p_stars = gtools::stars.pval(adj.p.value))


# sub = od750 %>% filter(Strain %in% c('Control', 'gltA'))

# od595.sum = tsum %>% filter(Data == '595nm_f', Time_h == 18)
# od750.sum = tsum %>% filter(Data == '750nm_f', Time_h == 18)

# ## barplot
# od750.sum %>% 
#     ggplot(aes(x = Strain, y = Mean)) +
#     geom_bar(stat = 'identity') +
#     geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD))







##############
### corrected values for end point


corr_val = c(  0.434, 0.916, 0.687, 0.516, 0.397, 0.200, 0.042, 0.042, 
    0.357, 0.925, 0.665, 0.533, 0.448, 0.217, 0.043, 0.043,
    0.413, 0.925, 0.653, 0.547, 0.402, 0.146, 0.042, 0.042,
    0.223, 0.301, 0.044, 0.632, 0.377, 0.298,
    0.242, 0.313, 0.044, 0.642, 0.418, 0.284,
    0.213, 0.328, 0.044, 0.645, 0.404, 0.290)



od750rep = od750
od750rep$OD = corr_val

model = aov(OD ~ Strain, data = od750rep)

results = tidy(TukeyHSD(model))

contrasts = c('Met-Control', 'MR-Met', 'MR-Control', 'MR_Met-Control', 'MR_Met-Met', 
    'TMP-Control', 'TMPR-Control', 'TMPR_TMP-Control', 'TMPR_TMP-TMPR', 'TMPR_TMP-TMP', 
    'Control-Carb', 'CtrSoy-Control', 'CtrTryp-Control', 'gltA-Control', 'HK-Control', 
    'Met-Control', 'sucA-Control', 'UV-Control')


res.750 = results %>% 
    filter(comparison %in% c(contrasts)) %>%
    mutate(p_stars = gtools::stars.pval(adj.p.value))




# save all the summary statistics in one file

list_of_datasets <- list('raw_data_595' = od595, 'summary_stats_595' = od595.sum, "anova_595" = res.595, 
                         'raw_data_750' = od750rep, 'summary_stats_750' = od750p.sum, "anova_750" = res.750)


write.xlsx(list_of_datasets, 'stats_OD.xlsx', colNames = T, rowNames = F) 




################
#### barplot ###
################


od750p.sum = od750rep %>% group_by(Strain) %>%
    summarise(Mean = mean(OD), SD = sd(OD)) 

od750p.sum$Strain <- factor(od750p.sum$Strain, 
    levels = c('Control', 'Met', 'MR', 'MR_Met', 'TMP', 'TMPR', 'TMPR_TMP',
        'CtrSoy', 'CtrTryp', 'sucA', 'gltA', 'HK', 'Carb', 'UV'))

od750p.sum %>%
    ggplot(aes(x = Strain, y = Mean)) +
    geom_bar(stat = 'identity', aes(fill = Strain), colour = 'black') +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
    ylab("OD") +
    xlab("Sample") +
    ylim(0,1) +
    scale_colour_manual(values = colsel(14, palette = 'total')) +
    scale_fill_manual(values =   colsel(14, palette = 'total')) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'Barplot_750_OD.pdf'), 
  width = 14, height = 9, family = 'Arial')




od595.sum = od595 %>% group_by(Strain) %>%
    summarise(Mean = mean(OD), SD = sd(OD)) 

od595.sum$Strain <- factor(od595.sum$Strain, 
    levels = c('Control', 'Met', 'MR', 'MR_Met', 'TMP', 'TMPR', 'TMPR_TMP',
        'CtrSoy', 'CtrTryp', 'sucA', 'gltA', 'HK', 'Carb', 'UV'))

od595.sum %>%
    ggplot(aes(x = Strain, y = Mean)) +
    geom_bar(stat = 'identity', aes(fill = Strain), colour = 'black') +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
    ylab("OD") +
    xlab("Sample") +
    scale_y_continuous(limits = c(-0.05, 0.70), breaks = seq(0,0.7, by = 0.1)) +
    scale_colour_manual(values = colsel(14, palette = 'total')) +
    scale_fill_manual(values =   colsel(14, palette = 'total')) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('Plots', 'Barplot_595_OD.pdf'), 
  width = 14, height = 9, family = 'Arial')


list_of_datasets <- list('OD_595' = od595, 'OD_750' = od750rep)

write.xlsx(list_of_datasets, 'raw_OD_values.xlsx', colNames = T, rowNames = F) 




##########################
### interaction effect ###
##########################


## Metformin
## OD 750
# data preparation
inter = od750rep %>% 
    filter(Strain %in% c('Control', 'Met', 'MR', 'MR_Met')) %>%
    mutate(Strain = as.character(Strain))

inter['Metformin'] = 1

inter[inter$Strain == 'Control' | inter$Strain == 'MR',]$Metformin = 0
inter[inter$Strain == 'Met',]$Strain = 'Control'
inter[inter$Strain == 'MR_Met',]$Strain = 'MR'


inter = inter %>% select(Strain, Metformin, Replicate, OD) %>%
    mutate_at(c('Strain', 'Metformin'), as.factor)

model = aov(OD ~ Strain * Metformin, data = inter)

met750.inter = tidy(model) %>% mutate(p_stars = gtools::stars.pval(p.value))


## OD 595
# data preparation
inter = od595 %>% 
    filter(Strain %in% c('Control', 'Met', 'MR', 'MR_Met')) %>%
    mutate(Strain = as.character(Strain))

inter['Metformin'] = 1

inter[inter$Strain == 'Control' | inter$Strain == 'MR',]$Metformin = 0
inter[inter$Strain == 'Met',]$Strain = 'Control'
inter[inter$Strain == 'MR_Met',]$Strain = 'MR'


inter = inter %>% select(Strain, Metformin, Replicate, OD) %>%
    mutate_at(c('Strain', 'Metformin'), as.factor)

model = aov(OD ~ Strain * Metformin, data = inter)

met595.inter = tidy(model) %>% mutate(p_stars = gtools::stars.pval(p.value))


## TMP
## OD 750

inter = od750rep %>% 
    filter(Strain %in% c('Control', 'TMP', 'TMPR', 'TMPR_TMP')) %>%
    mutate(Strain = as.character(Strain))

inter['TMP'] = 1

inter[inter$Strain == 'Control' | inter$Strain == 'TMPR',]$TMP = 0
inter[inter$Strain == 'TMP',]$Strain = 'Control'
inter[inter$Strain == 'TMPR_TMP',]$Strain = 'TMPR'


inter = inter %>% select(Strain, TMP, Replicate, OD) %>%
    mutate_at(c('Strain', 'TMP'), as.factor)

model = aov(OD ~ Strain * TMP, data = inter)

tmp750.inter = tidy(model) %>% mutate(p_stars = gtools::stars.pval(p.value))



# model = lm(OD ~ Strain * TMP, data = inter)



## TMP
## OD 595

inter = od595 %>% 
    filter(Strain %in% c('Control', 'TMP', 'TMPR', 'TMPR_TMP')) %>%
    mutate(Strain = as.character(Strain))

inter['TMP'] = 1

inter[inter$Strain == 'Control' | inter$Strain == 'TMPR',]$TMP = 0
inter[inter$Strain == 'TMP',]$Strain = 'Control'
inter[inter$Strain == 'TMPR_TMP',]$Strain = 'TMPR'


inter = inter %>% select(Strain, TMP, Replicate, OD) %>%
    mutate_at(c('Strain', 'TMP'), as.factor)

model = aov(OD ~ Strain * TMP, data = inter)

tmp595.inter = tidy(model) %>% mutate(p_stars = gtools::stars.pval(p.value))






# save all the summary statistics in one file

list_of_datasets <- list('Met_595' = met595.inter, 'Met_750' = met750.inter, 
                         'TMP_595' = tmp595.inter, 'TMP_750' = tmp750.inter)


write.xlsx(list_of_datasets, 'Met_TMP_interactions.xlsx', colNames = T, rowNames = F) 









###################################
### extra Bacterial Growh Assay ###
###   final data for paper      ###   
###################################



# path: "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Filipe/Metab_natcom"



library(tidyverse)
library(readxl)
library(ggrepel)
library(here)
library(openxlsx)

source('/Users/dmarti14/Documents/MRC_Postdoc/scripts/R_functions/all_functions.R')

# session options
options(width = 220)


time.data = read_xlsx('data/Bacterial growth curve.xlsx', sheet = 'Summary') %>%
    gather(Time_m, OD, matches('\\d')) %>%
    mutate(
           Time_m = as.numeric(Time_m),
           Time_h = Time_m/60) %>%
    mutate_at(c('Day', 'Strain'), as.factor) %>%
    select(Day, Strain, Replicate, Time_m, Time_h, OD)

# test
cosa = time.data %>% filter(Day == 'Day1', Strain == 'N2', Replicate == 1) %>% data.frame
auc(cosa$Time_h, cosa$OD)

auc_data = time.data %>%
    group_by(Day, Strain, Replicate) %>%
    summarise(AUC = auc(Time_h, OD),
        AUC_log = log(AUC, base = 2))


# save data into excel file
list_of_datasets <- list('AUC' = auc_data)

write.xlsx(list_of_datasets, 'BGA_AUC.xlsx', colNames = T, rowNames = F) 


# summary
tsum = time.data %>%
  group_by(Day, Strain, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
  ungroup 





# complete plot
tsum %>%
  ggplot(aes(x = Time_h, y = Mean, fill = Strain, color = Strain)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Strain") +
  #facet_wrap(~MetaboliteU,ncol = 2)+
  facet_wrap(vars(Day)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.text.y = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.x = element_text(family = 'Arial', colour = 'black', size = 12),
        axis.title.y = element_text(family = 'Arial', colour = 'black', size = 12)) 

quartz.save(type = 'pdf', 
  file = here('analysis', 'BGA_extra.pdf'), 
  width = 18, height = 8, family = 'Arial')








###################################
###   Bacterial Growth curves   ###
###   final data for paper      ###   
###################################


library(growthrates)


time.data

## data preparation: delete some time-points that are a bit shitty
# try to delete the first and last timepoints 
time.data2 = time.data[time.data$Time_h > 0.4 & time.data$Time_h < 15,]
# delete TMP rep 3, and replace it with a corrected one
TMP.corrected = time.data2[time.data2$Strain == 'TMP' & time.data2$Replicate == 3 & time.data2$Time_h >= 2,]
time.data2 = time.data2 %>% filter(!(Strain == 'TMP' & Replicate == 3)) %>% rbind(TMP.corrected)
# correct first value of CtrSoy
time.data2[time.data2$Strain == 'CtrSoy' & time.data2$Replicate == 3, ][1:2,]$OD = 0.00120

# delete first 3 timepoints from Strain == 'CtrTryp', Replicate == 1
CtrTryp.corrected = time.data2[time.data2$Strain == 'CtrTryp' & time.data2$Replicate == 1 & time.data2$Time_h >= 1,]
time.data2 = time.data2 %>% filter(!(Strain == 'CtrTryp' & Replicate == 1)) %>% rbind(CtrTryp.corrected)

# some tests
test = time.data2 %>% filter(Data == '595nm_f',Strain == 'Control', Replicate == 1)
fit = fit_easylinear(test$Time_h, test$OD, h = 10, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)


test2 = time.data2 %>% filter(Data == '595nm_f',Strain == 'CtrSoy', Replicate == 3)
fit = fit_easylinear(test2$Time_h, test2$OD, h = 10, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)




easyfit = function(df, h){
  fit = fit_easylinear(df$Time_h, df$OD, h = h)
  return(coef(fit)[3])
}

# nested data
nested = time.data2 %>%
  filter(Data == '595nm_f') %>%
  mutate(OD = ifelse(OD == 0, 0.001, OD)) %>%
  filter(!Strain %in% c('HK', 'UV', 'Carb', 'MR', 'MR_Met', 'TMPR', 'TMPR_TMP')) %>%
  group_by(Strain, Replicate) %>%
  nest() 

fits = nested %>%
  mutate(fit = map(data, easyfit, 15))

fits = fits %>% 
  select(Strain, Replicate, mu_max = fit) %>% 
  unnest() %>%
  arrange(Strain)

fits %>%
  ggplot(aes(x = Strain, y = fit)) +
  geom_point() +
  theme_classic()


fits %>%
  group_by(Strain) %>%
  summarise(Mean = mean(mu_max),
            SD = sd(mu_max)) %>%
  ggplot(aes(Strain, Mean)) +
  geom_pointrange(aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_classic()

fits.sum = fits %>%
  group_by(Strain) %>%
  summarise(Mean = mean(mu_max),
            SD = sd(mu_max)) 

# save statistics list
list_of_datasets = list('Growth rate' = fits, 'Growth rate summary' = fits.sum)

write.xlsx(list_of_datasets, 'Growth_rate.xlsx', colNames = T, rowNames = F) 



###############
### TESTING ZONE



# # check that everything works as usual
# for (i in 1:42){
#   print(paste0('Analysing dataset number: ', i))
#   print(paste0('Bacteria ', nested[i,]$Strain, ', Replicate ', nested[i,2]))
#   easyfit(nested$data[[i]])
# }



# Control plots
test = time.data2 %>% filter(Data == '595nm_f',Strain == 'CtrSoy', Replicate == 1)
fit = fit_easylinear(test$Time_h, test$OD, h = 12, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)

test = time.data2 %>% filter(Data == '595nm_f',Strain == 'CtrSoy', Replicate == 2)
fit = fit_easylinear(test$Time_h, test$OD, h = 15, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)

test = time.data2 %>% filter(Data == '595nm_f',Strain == 'CtrSoy', Replicate == 3)
fit = fit_easylinear(test$Time_h, test$OD, h = 15, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)




# CtrSoy plots
test = time.data2 %>% filter(Data == '595nm_f',Strain == 'CtrTryp', Replicate == 1)
fit = fit_easylinear(test$Time_h, test$OD, h = 13, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)

test = time.data2 %>% filter(Data == '595nm_f',Strain == 'CtrTryp', Replicate == 2)
fit = fit_easylinear(test$Time_h, test$OD, h = 13, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)

test = time.data2 %>% filter(Data == '595nm_f',Strain == 'CtrTryp', Replicate == 3)
fit = fit_easylinear(test$Time_h, test$OD, h = 13, quota = 0.95)
summary(fit)
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)












########################
# manual way of calculating growth rate


gr = function(df, window = 5, tp = 4){
    dim(test)[1] / window
} 
























