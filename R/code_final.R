library(smatr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(car)
library(rotl)
library(MCMCglmm)
library(postMCMCglmm)
library(phytools)
library(nlme)
library(geiger)
library(phylolm)
library(MuMln)
library(lme4)
library(stringr)
library(TreePar)
library(HDInterval)

#Read the data frame in. The data file can be found as a .csv in the data folder. One could also download the data from the appendix of the article, just make sure you use the right variables. Here are what the variables in "survival_size_full_data_March_6_20.csv" are called in the appendix:
#"mass.at.hatch..g." = "mo"
#"Size.at.maturity..g." = "ma"
#"Clutch.size*Clutches.yr" = "N" (Note that you need to multiple these two variables to get N, as done below)
#"juvenile.survival.to.adulthood" = "Sa"

surv <- read.csv("~/Documents/research/completed/relative offspring size:survival/survival_size_full_data_March_6_20.csv", stringsAsFactor = F) 
head(surv)


#These are nested ANOVAs to calculate proportions of variance associated with different taxonomic levels.
mod_nanova_mo <- aov(log10(mass.at.hatch..g.) ~ group + (order%in%group) + (fam_%in%order) , data = surv)
summary(mod_nanova_order)
sstot_mo <- sum(summary(mod_nanova_mo)[[1]][,2])
ssclass_mo<- summary(mod_nanova_mo)[[1]][1,2]
ssorder_mo<- summary(mod_nanova_mo)[[1]][2,2]
ssfam_mo<- summary(mod_nanova_mo)[[1]][3,2]
ssresid_mo <- summary(mod_nanova_mo)[[1]][4,2]
ssclass_mo/sstot_mo
ssorder_mo/sstot_mo
ssfam_mo/sstot_mo
ssresid_mo/sstot_mo

mod_nanova_ma <- aov(log10(Size.at.maturity..g.) ~ group + (order%in%group) + (fam_%in%order) , data = surv)
summary(mod_nanova_ma)
sstot_ma <- sum(summary(mod_nanova_ma)[[1]][,2])
ssclass_ma<- summary(mod_nanova_ma)[[1]][1,2]
ssorder_ma<- summary(mod_nanova_ma)[[1]][2,2]
ssfam_ma<- summary(mod_nanova_ma)[[1]][3,2]
ssresid_ma <- summary(mod_nanova_ma)[[1]][4,2]
ssclass_ma/sstot_ma
ssorder_ma/sstot_ma
ssfam_ma/sstot_ma
ssresid_ma/sstot_ma

mod_nanova_N <- aov(log10(Clutch.size*Clutches.yr) ~ group + (order%in%group) + (fam_%in%order) , data = surv)
summary(mod_nanova_N)
sstot_N <- sum(summary(mod_nanova_N)[[1]][,2])
ssclass_N<- summary(mod_nanova_N)[[1]][1,2]
ssorder_N<- summary(mod_nanova_N)[[1]][2,2]
ssfam_N<- summary(mod_nanova_N)[[1]][3,2]
ssresid_N <- summary(mod_nanova_N)[[1]][4,2]
ssclass_N/sstot_N
ssorder_N/sstot_N
ssfam_N/sstot_N
ssresid_N/sstot_N


mod_nanova_S <- aov(log10(juvenile.survival.to.adulthood) ~ group + (order%in%group) + (fam_%in%order) , data = surv)
summary(mod_nanova_S)
sstot_S <- sum(summary(mod_nanova_S)[[1]][,2])
ssclass_S<- summary(mod_nanova_S)[[1]][1,2]
ssorder_S<- summary(mod_nanova_S)[[1]][2,2]
ssfam_S<- summary(mod_nanova_S)[[1]][3,2]
ssresid_S <- summary(mod_nanova_S)[[1]][4,2]
ssclass_S/sstot_S
ssorder_S/sstot_S
ssfam_S/sstot_S
ssresid_S/sstot_S

mod_nanova_order <- aov(log10(juvenile.survival.to.adulthood) ~ group + (order%in%group) + (fam_%in%order) , data = surv)
summary(mod_nanova_order)
etaSquared(mod_nanova_order)
mod_nanova_order <- aov(log10(mo.M) ~ group + (order%in%group) + (fam_%in%order) , data = surv)
summary(mod_nanova_order)


mod_nanova_order <- aov(log10(Clutch.size*Clutches.yr) ~ group + (order%in%group) + (fam_%in%order) , data = surv)
summary(mod_nanova_order)


#These are synonyms/typo fixes to link the life history data to species in the Open Tree of Life
surv$Species[surv$Species == "Geochelone gigantea"] <- "Aldabrachelys gigantea"
surv$Species[surv$Species == "Cnemidophorus sexlineatus"] <- "Aspidoscelis sexlineatus"
surv$Species[surv$Species == "Cnemidophorus uniparens"] <- "Aspidoscelis uniparens"
surv$Species[surv$Species == "Cnemidophorus tigris"] <- "Aspidoscelis tigris"
surv$Species[surv$Species == "Lacerta vivipara"] <- "Zootoca vivipara"
surv$Species[surv$Species == "Spermophilus lateralis"] <- "Callospermophilus lateralis"
surv$Species[surv$Species == "Spermophilus beldingi"] <- "Urocitellus beldingi"
surv$Species[surv$Species == "Clethrionomys glareolus"] <- "Myodes glareolus"
surv$Species[surv$Species == "Masticophis taeniatus"] <- "Coluber taeniatus"
surv$Species[surv$Species == "Hylarana guentheri"] <- "Sylvirana guentheri"
surv$Species[surv$Species == "Lithobates catesbeianus"] <- "Rana catesbeiana"
surv$Species[surv$Species == "Lithobates pipiens"] <- "Rana pipiens"
surv$Species[surv$Species == "Lithobates sevosus"] <- "Rana sevosa"
surv$Species[surv$Species == "Glandirana rugosa"] <- "Rugosa rugosa"
surv$Species[surv$Species == "Lithobates sylvaticus"] <- "Rana sylvaticus"
surv$Species[surv$Species == "Crotaphytus wislizeni"] <- "Gambelia sila"
surv$Species[surv$Species == "Spermophilus arniatus"] <- "Spermophilus relictus"
surv$Species[surv$Species == "Aepyceros metampus"] <- "Aepyceros melampus"
surv$Species[surv$Species == "Ailuropoda melanokuca"] <- "Ailuropoda melanoleuca"
surv$Species[surv$Species == "Myocaster coypus"] <- "Myocastor coypus"
surv$Species[surv$Species == "Sus scrota"] <- "Sus scrofa"
surv$Species[surv$Species == "Phoca hispida"] <- "Pusa hispida"
surv$Species[surv$Species == "Lutra canadensis"] <- "Lontra canadensis"
surv$Species[surv$Species == "Felts catus"] <- "Felis catus"
surv$Species[surv$Species == "Hyla versicolor"] <- "Dryophytes versicolor"
surv$Species[surv$Species == "Coluber taeniatus"] <- "Masticophis taeniatus"
surv$Species[surv$Species == "Rana sylvaticus"] <- "Rana sylvatica"
surv$Species[surv$Species == "Engraulis encrasicolus"] <- "Engraulis encrasicolus"
surv$Species[surv$Species == "Peromyscus maniculatus"] <- "Peromyscus maniculatus"
surv$Species[surv$Species == "Myocastor coypus"] <- "Myocastor coypus"
surv$Species[surv$Species == "Zapus hudsonius"] <- "Zapus hudsonius"
surv$Species[surv$Species == "Urosaurus ornatus"] <- "Urosaurus ornatus"
surv$Species[surv$Species == "Chrysemys picta"] <-  "Chrysemys picta"
surv$Species[surv$Species == "Kinosternon subrubum"] <-  "Kinosternon subrubrum"
surv$Species[surv$Species == "Phrynosoma douglasi"] <-  "Phrynosoma douglasii"
surv$Species[surv$Species == "Podocnemis voglii"] <-  "Podocnemis vogli"
surv$Species[surv$Species == "Sceloporus jarrovi"] <-  "Sceloporus jarrovii"
surv$Species[surv$Species == "Sceloporus poinsetti"] <-  "Sceloporus poinsettii"

surv_pgls_order <- surv  %>% group_by(order) %>% filter(n()>2) %>% ungroup() %>% mutate(animal = tolower(Species), log10_survival = log10(juvenile.survival.to.adulthood), log10_mo_M = log10(mo.M), log10_mo = log10(mass.at.hatch..g.), log10_M = log10(Size.at.maturity..g.), log10_N = log10(Clutch.size*Clutches.yr), log10_R = log10(Clutch.size*Clutches.yr*mass.at.hatch..g.), log10_L = log10(adult.lifespan)) %>% filter(!is.na(juvenile.survival.to.adulthood), !is.na(mo.M)) %>% arrange(order) %>% data.frame()

surv_pgls_order_lvl_data <-  surv%>% group_by(order)%>% summarize(animal = tolower(first(order)), log10_survival = log10(mean(juvenile.survival.to.adulthood)), log10_mo_M = log10(mean(mo.M)), log10_mo = log10(mean(mass.at.hatch..g.)), log10_M = log10(mean(Size.at.maturity..g.)), log10_N = log10(mean(Clutch.size*Clutches.yr)), log10_R = log10(mean(Clutch.size*Clutches.yr*mass.at.hatch..g.)), log10_L = log10(mean(adult.lifespan)), log10_C = log10(mean(juvenile.survival.to.adulthood*Clutch.size*Clutches.yr)), group = first(group), fam_ = first(fam_)) %>% data.frame()


##Build phylogenies.
#This is the code I used to construct the phylogenies from the open tree of life in August of 2020. The phylogenies used in analyses can be found in the data folder on github. Use read.nexus(~/survival_tree_all_species_Nov_2020_3.nex) to read in the tree.
#tree_order<-read.nexus("~/Documents/survival_tree_all_species_Aug_2020 (1).nex")
#resolved_names_order <- tnrs_match_names(unique(surv_pgls_order$Species))
#tree_order <- tol_induced_subtree(ott_ids = resolved_names_order$ott_id)
#tree_order <- compute.brlen(multi2di(tree_order))
#write.nexus(tree_order, file = "survival_tree_all_species_Nov_2020_3.nex")

#Note - I remade this tree in Nov 2020. When I first wrote the tree out in August 2020 (with write.nexus), I did it after mapping the tip labels to the phylo object. This makes it difficult to read-in the phylogeny with read.nexus. So, I remade the tree and wrote it out in a form that is easy to read back in. 
tree_order <-read.nexus("~/Documents/research/completed/relative offspring size:survival/survival_tree_all_species_Nov_2020_3.nex")

#Map the tip labels to the data frame
taxon_map_order_ <- structure( unique(surv_pgls_order$animal), names =  unique(surv_pgls_order$Species))
otl_tips_order <- strip_ott_ids(tree_order$tip.label, remove_underscores = T)
otl_tips_order[otl_tips_order == "Gadus morhua (species in domain Eukaryota)"] <- "gadus morhua"
names(taxon_map_order_)[taxon_map_order_ == "gadus morhua"] <- "gadus morhua"
#cbind(sort(taxon_map_order_[otl_tips_order]), sort(otl_tips_order))#check
#cbind(taxon_map_order_, otl_tips_order ,tree_order$tip.label) #check
tree_order$tip.label <- taxon_map_order_[ otl_tips_order ]
#cbind(surv_pgls_order$animal, tree_order$tip.label)#check

invA_order <- inverseA(tree_order, nodes = "TIPS")


#plot(tree_order, type="fan")

###Code used to create phylogeny plot for supplement.
#fish_list <- filter(surv_pgls_order, group == "Fish")$animal
#mamm_list <- filter(surv_pgls_order, group == "Mammal")$animal
#rept_list <- filter(surv_pgls_order, group == "Reptile")$animal
#amph_list <- filter(surv_pgls_order, group == "Amphibian")$animal
#tree_order$node.label
#tree_order_plot <- paintSubTree(tree_order, findMRCA(tree_order, fish_list), state="Fish",  stem=T, ancestral.state="0")
#tree_order_plot <- paintSubTree(tree_order_plot, findMRCA(tree_order, mamm_list), state="Mammal",  stem=T, ancestral.state="1")
#tree_order_plot <- paintSubTree(tree_order_plot, findMRCA(tree_order, rept_list), state="Reptile", stem=T, ancestral.state="2")
#tree_order_plot <- paintSubTree(tree_order_plot, findMRCA(tree_order, amph_list), state="Amphibian",  stem=T, ancestral.state="3")
#cols <- setNames(c("black", "black","black", "black", "blue", "red", "brown", "green"), c("0", "1", "2", "3", "Fish", "Mammal", "Reptile", "Amphibian"))
#tree_order_plot$tip.label<-capitalize(tree_order_plot$tip.label)
#plot(tree_order_plot, type="fan", colors = cols, ftype="i", split.vertical=T, fsize=0.55)

############
#Code used to produce order-level phylogeny.
#resolved_names_order_lvl <- tnrs_match_names(unique(surv_pgls_order_lvl_data$animal), context_name = "Animals")
#tree_order_lvl <- tol_induced_subtree(ott_ids = resolved_names_order_lvl$ott_id)
#tree_order_lvl <- compute.brlen(multi2di(tree_order_lvl))
#plot(tree_order_lvl)
#cbind(sort(tree_order_lvl$tip.label), surv_pgls_order_lvl_data$animal)
#tree_order_lvl <- drop.tip(tree_order_lvl, 12:14)
#fish<-studies_find_studies(property="ot:focalCladeOTTTaxonName", value = "Actinopterygii")
#fish_phy <- get_study_tree(study_id = "ot_1699", tree_id="tree1")
#fish_fams <- tolower(filter(surv_pgls_order_lvl_data, group == "Fish")$fam_)
#fish_phy$tip.label <- str_split(tolower(fish_phy$tip.label), "_", simplify=T)[,1] 
#to_drop <- fish_phy$tip.label[!(fish_phy$tip.label %in% fish_fams)]
#fish_phy <- drop.tip(fish_phy, to_drop)
#fish_phy <- drop.tip(fish_phy, 1:8)
#fish_phy <- drop.tip(fish_phy, 2:4)
#fish_phy <- drop.tip(fish_phy,3:5)
#fish_phy <- drop.tip(fish_phy,5:7)
#fish_phy$tip.label <- c("Gadiformes", "Scombriformes", "Perciformes", "Pleuronectiformes", "Clupeiformes")
#tree_order_lvl <- addroot(tree_order_lvl, rootlength=1)
#fish_phy <- addroot(fish_phy, rootlength=1)
#tree_order_lvl <- bind.tree(fish_phy, tree_order_lvl, where="root")
#tree_order_lvl <- compute.brlen(multi2di(tree_order_lvl))
#plot(tree_order_lvl, show.node.label=F)
#tree_order_lvl$node.label <- seq(1:length(tree_order_lvl$node.label))
#write.nexus(tree_order_lvl, file="survival_tree_order_lvl_Aug_2020.nex")
tree_order_lvl<-read.nexus("~/Documents/research/completed/relative offspring size:survival/survival_tree_order_lvl_Aug_2020 (1).nex")

#~~# Note - I wrote the order-level tree out after mapping the tips to the names in the data frame. I still need to write some code to allow you to read the tree back in.

#taxon_map_order_lvl <- structure( resolved_names_order_lvl$search_string, names =  resolved_names_order_lvl$unique_name)
#otl_tips_order_lvl <- strip_ott_ids(tree_order_lvl$tip.label, remove_underscores = T)
#otl_tips_order_lvl[otl_tips_order_lvl == "Proboscidea (order in Deuterostomia)"] <- "Proboscidea"
#names(taxon_map_order_lvl)[taxon_map_order_lvl == "proboscidea"] <- "Proboscidea"
#otl_tips_order_lvl[otl_tips_order_lvl == "Squamata (order in Deuterostomia)"] <- "Squamate"
#names(taxon_map_order_lvl)[taxon_map_order_lvl == "squamate"] <- "Squamate"
#cbind(taxon_map_order_lvl[otl_tips_order_lvl], tree_order_lvl$tip.label) #check
#tree_order_lvl$tip.label <- taxon_map_order_lvl[ (otl_tips_order_lvl) ]
#cbind(sort(unique(surv_pgls_order_lvl_data$animal)), sort(tree_order_lvl$tip.label))
invA_order_lvl <- inverseA(tree_order_lvl, nodes = "TIPS")

#tree_order_lvl_plot <- tree_order_lvl
#tree_order_lvl_plot$tip.label <- capitalize(tree_order_lvl_plot$tip.label)
#plot(tree_order_lvl_plot)

##########################################
#First analysis - survivorship to relative offspring size

#Fit the models.

mcmc_mod1_order <- MCMCglmm(log10_survival ~ log10_mo_M *order, random = ~animal+Species, ginverse = list(animal=invA_order$Ainv), nitt = 100000, thin = 10, burnin = 5000, data = surv_pgls_order, pr=T, verbose=T)
summary(mcmc_mod1_order)
#plot(mcmc_mod1_order)



mcmc_mod1_order_lvl <- MCMCglmm(log10_survival ~ log10_mo_M ,random = ~animal, ginverse = list(animal=invA_order_lvl$Ainv),  data = surv_pgls_order_lvl_data, pr=T) 
summary(mcmc_mod1_order_lvl)
#plot(mcmc_mod1_order_lvl)

#OLS analyses
ols_mod1<-lm(log10_survival~log10_mo_M*order,data = surv_pgls_order)
summary(ols_mod1)

confint(lm(log10_survival~log10_mo_M,data = filter(surv_pgls_order, order=="Testudines")))
summary(lm(log10_survival~log10_mo_M,data = filter(surv_pgls_order, order=="Testudines")))

olsmod1order <- lm(log10_survival ~ log10_mo_M ,  data = surv_pgls_order)
summary(olsmod1order)
confint(olsmod1order)


#Extract slopes, intercepts, and 95%CIs.
head(surv_pgls)
slopes_1 <- cbind(mcmc_mod1_order$Sol[,2], mcmc_mod1_order$Sol[,2] + mcmc_mod1_order$Sol[,12:20]) %>% data.frame()
colnames(slopes_1) <- as.character(unique(sort(surv_pgls_order$order)))
slopes_1 <- tidyr::gather(slopes_1) %>% group_by(key) %>% summarize(mean_slope = mean(value), slope_low_ci = quantile(value, probs=c(0.025)), slope_high_ci = quantile(value, probs=c(0.975))) %>% data.frame()
head(slopes_1)
intercepts_1 <- cbind(mcmc_mod1_order$Sol[,1], mcmc_mod1_order$Sol[,1] + mcmc_mod1_order$Sol[,3:11]) %>% data.frame()
colnames(intercepts_1) <- as.character(unique(sort(surv_pgls_order$order)))
intercepts_1 <- tidyr::gather(intercepts_1) %>% group_by(key) %>% summarize(mean_intercept = mean(value), intercept_low_ci = quantile(value, probs=c(0.025)), intercept_high_ci = quantile(value, probs=c(0.975))) %>% data.frame()

all_params_1 <- inner_join(slopes_1, intercepts_1)
head(all_params_1)
all_params_1$group <- c("Amphibian", "Mammal", "Mammal", "Amphibian", "Fish", "Fish", "Mammal", "Mammal", "Reptile", "Reptile")
all_params_1$order <- as.factor(all_params_1$key)
all_params_1$order <- factor(all_params_1$order, levels = all_params_1$order[order(c(3,9,8,4,2,1,7,10,6,5))])


#Create data frame with data points adjusted for phylogeny.
fig1_data <- left_join(surv_pgls_order, all_params_1, by = "order") %>% mutate(predicted_log10_Sa = predict2(mcmc_mod1_order,  type = "lp", use = "mean")[1,1:nrow(surv_pgls_order)])  %>% mutate(residual_log10_Sa = log10_survival-predicted_log10_Sa) %>% rowwise() %>% mutate(corrected_log10_Sa = residual_log10_Sa + (mean_intercept + mean_slope*log10_mo_M))

fig1_data_order_lvl <- surv_pgls_order_lvl_data %>% mutate(predicted_log10_Sa = predict2(mcmc_mod1_order_lvl,  type = "lp", use = "mean")[1,1:nrow(surv_pgls_order_lvl_data)])  %>% mutate(residual_log10_Sa = log10_survival-predicted_log10_Sa) %>% rowwise() %>% mutate(corrected_log10_Sa = residual_log10_Sa + (summary(mcmc_mod1_order_lvl)$solutions[1,1] + summary(mcmc_mod1_order_lvl)$solutions[2,1]*log10_mo_M)) %>% mutate(dd = ifelse(order %in% unique(surv_pgls_order$order), "yes", "no"))


#Create plots.
all<-fig1_data %>% ggplot() + geom_abline(slope = 1, intercept = 0, linetype = 3) + geom_point(aes(x = log10_mo_M, y = corrected_log10_Sa, color = order, shape = group.x), size = 3, stroke = 1) + geom_smooth(aes(x = log10_mo_M, y = corrected_log10_Sa, linetype = order), size = 1.3,color = "black", method = "lm", se = F) + geom_smooth(aes(x = log10_mo_M, y = corrected_log10_Sa, color = order), size = 1.1, method = "lm", se = F)  + scale_shape_manual(values = c(21:24)) + scale_linetype_manual(values = rep(1,10)) + ylab(expression(paste("adjusted ",log[10]," survivorship (% surviving)")))+ xlab(expression(paste(log[10]," offspring mass/adult mass (g/g)"))) + xlim(c(-7, 0)) + ylim(c(-7, 0))  + scale_color_manual(values = c("green4", "red1", "red2", "green2", "blue4", "blue", "red3", "red4", "orange4", "orange")) + 
guides( color = guide_legend(override.aes = list(size = 0.5), title = NULL,  ncol=2), linetype = "none", shape = "none")  + theme(legend.position = c(0,0.9))  + scale_y_continuous(breaks = c(-6, -4, -2))+ xlab("") + ylab("")

inset <-all_params_1 %>% arrange(group) %>% ggplot() + geom_abline(slope = 0, intercept = 1, linetype = 3)  + geom_abline(slope = 0, intercept = 0) + geom_point(aes(x = order, y = mean_slope, color = order, shape = group), size = 3, stroke = .6) + geom_segment(aes(x = order, xend = order, y = slope_low_ci, yend = slope_high_ci))+ scale_y_continuous(limits = c(-3, 3)) + theme(legend.position = "none",  axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) +scale_shape_manual(values = c(21:24)) + scale_color_manual(values = c("blue", "blue4", "green4", "green2","orange", "orange4", "red3", "red2",   "red1", "red4" ))


means<-ggplot() + geom_abline(slope = 1, intercept = 0, linetype = 3) + geom_point(data = filter(fig1_data_order_lvl, dd == "no"), aes(x = log10_mo_M, y = corrected_log10_Sa, shape = group), size = 5, stroke = 1.2, color = "black") + geom_point(data = filter(fig1_data_order_lvl, dd == "yes"), aes(x = log10_mo_M, y = corrected_log10_Sa, shape = group, color = order), size = 5, stroke = 1.2)  + geom_smooth(data = fig1_data_order_lvl, aes(x = log10_mo_M, y = corrected_log10_Sa), method = "lm", se = F, color = "black")   + theme(legend.position = "none") + scale_shape_manual(values = c(21:24))  + xlim(c(-7, 0)) + ylab(expression(paste("mean ",log[10]," survivorship (% surviving)")))+ xlab(expression(paste("mean ",log[10]," offspring mass/adult mass (g/g)")))+ scale_color_manual(values = c("green4", "red1", "red2", "green2", "blue4", "blue", "red3", "red4", "orange4", "orange")) + guides(  linetype = guide_legend(title = NULL), shape = guide_legend(override.aes = list(size = 3), title = NULL,  ncol=1), fill = guide_legend(title = NULL), color = "none")  + 
theme(legend.position = c(0,0.9)) + scale_y_continuous(breaks = c(-6, -4, -2), limits = c(-7, 0))  + ylab("")+ xlab(expression(paste(log[10]," offspring mass/adult mass ", m[o]/m[a], " (g/g)")))


#Put plots together
p<-plot_grid(all, means, nrow=2, labels = "AUTO", label_x=0.09)
ggdraw(add_sub(p, expression(paste("adjusted ",log[10]," survivorship ", S[a], " (proportion surviving)")), y = -66, x = 0.04, angle = 90, vpadding=(grid::unit(-22.5, "lines")))) + draw_plot(inset, height = 0.8, width = 2,x = -0.25, y = 0.275, scale = 0.25) 



##########################################
#Second analysis - Fecundity per productivity vs offspring size


#First estimate relationships of productivity (rate of offspring biomass production) to adult size (Eqn 3 in main text)
mcmc_mod2_n <- MCMCglmm(log10_R ~ log10_M*order, random = ~animal+Species, ginverse = list(animal=invA_order$Ainv), data = surv_pgls_order,nitt = 100000, thin = 10, burnin = 5000, pr=T, verbose=F)
summary(mcmc_mod2_n)

mcmc_olsmod2_n <- lm(log10_R ~ log10_M*order,  data = surv_pgls_order)
summary(mcmc_olsmod2_n)

mcmc_mod2_n_order_lvl <- MCMCglmm(log10_R ~ log10_M, random = ~animal, ginverse = list(animal=invA_order_lvl$Ainv),prior=list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0))), data = surv_pgls_order_lvl_data, pr=T)
summary(mcmc_mod2_n_order_lvl)

parameters_n <- summary(mcmc_mod2_n)$solutions[,1]
parameters_n[3:11] <- parameters_n[3:11]+parameters_n[1]
parameters_n[12:length(parameters_n)] <- parameters_n[12:length(parameters_n)]+parameters_n[2]
parameters_n <- data.frame(order = unique(sort(surv_pgls_order$order)), intercepts = c(parameters_n[1], parameters_n[3:11]), slopes = c(parameters_n[2], parameters_n[12:length(parameters_n)]))


#Now, correct fecundity for differences in productivity using relationships of productivity to adult mass (as in Charnov & Ernest 2006 Am Nat)
fig3_data <- left_join(surv_pgls_order, parameters_n, by = "order") %>% mutate(predicted_R = predict2(mcmc_mod2_n,  type = "lp", use = "mean")[1, 1:nrow(surv_pgls_order)])  %>% mutate(residual_R = log10_R-predicted_R) %>% rowwise() %>% mutate(corrected_R = residual_R + (intercepts + slopes*log10_M), corrected_N = log10_N - (intercepts + slopes*log10_M) ) %>% data.frame()

fig3_data_order_lvl <- surv_pgls_order_lvl_data %>% mutate(predicted_R = predict2(mcmc_mod2_n_order_lvl,  type = "lp", use = "mean")[1, 1:nrow(surv_pgls_order_lvl_data)])  %>% mutate(residual_R = log10_R-predicted_R) %>% rowwise() %>% mutate(corrected_R = residual_R + (summary(mcmc_mod2_n_order_lvl)$solutions[1,1] + summary(mcmc_mod2_n_order_lvl)$solutions[2,1]*log10_M), corrected_N = log10_N - (summary(mcmc_mod2_n_order_lvl)$solutions[1,1] + summary(mcmc_mod2_n_order_lvl)$solutions[2,1]*log10_M) ) %>% data.frame()

#Finally, estimate relationships of fecundity to offspring mass, with fecundity corrected for differences in productivity.
mcmc_mod3_n <- MCMCglmm(corrected_N ~  log10_mo * order, random = ~animal+Species, ginverse = list(animal= invA_order$Ainv),nitt = 100000, thin = 10, burnin = 5000,  data = fig3_data, pr=T, verbose=F)
summary(mcmc_mod3_n)

mcmc_mod3_n_order_lvl <- MCMCglmm(corrected_N ~  log10_mo , random = ~animal, ginverse = list(animal= invA_order_lvl$Ainv),prior=list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0))),  data = fig3_data_order_lvl, pr=T)
summary(mcmc_mod3_n_order_lvl)

#Extract parameters from mcmc object to construct plots.
slopes_nc <- cbind(mcmc_mod3_n$Sol[,2], mcmc_mod3_n$Sol[,2] + mcmc_mod3_n$Sol[,12:20]) %>% data.frame()
colnames(slopes_nc) <- as.character(unique(sort(surv_pgls_order$order)))
slopes_nc <- tidyr::gather(slopes_nc) %>% group_by(key) %>% summarize(mean_slope = mean(value), slope_low_ci = quantile(value, probs=c(0.025)), slope_high_ci = quantile(value, probs=c(0.975))) %>% data.frame()
head(slopes_nc)
intercepts_nc <- cbind(mcmc_mod3_n$Sol[,1], mcmc_mod3_n$Sol[,1] + mcmc_mod3_n$Sol[,3:11]) %>% data.frame()
colnames(intercepts_nc) <- as.character(unique(sort(surv_pgls_order$order)))
intercepts_nc <- tidyr::gather(intercepts_nc) %>% group_by(key) %>% summarize(mean_intercept = mean(value), intercept_low_ci = quantile(value, probs=c(0.025)), intercept_high_ci = quantile(value, probs=c(0.975))) %>% data.frame()
head(slopes_C)
all_params_nc <- inner_join(slopes_nc, intercepts_nc)
head(all_params_nc)
all_params_nc$group <- c("Amphibian", "Mammal", "Mammal", "Amphibian", "Fish", "Fish", "Mammal", "Mammal", "Reptile", "Reptile")
all_params_nc$order <- as.factor(all_params_nc$key)
all_params_nc$order <- factor(all_params_nc$order, levels = all_params_nc$order[order(c(3,9,8,4,2,1,7,10,6,5))])



fig3_data <- left_join(fig3_data, all_params_nc, by = "order") %>% mutate(predicted_Nc = predict2(mcmc_mod3_n,  type = "lp", use = "mean")[1,1:nrow(surv_pgls_order)])  %>% mutate(residual_Nc = corrected_N-predicted_Nc) %>% rowwise() %>% mutate( corrected_N_phy = residual_Nc + (mean_intercept + mean_slope*log10_mo) ) %>% data.frame()


fig3_data_order_lvl <- fig3_data_order_lvl %>% mutate(predicted_Nc = predict2(mcmc_mod3_n_order_lvl,  type = "lp", use = "mean")[1,1:nrow(fig3_data_order_lvl)])  %>% mutate(residual_Nc = corrected_N-predicted_Nc) %>% rowwise() %>% mutate( corrected_N_phy = residual_Nc + (summary(mcmc_mod3_n_order_lvl)$solutions[1,1] + summary(mcmc_mod3_n_order_lvl)$solutions[2,1]*log10_mo) ) %>% data.frame()%>% mutate(dd = ifelse(order %in% unique(surv_pgls_order$order), "yes", "no"))



#Construct plots.
fig3b<-
fig3_data %>% ggplot() + geom_abline(slope = -1, intercept = 0, linetype = 3) + geom_point(aes(x = log10_mo, y = corrected_N_phy, color = order, shape = group.x), size = 3, stroke = 1) + geom_smooth(aes(x = log10_mo, y = corrected_N_phy, linetype = order), size = 1.3,color = "black", method = "lm", se = F)+ geom_smooth(aes(x = log10_mo, y = corrected_N_phy, linetype = order, color = order), size = 1, method = "lm", se = F)  + scale_shape_manual(values = c(21:24)) + scale_linetype_manual(values = rep(1,10))  +xlab("")+ ylab(expression(paste(log[10]," fecundity per productivity (#/yr/", g^{b}, ")")))   + scale_color_manual(values = c("green4", "red1", "red2", "green2", "blue4", "blue", "red3", "red4", "orange4", "orange")) + guides( color = guide_legend(override.aes = list(size = 0.5), title = NULL,  ncol=2), linetype ="none", shape =  "none")   + theme(legend.position = c(0.25,0.9)) + scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5)) + ylab("")


inset_fig3 <-
all_params_nc %>% arrange(group) %>% ggplot() + geom_abline(slope = 0, intercept = -1, linetype = 3)  + geom_abline(slope = 0, intercept = 0) + geom_point(aes(x = order, y = mean_slope, color = order, shape = group), size = 3, stroke = .6) + geom_segment(aes(x = order, xend = order, y = slope_low_ci, yend = slope_high_ci))+ scale_y_continuous(limits = c(-3, 3)) + theme(legend.position = "none",  axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) +scale_shape_manual(values = c(21:24)) + scale_color_manual(values = c("blue", "blue4", "green4", "green2","orange", "orange4", "red3", "red2",   "red1", "red4" )) + ylab("slopes")





fig3c<-
 ggplot() + geom_abline(slope = -1, intercept = 0, linetype = 3) + geom_point(data= filter(fig3_data_order_lvl, dd=="no"), aes(x = log10_mo, y = corrected_N_phy, shape = group), color ="black", size = 4, stroke=1.1)+ geom_point(data= filter(fig3_data_order_lvl, dd=="yes"), aes(x = log10_mo, y = corrected_N_phy, color = order, shape = group), size = 4, stroke=1.1) + geom_smooth(data= fig3_data_order_lvl, aes(x = log10_mo, y = corrected_N_phy), color = "black",  method = "lm", se=F) +  theme(legend.position = "none") + scale_shape_manual(values = c(21:24))  + ylab(expression(paste("mean ",log[10]," fecundity per productivity (#/yr/", g^{b}, ")")))+ xlab(expression(paste("mean ",log[10]," offspring mass (g)"))) + scale_color_manual(values = c("green4", "red1", "red2", "green2", "blue4", "blue", "red3", "red4", "orange4", "orange")) + guides( color = "none", linetype = guide_legend(title = NULL), fill = guide_legend(title = NULL), shape = guide_legend(override.aes = list(size = 3), title = NULL,  ncol=1))  + theme(legend.position = c(0.625,0.9)) + scale_y_continuous(breaks = c(-5, -2.5, 0, 2.5))+xlab(expression(paste(log[10]," offspring mass ", m[o], " (g)"))) + ylab("")


p_nc<-plot_grid(fig3b, fig3c, labels = "AUTO", nrow = 2, label_x = 0.12)
ggdraw(add_sub(p_nc, expression(paste("adjusted ",log[10]," fecundity per productivity ", N/I(m[a])^{x}, " (#/yr)")), y = -65, x = 0.04, angle = 90, vpadding=(grid::unit(-24, "lines")))) + draw_plot(inset_fig3, height = 0.8, width = 2,x = -0.57, y = 0.262, scale = 0.25)




#OLS analyses for supplement

ols_params_2 <- data.frame(order=sort(unique(surv_pgls_order$order)), slope = rep(0,10), intercept=rep(0,10))

for(i in 1:10){
	
	ols_params_2$slope[i] <- summary(lm(log10_R ~ log10_M, data = filter(surv_pgls_order, order == ols_params_2$order[i])))$coefficients[2,1]
	ols_params_2$intercept[i] <- summary(lm(log10_R ~ log10_M, data = filter(surv_pgls_order, order == ols_params_2$order[i])))$coefficients[1,1]
		
}

ols_corrected_n_data <- left_join(ols_params_2, surv_pgls_order, by="order") %>% mutate(corrected_N = log10_N - (intercept + slope*log10_M))

ols_corrected_2 <- data.frame(order=sort(unique(surv_pgls_order$order)), slope = rep(0,10), intercept=rep(0,10), slope_low_ci = rep(0,10), slope_high_ci = rep(0,10), intercept_low_ci = rep(0,10), intercept_high_ci = rep(0,10))

for(i in 1:10){
	
	ols_corrected_2$slope[i] <- summary(lm(corrected_N ~ log10_mo, data = filter(ols_corrected_n_data, order == ols_corrected_2$order[i])))$coefficients[2,1]
	ols_corrected_2$intercept[i] <- summary(lm(corrected_N ~ log10_mo, data = filter(ols_corrected_n_data, order == ols_corrected_2$order[i])))$coefficients[1,1]
		
	ols_corrected_2$slope_low_ci[i] <- confint(lm(corrected_N ~ log10_mo, data = filter(ols_corrected_n_data, order == ols_corrected_2$order[i])))[2,1]
	ols_corrected_2$slope_high_ci[i] <- confint(lm(corrected_N ~ log10_mo, data = filter(ols_corrected_n_data, order == ols_corrected_2$order[i])))[2,2]
	
	ols_corrected_2$intercept_low_ci[i] <- confint(lm(corrected_N ~ log10_mo, data = filter(ols_corrected_n_data, order == ols_corrected_2$order[i])))[1,1]
	ols_corrected_2$intercept_high_ci[i] <- confint(lm(corrected_N ~ log10_mo, data = filter(ols_corrected_n_data, order == ols_corrected_2$order[i])))[1,2]
}

lm(log10_R ~  log10_M ,   data = surv_pgls_order_lvl_data)
lm(corrected_N~log10_mo,data= mutate(surv_pgls_order_lvl_data, corrected_N=log10_N-(-0.3735+0.9022*log10_M)))
confint(lm(corrected_N~log10_mo,data= mutate(surv_pgls_order_lvl_data, corrected_N=log10_N-(-0.3735+0.9022*log10_M))))

##########################################
#Third (last) analysis - reproductive success C vs offspring size

head(surv_pgls_order)
surv_pgls_order <- surv_pgls_order %>% mutate(log10_C = log10((10^log10_survival)*(10^log10_N)))
mcmc_mod3_order <- MCMCglmm(log10_C ~ log10_mo *order, random = ~animal+Species, ginverse = list(animal=invA_order$Ainv), prior=list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0), G2=list(V=1,nu=0))), nitt = 100000, thin = 10, burnin = 5000, data = surv_pgls_order, pr=T)
summary(mcmc_mod3_order)
plot(mcmc_mod3_order)

confint(lm(log10_C ~ log10_mo,data = filter(surv_pgls_order, order=="Anura")))
summary(lm(log10_C ~ log10_mo,data = filter(surv_pgls_order, order=="Anura")))

mcmc_mod3_order_lvl <- MCMCglmm(log10_C ~ log10_mo ,random = ~animal, ginverse = list(animal=invA_order_lvl$Ainv),  data = surv_pgls_order_lvl_data, prior=list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0))),nitt = 100000, thin = 10, burnin = 5000, pr=T)
summary(mcmc_mod3_order_lvl)
#plot(mcmc_mod3_order_lvl)

olsmod3order <- lm(log10_C ~ log10_mo ,  data = surv_pgls_order_lvl_data)
summary(olsmod3order)
confint(olsmod3order)

###Here we show the relationship is not significant even if one removes the Orders with low sample sizes:
#tree_order_lvl_test <- drop.tip(tree_order_lvl_plot, c("proboscidea", "perissodactyla", "primate", "pleuronectiformes", "perciformes", "scombriformes"))
#plot(tree_order_lvl_test)
#invA_order_lvl_test <- inverseA(tree_order_lvl_test, nodes = "TIPS")
#surv_pgls_order_lvl_data_test <- filter(surv_pgls_order_lvl_data, !(order %in% c("Proboscidea", "Perissodactyla", "Primate", "Pleuronectiformes", "Perciformes", "Scombriformes")))

#mcmc_mod3_order_lvl_test <- MCMCglmm(log10_C ~ log10_mo ,random = ~animal, ginverse = list(animal= invA_order_lvl_test$Ainv),  data = surv_pgls_order_lvl_data_test, prior = list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002))), pr=T)
#summary(mcmc_mod3_order_lvl_test)


#Extract relevant parameters for plotting.
slopes_C <- cbind(mcmc_mod3_order$Sol[,2], mcmc_mod3_order$Sol[,2] + mcmc_mod3_order$Sol[,12:20]) %>% data.frame()
colnames(slopes_C) <- as.character(unique(sort(surv_pgls_order$order)))
slopes_C <- tidyr::gather(slopes_C) %>% group_by(key) %>% summarize(mean_slope = mean(value), slope_low_ci = quantile(value, probs=c(0.025)), slope_high_ci = quantile(value, probs=c(0.975))) %>% data.frame()
head(slopes_C)
intercepts_C <- cbind(mcmc_mod3_order$Sol[,1], mcmc_mod3_order$Sol[,1] + mcmc_mod3_order$Sol[,3:11]) %>% data.frame()
colnames(intercepts_C) <- as.character(unique(sort(surv_pgls_order$order)))
intercepts_C <- tidyr::gather(intercepts_C) %>% group_by(key) %>% summarize(mean_intercept = mean(value), intercept_low_ci = quantile(value, probs=c(0.025)), intercept_high_ci = quantile(value, probs=c(0.975))) %>% data.frame()
head(slopes_C)
all_params_C <- inner_join(slopes_C, intercepts_C)
head(all_params_C)
all_params_C$group <- c("Amphibian", "Mammal", "Mammal", "Amphibian", "Fish", "Fish", "Mammal", "Mammal", "Reptile", "Reptile")
all_params_C$order <- as.factor(all_params_C$key)
all_params_C$order <- factor(all_params_C$order, levels = all_params_C$order[order(c(3,9,8,4,2,1,7,10,6,5))])


head(fig4_data)
fig4_data <- left_join(surv_pgls_order, all_params_C, by = "order") %>% mutate(predicted_log10_C = predict2(mcmc_mod3_order,  type = "lp", use = "mean")[1,1:nrow(surv_pgls_order)])  %>% mutate(residual_log10_C = log10_C-predicted_log10_C) %>% rowwise() %>% mutate(corrected_log10_C = residual_log10_C + (mean_intercept + mean_slope*log10_mo)) %>% data.frame()

fig4_data_order_lvl <- surv_pgls_order_lvl_data %>% mutate(predicted_log10_C = predict2(mcmc_mod3_order_lvl,  type = "lp", use = "mean")[1,1:nrow(surv_pgls_order_lvl_data)])  %>% mutate(residual_log10_C = log10_C-predicted_log10_C) %>% rowwise() %>% mutate(corrected_log10_C = residual_log10_C + (summary(mcmc_mod3_order_lvl)$solutions[1,1] + summary(mcmc_mod3_order_lvl)$solutions[2,1]*log10_mo)) %>% mutate(dd = ifelse(order %in% unique(surv_pgls_order$order), "yes", "no"))

fig4a<-
fig4_data %>% ggplot()  + geom_point(aes(x = log10_mo, y = corrected_log10_C, color = order, shape = group.x), size = 3, stroke = 1) + geom_smooth(aes(x = log10_mo, y = corrected_log10_C, linetype = order), size = 1.3,color = "black", method = "lm", se = F)+ geom_smooth(aes(x = log10_mo, y = corrected_log10_C, linetype = order, color = order), size = 1, method = "lm", se = F)  + scale_shape_manual(values = c(21:24)) + scale_linetype_manual(values = rep(1,10))  +xlab("")   + scale_color_manual(values = c("green4", "red1", "red2", "green2", "blue4", "blue", "red3", "red4", "orange4", "orange")) + guides( color = guide_legend(override.aes = list(size = 0.5), title = NULL,  ncol=2), linetype ="none", shape =  "none")  + theme(legend.position = c(0,0.9)) + scale_y_continuous(limits = c(-2.5, 2), breaks = c(-2, -1, 0, 1)) + scale_x_continuous(limits = c(-5.5, 6)) + ylab("")

#ggplot()  + geom_point(data=filter(fig4_data, order=="Caudata"),aes(x = log10_mo, y = corrected_log10_C, color = order, shape = group), size = 3, stroke = 1) + geom_smooth(data=filter(fig4_data, order=="Caudata"),aes(x = log10_mo, y = corrected_log10_C, linetype = order), size = 1.3,color = "black", method = "lm", se = F)+ geom_smooth(data=filter(fig4_data, order=="Caudata"),aes(x = log10_mo, y = corrected_log10_C, linetype = order, color = order), size = 1, method = "lm", se = F)  +
#geom_point(data=filter(fig4_data, order!="Caudata"),aes(x = log10_mo, y = log10_C, color = order, shape = group), size = 3, stroke = 1) + geom_smooth(data=filter(fig4_data, order!="Caudata"),aes(x = log10_mo, y = log10_C, linetype = order), size = 1.3,color = "black", method = "lm", se = F)+ geom_smooth(data=filter(fig4_data, order!="Caudata"),aes(x = log10_mo, y = log10_C, linetype = order, color = order), size = 1, method = "lm", se = F)  +
 #scale_shape_manual(values = c(21:24)) + scale_linetype_manual(values = rep(1,10))  +xlab("")   + scale_color_manual(values = c("green4", "red1", "red2", "green2", "blue4", "blue", "red3", "red4", "orange4", "orange")) + guides( color = guide_legend(override.aes = list(size = 0.5), title = NULL,  ncol=2), linetype ="none", shape =  "none")  + theme(legend.position = c(0,0.9)) + scale_y_continuous(limits = c(-2.5, 2), breaks = c(-2, -1, 0, 1)) + scale_x_continuous(limits = c(-5.5, 6)) + ylab("")



inset_fig4 <-
all_params_C %>% arrange(group) %>% ggplot()  + geom_abline(slope = 0, intercept = 0) + geom_point(aes(x = order, y = mean_slope, color = order, shape = group), size = 3, stroke = .6) + geom_segment(aes(x = order, xend = order, y = slope_low_ci, yend = slope_high_ci))+ scale_y_continuous(limits = c(-1.5, 1.5)) + theme(legend.position = "none",  axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10)) +scale_shape_manual(values = c(21:24)) + scale_color_manual(values = c("blue", "blue4", "green4", "green2","orange", "orange4", "red3", "red2",   "red1", "red4" )) + ylab("slopes")




fig4b<-
 ggplot()  + geom_point(data= filter(fig4_data_order_lvl, dd=="no"), aes(x = log10_mo, y = corrected_log10_C, shape = group), color ="black", size = 4, stroke=1.1)+ geom_point(data= filter(fig4_data_order_lvl, dd=="yes"), aes(x = log10_mo, y = corrected_log10_C, color = order, shape = group), size = 4, stroke=1.1) + geom_smooth(data= fig4_data_order_lvl, aes(x = log10_mo, y = corrected_log10_C), color = "black",  method = "lm", se=F) +  theme(legend.position = "none") + scale_shape_manual(values = c(21:24))  + ylab(expression(paste("mean ",log[10]," fecundity per productivity (#/yr/", g^{b}, ")")))+ xlab(expression(paste("mean ",log[10]," offspring mass (g)"))) + scale_color_manual(values = c("green4", "red1", "red2", "green2", "blue4", "blue", "red3", "red4", "orange4", "orange")) + guides( color = "none", linetype = guide_legend(title = NULL), fill = guide_legend(title = NULL), shape = guide_legend(override.aes = list(size = 3), title = NULL,  ncol=1))  + theme(legend.position = c(0,0.9)) + xlab(expression(paste(log[10]," offspring mass ", m[o], " (g)"))) + ylab("") + scale_y_continuous(limits = c(-2, 1), breaks = c(-2, -1, 0))+ scale_x_continuous(limits = c(-5.5, 6))


p_4<-plot_grid(fig4a, fig4b, nrow=2, labels = "AUTO", label_x=0.09)

ggdraw(add_sub(p_4, expression(paste("adjusted ",log[10]," reproductive success C (# surviving/yr)")), y = 1200, x = 0.04, angle = 90, vpadding=(grid::unit(-23, "lines")))) + draw_plot(inset_fig4, height = 0.8, width = 2,x = -0.2, y = 0.25, scale = 0.2) 


