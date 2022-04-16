#clears
rm(list = ls())
dev.off()

#packages
install.packages(c("tidyverse","GGally","lmtest","moments","forecast","nlme","devtools","factoextra","pls"))
my_packages <- c("tidyverse","GGally","lmtest","moments","forecast","nlme","devtools","factoextra","pls")
lapply(my_packages, require, character.only = TRUE)

#code
decomp <- read.csv("Inputs/decomp_data.csv")
decomp_sub <- decomp[,-c(1:3)]
ggpairs(decomp_sub[,-8])

decomp_pca <- prcomp(decomp_sub[,-8], center = TRUE, scale. = TRUE)
print(decomp_pca)
eig.val <- get_eigenvalue(decomp_pca)
eig.val               
fviz_screeplot(decomp_pca, addlabels = TRUE)

# PCA with factoextra -----------------------------------------------------
humus_type <- as.factor(decomp_sub$humus_type)
fviz_pca_biplot(decomp_pca,
                col.ind = humus_type,
                palette = c("gold2","cyan3","forestgreen"),
                mean.point = FALSE,
                addEllipses = TRUE,
                ellipse.type = "euclid",
                legend.title = "Humus types",
                repel = TRUE,
                title = " ")

fviz_contrib(decomp_pca, 
             choice = "var", 
             fill = "cornflowerblue",
             color = "cornflowerblue",
             axes = 1, 
             top = 10)

fviz_contrib(decomp_pca, 
             choice = "var",
             fill = "cornflowerblue",
             color = "cornflowerblue",
             axes = 2, 
             top = 10)

fviz_contrib(decomp_pca, 
             choice = "var", 
             fill = "cornflowerblue",
             color = "cornflowerblue",
             axes = 1:2, 
             top = 10)

# PCA with ggbiplot -------------------------------------------------------
install_github("vqv/ggbiplot", force = TRUE)
library(ggbiplot)
decomp_biplot <- ggbiplot(decomp_pca,
                          labels = rownames(decomp),
                          obs.scale = 1,
                          var.scale = 1,
                          varname.size = 4,
                          varname.adjust = 2,
                          groups = decomp_sub$humus,
                          ellipse = TRUE) +
  scale_colour_manual(name = "Humus types", values = c("gold2","cyan3","forestgreen")) +
  theme(legend.position = "right", legend.direction = "vertical", legend.key = )
print(decomp_biplot)

# PCR with pls  -----------------------------------------------------------
decomp_pcr <- pcr(litter_loss ~., data = decomp_sub[,-8], scale = TRUE, validation = "CV")
summary(decomp_pcr)
validationplot(decomp_pcr)

train <- decomp_sub[1:75,]
y_test <- decomp_sub[75:112, 1]
test <- decomp_sub[75:112, 2:8]

Decomp_pcr <- pcr(litter_loss ~., data = train, scale =TRUE, validation = "CV")

pcr_pred <- predict(Decomp_pcr, test, ncomp = 4)
mean((pcr_pred - y_test)^2)

# MLR with specific transformations ---------------------------------------
hist(decomp_sub[,1]) #variables 2 and 3 positively skewed, variable 5 super positively skewed 

skewness(decomp_sub$bacterivore_abund, na.rm = TRUE) #value of 1.6904 > 0
skewness(decomp_sub$fungivore_abund, na.rm = TRUE) #value of 1.2021 > 0
skewness(decomp_sub$Othick, na.rm = TRUE) #value of 1.3254 > 0

hist(sqrt(decomp_sub[,2])) #sqrt transformations adequate for variables 2 and 3
zeros <- colSums(decomp_sub == 0) #log + 1 transformation for variable 5 since log(0) = undef

hist(decomp_sub[,5])
hist(decomp_sub_lm[,8])

sqrt_bacteri <- sqrt(decomp_sub$bacterivore_abund)
decomp_sub_lm <- decomp_sub %>%
  add_column(bacteri_sqrt = sqrt_bacteri, .after = "bacterivore_abund")

sqrt_fungi <- sqrt(decomp_sub$fungivore_abund)
decomp_sub_lm <- decomp_sub_lm %>%
  add_column(fungi_sqrt = sqrt_fungi, .after = "fungivore_abund")

log_Othick <- log(decomp_sub$Othick + 1)
decomp_sub_lm <- decomp_sub_lm %>%
  add_column(Othick_log = log_Othick, .after = "Othick")

decomp_fit <-lm(litter_loss ~ sqrt_bacteri + sqrt_fungi + tree_div + Athick + Bthick + log_Othick, data = decomp_sub_lm)
summary(decomp_fit) #non-normal & homoscedastic; adj-R2 0.03871

decomp_fit <-lm(log_Othick ~ sqrt_bacteri + sqrt_fungi+ tree_div + Athick + Bthick + litter_loss, data = decomp_sub_lm)
summary(decomp_fit) #non-normal & heteroscedastic; adj-R2 0.3765

plot(decomp_fit, which = 2)
shapiro.test(residuals(decomp_fit))
plot(decomp_fit, which = 3)
bptest(decomp_fit)

decomp_gls <- gls(log_Othick ~ sqrt_bacteri*litter_loss + sqrt_fungi*litter_loss + Athick + litter_loss, weights = varPower(), data = decomp_sub_lm)
plot(decomp_gls)
summary(decomp_gls)

# MRL with Box-Cox transformation --------------------------------------------------
decomp_sub <- decomp[,-c(1:3,11)]

BoxCox.lambda(decomp_sub$bacterivore_abund)
bacti_BoxCox <- BoxCox(decomp_sub$bacterivore_abund, lambda = 1.005365)
BoxCox.lambda(decomp_sub$fungivore_abund)
fungi_BoxCox <- BoxCox(decomp_sub$bacterivore_abund, lambda = 0.08341178)
BoxCox.lambda(decomp_sub$Othick + 1)
Othick_BoxCox <- BoxCox(decomp_sub$Othick +1, lambda = -0.09093155)

decomp_fit <-lm(litter_loss ~ bacti_BoxCox + fungi_BoxCox + tree_div + Athick + Bthick + Othick_BoxCox, data = decomp_sub)
summary(decomp_fit) #normal & homoscedastic; adj-R2 0.04887

decomp_fit <-lm(Othick_BoxCox ~ bacti_BoxCox + fungi_BoxCox + tree_div + Athick + Bthick + litter_loss, data = decomp_sub)
summary(decomp_fit) #non-normal & heteroscedastic visually (by test opposite outcomes); adj-R2 0.3804



# MLR mixed transformations -----------------------------------------------
decomp_sub <- decomp[,-c(1:3,11)]

sqrt_bacteri <- sqrt(decomp_sub$bacterivore_abund)
decomp_sub_lm <- decomp_sub %>%
  add_column(bacteri_sqrt = sqrt_bacteri, .after = "bacterivore_abund")

sqrt_fungi <- sqrt(decomp_sub$fungivore_abund)
decomp_sub_lm <- decomp_sub_lm %>%
  add_column(fungi_sqrt = sqrt_fungi, .after = "fungivore_abund")

BoxCox.lambda(decomp_sub$Othick + 1)
Othick_BoxCox <- BoxCox(decomp_sub_lm$Othick +1, lambda = -0.09093155)

decomp_fit <-lm(litter_loss ~ sqrt_bacteri + sqrt_fungi + tree_div + Athick + Bthick + Othick_BoxCox, data = decomp_sub_lm)
summary(decomp_fit) #normal & homoscedastic; adj-R2 0.03767

decomp_fit <-lm(Othick_BoxCox ~ sqrt_bacteri + sqrt_fungi + tree_div + Athick + Bthick + litter_loss, data = decomp_sub_lm)
summary(decomp_fit) #non-normal & heteroscedastic visually (by test opposite outcomes); adj-R2 0.3783

plot(decomp_fit, which = 2)
shapiro.test(residuals(decomp_fit))
plot(decomp_fit, which = 3)
bptest(decomp_fit)


