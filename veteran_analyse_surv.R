## packages

library(survival)   ## Ajustement des modèles de durées de survie
library(survminer)  ## Plot des fonctions de survie
library(ggplot2)    ## Graphiques
library(knitr)       # Rendre les tableaux plus présentables
library(kableExtra)   #Rendre les tableaux plus présentables
library(skimr)
library(tidyverse)  ## Manipulation de jeu de données

#chargement du jeu de données
vet <- veteran
dim(vet)
head(vet)

## Valeurs manquantes ?
any(is.na(vet))


# Nature des variables

str(vet)

## Conversion en factor (variables catégorielles/nominales)
vet$trt <- factor(vet$trt, levels = c(1, 2), labels = c("Standard", "test drug"))
vet$prior <- factor(vet$prior)

# Résumé descriptif
summary(vet)

# analyse descriptive des variables continues
base_stat <- skim(vet)[c(4,6:8),
                       c("skim_variable","numeric.p0","numeric.p100",
                         "numeric.mean", "numeric.sd", "numeric.p25", 
                         "numeric.p50", "numeric.p75")]

tabstat_names <- c("Variables","Min","Max","Moyenne",
                   "Ecart type", "Quartile 1","Mediane","Quartile 3")

names(base_stat) <- tabstat_names

base_stat <- base_stat |> mutate(across(where(is.numeric),\(x) round(x,2)))

## rendu du tableau descriptif
kable(base_stat,caption = "Quelques statistiques")|>
  kable_styling(latex_options = c("striped","scale_down","repeat_header"))|>
  column_spec(1, bold = T,width = "3cm", latex_valign = "m") |>
  row_spec(0, bold = T, color = "white", background = "black")


## Proportion de données censurées
prop.table(table(factor(vet$status)))


## Boxplot du temps de survie en fonction du traitement
survie_trt_plot <- vet|>
  filter(status==1)|>
  ggplot(aes(x=trt, y = time, fill =trt))+
  geom_boxplot()+
  labs(
    y = "Temps de survie",
    x = "Type de traitement",
    title = ""
  )+
  theme(
    legend.position = "bottom"
  )
survie_trt_plot


## temps de survie en fonction du score de Karnofsky
survie_score_plot <- vet|>filter(status==1)|>ggplot(aes(x=as.factor(karno), y = time))+
  geom_boxplot()+
  labs(
    y = "Temps de survie",
    x = "score de karno",
    title = ""
  )
survie_score_plot


## Ajustement des modèles

### Modèle de Kaplan-Meier (Survie globale)
km_model <- survfit(Surv(time, status)~1, data = vet)
summary(km_model)

# Temps médian de survie(50%)
median_time <- as.numeric(summary(km_model)$table["median"])

## Courbe de kaplan Meier 
km_curve <- ggsurvplot(
  fit = km_model,
  data=vet,
  conf.int = T,
  surv.median.line = "hv",
  risk.table = T,
  break.time.by = median_time,
  ggtheme = theme_minimal(),
  
  title = "Courbe survie de Kaplan Meier",
  legend.position = "none"
  
)

km_curve


##Survie suivant le traitement
km_model_trt <- survfit(Surv(time, status)~trt, data = vet)
summary(km_model_trt)

## Temps median
mt_trt <- as.numeric(summary(km_model_trt)$table[,"median"])
mt_trt

## courbe de survie
km_curve_trt <- ggsurvplot(
  fit = km_model_trt,
  data=vet,
  conf.int = T,
  surv.median.line = "hv",
  risk.table = T,
  ggtheme = theme_minimal(),
  title = "Courbe de survie de Kaplan Meier",
  pval = T
  
)

km_curve_trt

## Test de log rank 
survdiff(Surv(time, status)~trt, data = vet)


## Modèle de cox 
model_cox <- coxph(
  Surv(time, status)~., data = vet)

broom::tidy(model_cox) |>
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    p.value = signif(p.value, 3)
  )

## effet global des facteurs catégoriels
anova(model_cox)


## Diagnostic du modèle de cox
### Hypothèse de risque proportionnel
prop_test <- cox.zph(model_cox)
prop_test

## Courbe des résidus de Schoenfeld
ggcoxzph(cox.zph(model_cox), ggtheme=theme())


## Récupération des beta(t) pour la variable karno afin d'illustrer 
## la variation du risque relatif dans le temps


beta_t_df <- data.frame(
  time = prop_test$time,       
  beta_t = prop_test$y[, "karno"]
)

## Risque relatif pour Q1 et Q3 (Deux groupes différents pour servir de comparaison des HR)
Q1 <- 40
Q3 <- 75
beta_t_df$HR_Q1_Q3 <- exp(beta_t_df$beta_t * (Q1 - Q3))

## Courbe 

ggplot(beta_t_df, aes(x = time, y = HR_Q1_Q3)) +
  geom_line(color = "blue") +
  labs(x = "Temps", y = "HR(Q1 vs Q3)", 
       title = "Evolution du hazard ratio pour Karno") +
  theme_minimal()



## Modèle Stratifié(celltype)

model_strat <- coxph(Surv(time, status)~ karno + strata(celltype) , data = vet)

cox.zph(model_strat)
ggcoxzph(cox.zph(model_strat), ggtheme=theme())


## Modélisation de la non proportionnalité

## Découpage des observations pour les différents intervalles
vet_new <- survSplit(Surv(time, status)~., data = vet, 
                     cut=c(90,200), 
                     episode = "timegroup")


## Différents modèles testés
model_cox2 <- coxph(Surv(tstart, time, status)~trt+prior + diagtime + 
                      age+ karno*strata(timegroup)+celltype, data = vet_new )

model_cox3 <- coxph(Surv(tstart, time, status)~ karno*strata(timegroup)+
                      celltype*strata(timegroup), data = vet_new )

model_cox4 <- coxph(Surv(tstart, time, status)~ celltype + 
                      karno*strata(timegroup), data = vet_new )

model_cox5 <- coxph(Surv(tstart, time, status)~ trt + 
                      karno*strata(timegroup), data = vet_new )


## choix du meilleur modèle

AIC(model_cox2)

AIC(model_cox3)

AIC(model_cox4)

AIC(model_cox5)


## Validation modèle final
cox.zph(model_cox4) ## proportionnalité ok

ggcoxdiagnostics(model_cox4, type="martingale", linear.predictions = F) ## linéarité ok

ggcoxdiagnostics(model_cox4, type="deviance", linear.predictions = F)

## coefficients estimés 
broom::tidy(model_cox4, exponentiate = T )|>
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 3),
    p.value = signif(p.value, 3)
  )

## Effet du score de Karnofsky sur chaque période
b <- coef(model_cox4)
inter <- b[grep("karno:", names(b))]
hr_periode <- exp(b["karno"] + c(0, inter))
names(hr_periode) <- c("0-90 jours", "90-200 jours", "> 200 jours")
hr_periode ## HR pour 1 point
hr_periode^10 ## HR pour 10 points

## Qualité du modèle

broom::glance(model_cox4)|>
  pivot_longer(cols=-n)|>
  select(name, value)|>
  slice(3,7,12,15)

