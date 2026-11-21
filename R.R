cancercell<-read.csv(choose.files())
dim(cancercell)
names(cancercell)
summary(cancercell)
str(cancercell)
analysis<-cancercell%>%
  select(
    rownames,
    trt,
    celltype,
    time,
    status,
    karno,
    diagtime,
    age,
    prior
  )
names(analysis)
dim(analysis)
table(analysis$status)
table(analysis$time)
t.test(analysis$age)
t.test(analysis$time)
t.test(analysis$diagtime)
t.test(analysis$karno)
t.test(analysis$prior)
t.test(analysis$rownames)
t.test(age~trt,data=cancercell)
chisq.test(analysis$celltype,analysis$trt)
ggplot(analysis,aes(x=age))+
  geom_histogram()
ggplot(analysis,aes(x=status))+
  geom_bar()
ggplot(analysis,aes(x=time))+
  geom_bar()
ggplot(analysis,aes(x=diagtime))+
  geom_histogram()
chisq.test(table(cancercell$celltype,cancercell$trt))
km_model<-survfit(
  Surv(time,status)~1,
  data=analysis)
summary(km_model)
# Extract survival results

km_data <- data.frame(
  time = km_results$time,
  survival = km_results$surv,
  lower = km_results$lower,
  upper = km_results$upper
)
km_plot <- ggplot(
  km_data,
  aes(x = time, y = survival)
) +
  
  geom_step(
    linewidth = 1,
    colour = "blue"
  ) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2
  ) +
  
  labs(
    title = "Kaplan-Meier Survival Curve",
    x = "Survival Time (Days)",
    y = "Survival Probability"
  ) +
  
  theme_minimal()

print(km_plot)
logrank<-survdiff(
  Surv(time,event)~trt,
  data=analysis)
logrank
p_value<-1-pchisq(
  logrank$chisq,
  df=length(logrank$n)-1,
  lower.tail=FALSE)

p_value
cox_trt<-coxph(
  Surv(time,status)~factor(trt),
  data=analysis)
summary(cox_trt)
cox_age<-coxph(
  Surv(time,status)~age,
  data=analysis)

summary(cox_age)
cox_celltype<-coxph(
  Surv(time,status)~factor(celltype),
  data=analysis)
summary(cox_celltype)
cox_karno<-coxph(
  Surv(time,status)~karno,
  data=analysis)
summary(cox_karno)
cox_prior<-coxph(
  Surv(time,status)~factor(prior),
  data=analysis)
summary(cox_prior)
tidy(
  cox_trt,
  exponentiate=TRUE,
  conf.int=TRUE)
cox_results <- bind_rows(
  
  tidy(cox_trt, exponentiate = TRUE, conf.int = TRUE),
  tidy(cox_age, exponentiate = TRUE, conf.int = TRUE),
  tidy(cox_celltype, exponentiate = TRUE, conf.int = TRUE),
  tidy(cox_karno, exponentiate = TRUE, conf.int = TRUE),
  tidy(cox_prior, exponentiate = TRUE, conf.int = TRUE)
  
)

cox_results
cox_multi <- coxph(
  
  Surv(time, status) ~
    
    factor(trt) +
    factor(celltype) +
    age +
    karno +
    diagtime +
    factor(prior),
  
  data = analysis
)
summary(cox_multi)
exp(coef(cox_multi))
summary(cox_multi)$coefficients
model_summary <- summary(cox_multi)

cox_multi_results <- data.frame(
  
  Variable = rownames(model_summary$coefficients),
  
  Hazard_Ratio =
    model_summary$coefficients[, "exp(coef)"],
  
  Lower_CI =
    exp(confint(cox_multi))[, 1],
  
  Upper_CI =
    exp(confint(cox_multi))[, 2],
  
  P_value =
    model_summary$coefficients[, "Pr(>|z|)"]
  
)
cox_multi_results
summary(cox_multi)
HR <- exp(coef(cox_multi))
CI <- exp(confint(cox_multi))
p_values <- summary(cox_multi)$coefficients[, "Pr(>|z|)"]
forest_data <- data.frame(
  
  Variable = names(HR),
  
  HR = as.numeric(HR),
  
  Lower_CI = CI[, 1],
  
  Upper_CI = CI[, 2],
  
  P_value = p_values
  
)
forest_data

forest_plot <- ggplot(
  forest_data,
  aes(
    x = HR,
    y = reorder(Variable, HR)
  )
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_errorbar(
    aes(
      xmin = Lower_CI,
      xmax = Upper_CI
    ),
    width = 0.2
  ) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  
  scale_x_log10() +
  
  labs(
    title = "Forest Plot of Adjusted Hazard Ratios",
    x = "Adjusted Hazard Ratio (Log Scale)",
    y = "Clinical Predictors"
  ) +
  
  theme_minimal()

print(forest_plot)
ggsave(
  "outputs/figures/forest_plot.png",
  plot = forest_plot,
  width = 10,
  height = 7,
  dpi = 300
)
model_summary <- summary(cox_multi)
model_summary$concordance
concordance_result <- model_summary$concordance
concordance_result
concordance_result <- model_summary$concordance
concordance_result
summary(cox_multi)
AIC(cox_multi)
model_evaluation <- data.frame(
  
  C_index = model_summary$concordance[1],
  
  C_index_SE = model_summary$concordance[2],
  
  Likelihood_ratio =
    model_summary$logtest[1],
  
  Likelihood_ratio_p =
    model_summary$logtest[3],
  
  Wald_statistic =
    model_summary$waldtest[1],
  
  Wald_p =
    model_summary$waldtest[3],
  
  Score_statistic =
    model_summary$sctest[1],
  
  Score_p =
    model_summary$sctest[3],
  
  AIC = AIC(cox_multi)
  
)

model_evaluation
write.csv(
  model_evaluation,
  "outputs/tables/model_evaluation.csv",
  row.names = FALSE
)



