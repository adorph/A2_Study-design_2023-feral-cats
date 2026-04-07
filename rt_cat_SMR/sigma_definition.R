set.seed(123)

n <- 50000

# Draw from distribution
beta0 <- rnorm(n, log(304), 0.5)
# Calculated SD at 0.45, increased to a conservative 0.7 to avoid overconfidence because of the small 
# female sample size from the area (n=2) and transferability of estiamtes from other studies
beta_sex <- rnorm(n, log(584) - log(304), 0.75)
# Draw from distribution
# beta0 <- rnorm(n, 0, 100)
# beta_sex <- rnorm(n, 0, 100)

# Female and male
sigma_female <- exp(beta0)
sigma_male <- exp(beta0 + beta_sex)

df_sim <- data.frame(
  sigma = c(sigma_female, sigma_male),
  sex = rep(c("Female", "Male"), each = n)
)

ggplot(df_sim, aes(x = sigma, fill = sex)) +
  annotate("rect",
           xmin = 0, xmax = 50,
           ymin = -Inf, ymax = Inf,
           fill = "grey") +
  annotate("rect",
           xmin = 2000, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "grey") +
  geom_density(alpha = 0.4) +
  geom_vline(aes(xintercept=median(subset(df_sim, df_sim$sex == "Female")$sigma), lty="median"), color="#E69F00", size=1) +
  geom_vline(aes(xintercept=median(subset(df_sim, df_sim$sex == "Male")$sigma), lty="median"), color="#009E73", size=1) +
  geom_vline(aes(xintercept=mean(subset(df_sim, df_sim$sex == "Female")$sigma), lty="mean"), color="#E69F00", alpha=0.5, size=1) +
  geom_vline(aes(xintercept=mean(subset(df_sim, df_sim$sex == "Male")$sigma), lty="mean"), color="#009E73", alpha=0.5, size=1) +   
  scale_fill_manual(values=c("#E69F00", "#009E73")) +
  scale_x_continuous(expand=c(0,0),limits=c(0, 4000), breaks=c(0, 1000, 2000, 3000, 4000)) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="sigma (m)", y= "Density")+
  labs(fill = "Sex", lty = "Statistic")+
  theme_bw()
min(df_sim$sigma)
max(df_sim$sigma)

# We conducted prior predictive checks by simulating σ from the specified priors. 
# These priors produced median values consistent with telemetry estimates (304 m and 584 m), 
# while allowing 2–5 fold variation. Extreme values (>2 km or <50 m) occurred with low probability,
# indicating the priors were weakly informative but biologically plausible.



ggsave(file="outputs/sigma_distribution.png", width=18, height=12, units="cm", dpi=300)
