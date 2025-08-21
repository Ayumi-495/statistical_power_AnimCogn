# power_firepower_all, firepower.M_all, firepower.S_all

# prep ----
prep <- function(df, value_col) {
  df|> 
    mutate(
      level = case_when(
        str_detect(effect, "Sampling") ~ "Primary study",
        str_detect(effect, "Meta-analysis") ~ "Meta-analysis level",
        TRUE ~ "Other"
      ),
      correction = if_else(str_starts(effect, "c"), "Bias-corrected", "Uncorrected")
    )
}

scatter_power <- prep(power_firepower_all,  power)    
scatter_M <- prep(firepower.M_all,  power)      
scatter_S <- prep(firepower.S_all,  power)
names(power_firepower_all)

scatter_power <- scatter_power %>%
  mutate(level = fct_relevel(level, "Primary study", "Meta-analysis level"))
scatter_M <- scatter_M %>%
  mutate(level = fct_relevel(level, "Primary study", "Meta-analysis level"))
scatter_S <- scatter_S %>%
  mutate(level = fct_relevel(level, "Primary study", "Meta-analysis level"))

colour <- c("Bias-corrected" = "#52595D",  
         "Uncorrected"   = "#B6B6B6")

## power ----
# head(scatter_power[order(-scatter_power$power), ], 10)
# head(scatter_power[order(scatter_power$power), ], 10)

scatter_power$correction <- factor(
  scatter_power$correction,
  levels = c("Uncorrected", "Bias-corrected")  
  )

p_scatter_power <- ggplot(
  scatter_power,
  aes(x = correction, y = power, fill = correction)
  ) +
    geom_violin(alpha = 0.6, width = 0.9, trim = FALSE, color = "#DC143C",
              position = position_dodge(width = 0.8)) +
  geom_point(
    aes(fill = correction),
    shape = 21,
    color = "#DC143C",
    size = 2.5, alpha = 0.6,
    position = position_dodge(width = 0.8,
                              # jitter.width = 0.12
                              )
  ) +
  # geom_dotplot(
  #   aes(fill = correction),
  #   binaxis = "y", stackdir = "center",
  #   dotsize = 0.5, color = "deeppink3",
  #   position = position_dodge(width = 0.8)
  # ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    aes(fill = correction),
    color = "#DC143C",
    alpha = 0.8,
    width = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ level, 
             nrow = 1, 
             scales = "free_y") +
  labs(x = NULL, y = "Proportion",
       fill = "Correction", 
       color = "Correction",
       title = "Power") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = colour,
                    breaks = rev(names(colour)),  
                    limits = rev(names(colour))) +
  scale_colour_manual(values = colour,
                      breaks = rev(names(colour)), 
                      limits = rev(names(colour))) + 
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA)
  )

p_scatter_power

## type M error ----
# summary(scatter_M$power)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.997     1.156     2.110   107.366     4.991 14841.502 

# head(scatter_M[order(-scatter_M$power), ], 10) 
# MA_case            es     power         effect         es_cat               level     correction
# 51   SMD 31  6.492139e-05 14841.502      cSampling a.Intercept_es       Primary study Bias-corrected
# 115  SMD 31  6.492139e-05  4233.801 cMeta-analysis     b.MA.power Meta-analysis level Bias-corrected
# 43   SMD 23 -1.757052e-03   616.574      cSampling a.Intercept_es       Primary study Bias-corrected
# 107  SMD 23 -1.757052e-03    77.506 cMeta-analysis     b.MA.power Meta-analysis level Bias-corrected
# 158   Zr 10  8.735408e-03    67.069      cSampling a.Intercept_es       Primary study Bias-corrected
# 22    SMD 2 -2.913212e-02    42.823      cSampling a.Intercept_es       Primary study Bias-corrected
# 34   SMD 14 -2.989572e-02    41.926      cSampling a.Intercept_es       Primary study Bias-corrected
# 98   SMD 14 -2.989572e-02    24.487 cMeta-analysis     b.MA.power Meta-analysis level Bias-corrected
# 84   SMD 32  5.200886e-02    22.752       Sampling a.Intercept_es       Primary study    Uncorrected
# 52   SMD 32  5.200886e-02    22.531      cSampling a.Intercept_es       Primary study Bias-corrected

# head(scatter_M[order(scatter_M$power), ], 10)
# MA_case         es power         effect     es_cat               level     correction
# 113  SMD 29 -0.2360722 0.997 cMeta-analysis b.MA.power Meta-analysis level Bias-corrected
# 118   SMD 2 -1.3847886 0.998  Meta-analysis b.MA.power Meta-analysis level    Uncorrected
# 101  SMD 17 -0.9161862 0.999 cMeta-analysis b.MA.power Meta-analysis level Bias-corrected
# 117   SMD 1 -1.4714368 0.999  Meta-analysis b.MA.power Meta-analysis level    Uncorrected
# 125   SMD 9  0.3359196 0.999  Meta-analysis b.MA.power Meta-analysis level    Uncorrected
# 133  SMD 17 -0.9161862 0.999  Meta-analysis b.MA.power Meta-analysis level    Uncorrected
# 137  SMD 21  0.4057026 0.999  Meta-analysis b.MA.power Meta-analysis level    Uncorrected
# 138  SMD 22  0.6633060 0.999  Meta-analysis b.MA.power Meta-analysis level    Uncorrected
# 146  SMD 30  1.1628672 0.999  Meta-analysis b.MA.power Meta-analysis level    Uncorrected
# 178    Zr 8  0.4936675 0.999 cMeta-analysis b.MA.power Meta-analysis level Bias-corrected

# I need to transform the power to a log scale as five powers
scatter_M$correction <- factor(
  scatter_M$correction,
  levels = c("Uncorrected", "Bias-corrected")  
)

p_scatter_M <- ggplot(
  scatter_M,
  aes(x = correction, y = log(power), fill = correction)
) +
  geom_violin(alpha = 0.6, 
              width = 0.9, 
              trim = FALSE, 
              color = "#1B8A6B",
              position = position_dodge(width = 0.8)
              ) +
  geom_point(
    aes(fill = correction),
    shape = 21,
    color = "#1B8A6B",
    size = 2.5, alpha = 0.6,
    position = position_dodge(width = 0.8,
                              # jitter.width = 0.12
    ) 
  )+
  stat_summary(
    fun = mean,
    geom = "crossbar",
    aes(fill = correction),
    color = "#1B8A6B",
    alpha = 0.8,
    width = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ level,
             nrow = 1, 
             scales = "free_y") +
  labs(x = NULL, y = "Log-transformed magnitude",
       fill = "Correction", 
       color = "Correction",
       title = "Type M error") +
  scale_y_continuous(limits = c(-1, 10)) +
  # coord_cartesian(ylim = c(min(log(scatter_M$power), na.rm = TRUE), 10)) +
  scale_fill_manual(values = colour,
                    breaks = rev(names(colour)),  
                    limits = rev(names(colour))) +
  scale_colour_manual(values = colour,
                      breaks = rev(names(colour)), 
                      limits = rev(names(colour))) + 
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA)
  )

p_scatter_M

## type S error ----
head(scatter_S[order(-scatter_S$power), ], 10)
head(scatter_S[order(scatter_S$power), ], 10)

scatter_S$correction <- factor(
  scatter_S$correction,
  levels = c("Uncorrected", "Bias-corrected")  
)

p_scatter_S <- ggplot(
  scatter_S,
  aes(x = correction, y = power, fill = correction)
) +
  geom_violin(alpha = 0.6, width = 0.9, trim = FALSE, color = NA,
              position = position_dodge(width = 0.8)) +
  
  geom_violin(alpha = 0.6, 
              width = 0.9, 
              trim = FALSE, 
              color = "#F5E216",
              position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(fill = correction),
    shape = 21,
    color = "#F5E216",
    size = 2.5, alpha = 0.6,
    position = position_dodge(width = 0.8,
                              # jitter.width = 0.12
    ) 
  )+
  stat_summary(
    fun = mean,
    geom = "crossbar",
    aes(fill = correction),
    color = "#F5E216",
    alpha = 0.8,
    width = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ level, nrow = 1, scales = "free_y") +
  labs(x = NULL, 
       y = "Proportion",
       fill = "Correction", 
       color = "Correction",
       title = "Type S error") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = colour,
                    breaks = rev(names(colour)),  
                    limits = rev(names(colour))) +
  scale_colour_manual(values = colour,
                      breaks = rev(names(colour)), 
                      limits = rev(names(colour))) + 
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = NA)
  )


p_scatter_S

# all plots ----
(p_scatter_power + p_scatter_M + p_scatter_S) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "none")
