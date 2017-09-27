library(dplyr)
library(readr)
library(cowplot)

my_data <- read_csv("../map_structural_features/3rze.map.rates_features.csv") %>%
  mutate(r4s_rate_norm = r4s_rate/mean(r4s_rate))  # Normalize rates to a mean of 1

codon_wcn <- ggplot(my_data, aes(x = wcn_sc, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('side-chain WCN')
codon_rsa <- ggplot(my_data, aes(x = rsa, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('RSA')

aa_wcn <- ggplot(my_data, aes(x = wcn_sc, y = `r4s_rate_norm`)) + 
  geom_point(alpha = 0.5) +
  labs(x = 'side-chain WCN', y = 'rate4site score')
aa_rsa <- ggplot(my_data, aes(x = rsa, y = `r4s_rate_norm`)) + 
  geom_point(alpha = 0.5) +
  labs(x = 'RSA', y = 'rate4site score')

all_plots <- plot_grid(codon_wcn, codon_rsa, aa_wcn, aa_rsa, nrow = 2, ncol = 2, align='vh', labels = 'auto', scale = 0.95)
save_plot("../figures/scatterplots.pdf", all_plots, base_width = 8, base_height = 7)
