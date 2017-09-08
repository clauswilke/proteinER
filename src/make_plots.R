library(dplyr)
library(readr)
library(cowplot)

my_data <- read_csv("../aligning_structural_features/3rze.map.rates_features.csv")

codon_wcn <- ggplot(my_data, aes(x = wcn_sc, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('side-chain WCN')
codon_rsa <- ggplot(my_data, aes(x = rsa, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('RSA')

aa_wcn <- ggplot(my_data, aes(x = wcn_sc, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('side-chain WCN')
aa_rsa <- ggplot(my_data, aes(x = rsa, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('RSA')

all_plots <- plot_grid(codon_wcn, codon_rsa, aa_wcn, aa_rsa, nrow = 2, ncol = 2, labels = 'auto')
save_plot("../figures/scatterplots.pdf", all_plots, base_width = 8, base_height = 7)
