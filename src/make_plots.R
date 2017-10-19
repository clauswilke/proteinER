library(dplyr)
library(readr)
library(cowplot)

my_data <- read_csv("../map_structural_features/3rze.map.rates_features.csv") %>%
  mutate(Rate_norm = Rate/mean(Rate))  # Normalize rates to a mean of 1

codon_wcn <- ggplot(my_data, aes(x = wcn_sc, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('side-chain WCN')
codon_rsa <- ggplot(my_data, aes(x = rsa, y = `dN/dS`)) + 
  geom_point(alpha = 0.5) +
  xlab('RSA')

aa_wcn <- ggplot(my_data, aes(x = wcn_sc, y = `Rate_norm`)) + 
  geom_point(alpha = 0.5) +
  labs(x = 'side-chain WCN', y = 'Relative rate')
aa_rsa <- ggplot(my_data, aes(x = rsa, y = `Rate_norm`)) + 
  geom_point(alpha = 0.5) +
  labs(x = 'RSA', y = 'Relative rate')

all_plots <- plot_grid(codon_wcn, codon_rsa, aa_wcn, aa_rsa, nrow = 2, ncol = 2, align='vh', labels = 'auto', scale = 0.95)
save_plot("../figures/scatterplots.pdf", all_plots, base_width = 8, base_height = 8)

#correlation tests 
#WCN and dN/dS
cor.test(my_data$wcn_sc,my_data$`dN/dS`)
#RSA and dN/dS
cor.test(my_data$rsa,my_data$`dN/dS`)
#WCN and relative amino acid rate
cor.test(my_data$wcn_sc,my_data$Rate_norm)
#RSA and relative amino acid rate
cor.test(my_data$rsa,my_data$Rate_norm)