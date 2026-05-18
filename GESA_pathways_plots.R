library(tidyverse)

Downpath <- read.delim("Downregulated pathways.txt", header = T,sep = "\t")
Uppath <- read.delim("upregulated pathways.txt", header = T,sep = "\t")
#-------------------------------------------------------------------------------
#UPregulation pathway plot
#-------------------------------------------------------------------------------
up_data <- Uppath %>%
  arrange(desc(NES)) %>%
  slice_head(n = 12) %>%
  arrange(NES) %>%
  mutate(NAME = factor(NAME, levels = NAME))
ggplot(up_data, aes(x = NES, y = NAME)) +
  geom_point(aes(size = SIZE, color = FDR.q.val)) + 
  scale_color_gradient(low = "darkred", high = "mistyrose") + 
  theme_classic() +
  labs(
    size = "Gene Set Size",
    color = "FDR q-value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.text.x = element_text(size = 10, color = "black")
  )
#-------------------------------------------------------------------------------
#Downregulation pathway plot
#-------------------------------------------------------------------------------
down_data <- Downpath %>%
  arrange(NES) %>%
  slice_head(n = 12) %>%
  arrange(desc(NES)) %>%
  mutate(NAME = factor(NAME, levels = NAME))

ggplot(down_data, aes(x = NES, y = NAME)) +
  geom_point(aes(size = SIZE, color = FDR.q.val)) + 
  scale_color_gradient(low = "darkblue", high = "lightblue") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  labs(
    size = "Gene Set Size",
    color = "FDR q-value"
  ) +
  theme(
    legend.position = "left",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 10, color = "black")
  )
