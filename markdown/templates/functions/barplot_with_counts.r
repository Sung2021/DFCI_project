obj.srt@meta.data %>%
  group_by(!!sym(res), orig.ident) %>%
  summarise(count = n()) %>%
  group_by(!!sym(res)) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ggplot(aes(x = !!sym(res), y = percent, fill = orig.ident)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%", percent)), 
            position = position_stack(vjust = 0.5),
            show.legend = FALSE, size=3) +
  labs(y = "Percentage") +
  theme_minimal()
