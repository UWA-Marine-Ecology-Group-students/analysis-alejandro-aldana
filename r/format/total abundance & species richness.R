## Example for getting total abundance & species richness data from your
## Count data

## edit as needed

ta.sr <-  complete.count %>%
  select(-c(family, genus, species))%>%
  pivot_wider(names_from = "scientific", values_from = count, values_fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_abundance = rowSums(.[, 6:(ncol(.))], na.rm = TRUE),
                species_richness = rowSums(.[, 6:(ncol(.))] > 0)) %>%
  dplyr::select(campaignid, sample, total_abundance, species_richness) %>%
  pivot_longer(cols = c("total_abundance", "species_richness"),
               names_to = "response", values_to = "number") %>%
  glimpse()

saveRDS(ta.sr, file = here::here(paste0("data/tidy/",
                                        name, "_ta.sr.rds")))