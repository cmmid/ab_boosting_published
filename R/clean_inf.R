read.csv(here::here("data", "hcw", "2022_2023_Flu_Swabs.csv"), header = TRUE, row.names = NULL) -> inf_data

inf_data_trim <- inf_data %>% filter(year == 2022) %>% select(pid, year, samp_date, Combined_Subtype_result) 

load(file = here::here("outputs", "data_model", "h3only_hcwonly_df.RData"))

inf_data_trim %>% pull(pid) %>% length

inf_data_trim %>% left_join(h3only_df %>% separate(pid, c("pid", "year"), sep = "_") %>% mutate(year = as.numeric(year)) ) %>% 
    filter(!is.na(base_titre)) %>% pull(pid) %>% unique %>% length
# 0 individuals in the 2020 and 2021 flu season
# 34 indivudals in the 2022 flu season
