create_df_comparison <- function(fit_hX, stan_hX) {
    mean_boost_id <- fit_hX %>% as_draws_df %>% spread_draws(titre_est[s]) %>% 
    group_by(s) %>%
    summarise(titre_est_mean = mean(titre_est)) %>% mutate(titre_est_mean_round = round(titre_est_mean, 0)) 

    data.frame(
            pid = stan_hX$ind,
            data = stan_hX$boost,
            model = mean_boost_id$titre_est_mean_round
        ) %>% 
        left_join(    
            data.frame(
                pid = 1:stan_hX$N_ind,
                study = stan_hX$study
            )
    )
}

create_df_comparison_pop <- function(fit_hX, stan_hX) {
    mean_boost_id_ind <- fit_hX %>% as_draws_df %>% spread_draws(titre_est[s]) %>% 
    group_by(s) %>%
    summarise(titre_est_mean = mean(titre_est))

    mean_boost_id_pop <- fit_hX %>% as_draws_df %>% spread_draws(z_boost_ind[j, i], sigma_b_ind[j], z_wane_ind[j, i], sigma_w_ind[j]) %>% 
        group_by(i) %>%
        summarise(z_boost_ind = mean(z_boost_ind), sigma_b_ind = mean(sigma_b_ind), z_wane_ind = mean(z_wane_ind), sigma_w_ind = mean(sigma_w_ind)) 

    ids <- stan_hX$ind
    days <- stan_hX$day
    model_pred <- vector()
    for (s in 1:nrow(mean_boost_id_ind)) {
        model_pred[s] <- mean_boost_id_ind[[s, 2]] - mean_boost_id_pop[[ids[s], 2]] * mean_boost_id_pop[[ids[s], 3]] - (mean_boost_id_pop[[ids[s], 4]] *  mean_boost_id_pop[[ids[s], 5]]) * days[s]
    }



    data.frame(
            pid = stan_hX$ind,
            data = stan_hX$boost,
            model = round(model_pred, 0)
        ) %>% 
        left_join(    
            data.frame(
                pid = 1:stan_hX$N_ind,
                study = stan_hX$study
            )
    )
}

load(file = here::here("outputs", "data_model", "h1only_hcwonly_stan.RData"))
load(file = here::here("outputs", "data_model", "h3only_hcwonly_stan.RData"))
load(file = here::here("outputs", "data_model", "h1cell_hcwonly_stan.RData"))
load(file = here::here("outputs", "data_model", "h3cell_hcwonly_stan.RData"))

best_fit_h3only <- readRDS(here::here("outputs", "stan", "fit_h3only_hcwonly_base.RData"))
best_fit_h1only <- readRDS(here::here("outputs", "stan", "fit_h1only_hcwonly_base.RData"))
best_fit_h1cell <- readRDS(here::here("outputs", "stan", "fit_h1cell_hcwonly_base.RData"))
best_fit_h3cell <- readRDS(here::here("outputs", "stan", "fit_h3cell_hcwonly_base.RData"))

plot_data_h3vac <- create_df_comparison(best_fit_h3only, h3only_stan)
plot_data_h1vac <- create_df_comparison(best_fit_h1only, h1only_stan)
plot_data_h3cell <- create_df_comparison(best_fit_h3cell, h3cell_stan)
plot_data_h1cell <- create_df_comparison(best_fit_h1cell, h1cell_stan)


p1 <- plot_data_h1vac %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        geom_smooth(aes(data, model), method = "lm", size = 1.5, color = "white") +
        geom_smooth(aes(data, model), method = "lm", size = 0.5, color = "darkred") +  
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits for A(H1N1) vaccinating strains") + guides(size = "none")

p2 <- plot_data_h3vac %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        geom_smooth(aes(data, model), method = "lm", size = 1.5, color = "white") +
        geom_smooth(aes(data, model), method = "lm", size = 0.5, color = "darkred") +  
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits for A(H3N2) vaccinating strains")  + guides(size = "none")

p3 <- plot_data_h1cell %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        geom_smooth(aes(data, model), method = "lm", size = 1.5, color = "white") +
        geom_smooth(aes(data, model), method = "lm", size = 0.5, color = "darkred") +  
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits for A(H1N1) cell-grown strains") + guides(size = "none")

p4 <- plot_data_h3cell %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        geom_smooth(aes(data, model), method = "lm", size = 1.5, color = "white") +
        geom_smooth(aes(data, model), method = "lm", size = 0.5, color = "darkred") +  
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits for A(H3N2) cell-grown strains") + guides(size = "none")

p1 / p2 / p3 / p4 + plot_annotation(tag_levels = "A")
ggsave(file = here::here("outputs", "figs", "main", "fig0.pdf"), height = 12, width = 10)


######## ######## ######## ######## ######## 
######## Comparison when individual effects removed ######## 
######## ######## ######## ######## ######## 


plot_data_h3vac_pop <- create_df_comparison_pop(best_fit_h3only, h3only_stan)
plot_data_h1vac_pop <- create_df_comparison_pop(best_fit_h1only, h1only_stan)
plot_data_h3cell_pop <- create_df_comparison_pop(best_fit_h3cell, h3cell_stan)
plot_data_h1cell_pop <- create_df_comparison_pop(best_fit_h1cell, h1cell_stan)

p1 <- plot_data_h1vac_pop %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits (no ind-effects) for A(H1N1) vaccinating strains") + guides(size = "none")

p2 <- plot_data_h3vac_pop %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits (no ind-effects) for A(H3N2) vaccinating strains")  + guides(size = "none")

p3 <- plot_data_h1cell_pop %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits (no ind-effects) for A(H1N1) cell-grown strains") + guides(size = "none")

p4 <- plot_data_h3cell_pop %>%  mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + geom_count(aes(data, model), alpha = 0.5, shape = 21, fill = "red") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        theme_bw() + facet_grid(cols = vars(study)) + 
        scale_y_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        scale_x_continuous(breaks = seq(-2, 10, 2), labels = 2^seq(-2, 10, 2)) + 
        labs(x = "HAI fold-rise from data (GTR) at bleed date",
            y = "Model-predicted HAI fold-rise (GTR)\n at bleed date") + 
        ggtitle("Model fits (no ind-effects) for A(H3N2) cell-grown strains") + guides(size = "none")

p1 / p2 / p3 / p4 + plot_annotation(tag_levels = "A")
ggsave(file = here::here("outputs", "figs", "main", "fig0_pop.pdf"), height = 12, width = 10)


######## ######## ######## ######## ######## 
######## Metrics for manuscript ######## 
######## ######## ######## ######## ######## 

# Determine error rate (Within one-fold change unit)

plot_data_h1vac %>% mutate(diff = abs(data - model)) %>% 
    mutate(in_bounds = (diff <= 1)) %>% summarise(n = n() / nrow(.), .by = in_bounds)
plot_data_h3vac %>% mutate(diff = abs(data - model)) %>% 
    mutate(in_bounds = (diff <= 1)) %>% summarise(n = n() / nrow(.), .by = in_bounds)
plot_data_h3cell %>% mutate(diff = abs(data - model)) %>% 
    mutate(in_bounds = (diff <= 1)) %>% summarise(n = n() / nrow(.), .by = in_bounds)


plot_data_h1vac %>%  mutate(above_16 = (data >= 5)) %>% summarise(n = n() / nrow(.), .by = above_16)
plot_data_h3vac %>%  mutate(above_16 = (data >= 5)) %>% summarise(n = n() / nrow(.), .by = above_16)
plot_data_h3cell %>%  mutate(above_16 = (data >= 5)) %>% summarise(n = n() / nrow(.), .by = above_16)