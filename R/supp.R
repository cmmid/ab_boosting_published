
## ## ## ## ## ## ## ## ## ##
## ## Correlation plots ## ##
## ## ## ## ## ## ## ## ## ##


recode_vhist <- c("1" = "<2 vaccines in last 5 seasons", "2" = "<2 vaccines in last 5 seasons", 
    "3" = "2 or more vaccines in last 5 seasons", "4" = "2 or more vaccines in last 5 seasons",
    "5" = "2 or more vaccines in last 5 seasons", "6" = "2 or more vaccines in last 5 seasons")

get_ind_boosts <- function(best_fit_hXonly, hXonly_stan) {
    boostingwane_t <- best_fit_hXonly %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])

    meta_join_h1 <- data.frame(
        i = 1:hXonly_stan$N_ind,
        titre_i = hXonly_stan$titre_i,
        vac_hist = hXonly_stan$vac_hist
    )
    df_boost_uncert_meta_h1 <- boostingwane_t %>% left_join(meta_join_h1) %>% ungroup %>%
        summarise(boost_ind = mean(boost_ind), wane_ind = mean(wane_ind), .by = c(titre_i, i, vac_hist)) %>%
        mutate(vac_hist = recode(vac_hist, !!!recode_vhist)) %>% rename(t = titre_i) %>% add_titre_info_trunc
}

ind_traj_h1h1egg <- get_ind_boosts(best_fit_h1only, h1only_stan) %>% mutate(subtype = "A(H1N1) vaccinating")
ind_traj_h3h2egg <- get_ind_boosts(best_fit_h3only, h3only_stan) %>% mutate(subtype = "A(H3N2) vaccinating")
ind_traj_h1h1cell <- get_ind_boosts(best_fit_h1cell, h1cell_stan) %>% mutate(subtype = "A(H1N1) cell-grown")
ind_traj_h3h2cell <- get_ind_boosts(best_fit_h3cell, h3cell_stan) %>% mutate(subtype = "A(H3N2) cell-grown")



p1 <- ind_traj_h1h1egg %>% bind_rows(ind_traj_h3h2egg) %>% select(!c(titre_vals, t)) %>% 
    pivot_wider(names_from = "subtype", values_from = c("boost_ind", "wane_ind")) %>%
    ggplot() + 
        geom_abline(intercept = 0, slope = 1, color = "gray") +
        geom_point(aes(x = `boost_ind_A(H1N1) vaccinating`, y  = `boost_ind_A(H3N2) vaccinating`), alpha = 0.5) + 
        geom_smooth(method = "lm", aes(x = `boost_ind_A(H1N1) vaccinating`, y  = `boost_ind_A(H3N2) vaccinating`)) + 
        scale_y_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + theme_bw() + 
        scale_x_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + 
        labs(x = "Fold-rise, A(H1N1) vaccinating", y = "Fold-rise A(H3N2) vaccinating")

p2 <- ind_traj_h1h1cell %>% bind_rows(ind_traj_h3h2cell) %>% select(!c(titre_vals, t)) %>% 
    pivot_wider(names_from = "subtype", values_from = c("boost_ind", "wane_ind")) %>%
    ggplot() + 
        geom_abline(intercept = 0, slope = 1, color = "gray") +
        geom_point(aes(x = `boost_ind_A(H1N1) cell-grown`, y  = `boost_ind_A(H3N2) cell-grown`), alpha = 0.5) + 
        geom_smooth(method = "lm", aes(x = `boost_ind_A(H1N1) cell-grown`, y  = `boost_ind_A(H3N2) cell-grown`)) + 
        scale_y_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + theme_bw() + 
        scale_x_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + 
        labs(x = "Fold-rise, A(H1N1) cell-grown", y = "Fold-rise A(H3N2) cell-grown")

        
p1 + p2
ggsave(here::here("outputs", "figs", "supp", "comparison.pdf"))      


         +

df_boost_uncert_h1 <- best_fit_h1only %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])
meta_join_h1 <- data.frame(
    i = 1:h1only_stan$N_ind,
    titre_i = h1only_stan$titre_i,
    vac_hist = h1only_stan$vac_hist
)
df_boost_uncert_meta_h1 <- df_boost_uncert_h1 %>% left_join(meta_join_h1) %>% ungroup


df_boost_uncert_h3 <- best_fit_h3only %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])
meta_join_h3 <- data.frame(
    i = 1:h3only_stan$N_ind,
    titre_i = h3only_stan$titre_i,
    vac_hist = h3only_stan$vac_hist
)
df_boost_uncert_meta_h3 <- df_boost_uncert_h3 %>% left_join(meta_join_h3) %>% ungroup


df_boost_uncert_h3_cell <- best_fit_h3cell %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])
meta_join_h3_cell <- data.frame(
    i = 1:h3cell_stan$N_ind,
    titre_i = h3cell_stan$titre_i,
    vac_hist = h3cell_stan$vac_hist
)
df_boost_uncert_meta_h3_cell <- df_boost_uncert_h3_cell %>% left_join(meta_join_h3_cell) %>% ungroup

df_boost_uncert_meta_un <- bind_rows(
    df_boost_uncert_meta_h1 %>% mutate(subtype = "A(H1N1) vaccinating"),
    df_boost_uncert_meta_h3 %>% mutate(subtype = "A(H3N2) vaccinating"),
    df_boost_uncert_meta_h3_cell %>% mutate(subtype = "A(H3N2) circulating"),
)


individual_level_bw <- df_boost_uncert_meta_un %>%
    summarise(boost_ind = mean(boost_ind), wane_ind = mean(wane_ind), .by = c(titre_i, i, subtype, vac_hist)) %>%
    mutate(titre_vals = recode(titre_i, !!!titres_labels)) %>% mutate(titre_vals = factor(titre_vals, levels = titres_labels)) %>%
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating",  "A(H3N2) vaccinating", "A(H3N2) circulating")))

plot_compare <- individual_level_bw %>% select(!c(wane_ind, titre_vals, vac_hist, titre_i)) %>% 
    pivot_wider(names_from = "subtype", values_from = "boost_ind") 

fit_boost <- lm(`A(H3N2) vaccinating` ~ `A(H1N1) vaccinating`, data = plot_compare)
fit_wane <- lm(`A(H3N2) vaccinating` ~ `A(H3N2) circulating`, data = plot_compare)
summary(fit_boost) # 0.1299
summary(fit_wane) # 0.1323

p1 <- plot_compare %>% 
    ggplot() + 
        geom_abline(intercept = 0, slope = 1, color = "gray") +
        geom_point(aes(`A(H1N1) vaccinating`,  `A(H3N2) vaccinating`), alpha = 0.6) + 
        geom_smooth(aes(`A(H1N1) vaccinating`, `A(H3N2) vaccinating`), method = "lm", formula= y~x) + 
        scale_y_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + 
        scale_x_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + 
        theme_bw()

p2 <- plot_compare %>% 
    ggplot() + 
        geom_abline(intercept = 0, slope = 1, color = "gray") +
        geom_point(aes(`A(H3N2) vaccinating`, `A(H3N2) circulating`), alpha = 0.6) + 
        geom_smooth(aes(`A(H3N2) vaccinating`, `A(H3N2) circulating`), method = "lm", formula= y~x) + 
        scale_y_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + 
        scale_x_continuous(breaks = -1:8, labels = 2^c(-1:8) ) + 
        theme_bw()
p1 + p2
ggsave(here::here("outputs", "figs", "supp", "comparison.pdf"))