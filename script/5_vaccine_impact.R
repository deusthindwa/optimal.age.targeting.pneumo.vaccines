# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 1/08/2021-30/09/2021

# define efficacy waning using studies
VE_impact_by_age <- VE_by_Vac.age %>%
    dplyr::group_by(serogroup,
                    age_dep,
                    sim,
                    Vac.age,
                    country) %>%
    
    dplyr::summarise(Impact = sum(Impact)) 


# add uncertainty to VE impact
VE_impact_by_age_ <- VE_impact_by_age %>%
    filter(!is.na(Impact)) %>%
    nest(data = c(sim, Impact)) %>%
    mutate(Q = map(data, ~quantile(x     = .x$Impact,
                                   probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q) %>%
    select(-data)


make_grid_plot <- function(x, ylab = NULL, percent = FALSE, ylim = c(0,NA)){
    
    
    p <- ggplot(x,
                aes(x = Vac.age, y = `50%`,
                    color = factor(age_dep),
                    group = interaction(age_dep, serogroup, country))) +
        geom_line() + 
        geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, 
                        fill = factor(age_dep)), color = NA, alpha = 0.2) +
        geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, 
                        fill = factor(age_dep)), color = NA, alpha = 0.2) +
        facet_nested(country ~ serogroup, scales = "free_y") +
        theme_bw() +
        scale_x_continuous(breaks = seq(60, 90, 10)) +
        theme(axis.text=element_text(size=10, color="black")) +
        xlab("Vaccination Age (years)") +
        ylab(ylab) +
        theme_bw(base_size = 14, base_family = "Lato") +
        theme(axis.text        = element_text(face = "bold"),
              strip.background = element_rect(fill = "white"),
              panel.border     = element_rect(colour = "black", fill=NA, size=1)) +
        theme(legend.position = "bottom") +
        scale_color_brewer(name = "Age dependent vaccine efficacy", palette = "Set1") + 
        scale_fill_brewer(name = "Age dependent vaccine efficacy", palette = "Set1")
    
    if (percent){
        p <- p + scale_y_continuous(labels = function(x){sprintf("%g%%",x*100)},
                                    limits = ylim)
    } else {
        p <- p + scale_y_continuous(limits = ylim)
    }
    
    p
}

# plot vaccine impact as in expected number of cases averted

VE_A <- make_grid_plot(x = VE_impact_by_age_, 
                       ylab = "Vaccine impact (expected total cases averted)")

ggsave(filename = "output/Fig2_vaccine_impact_cohort.png", 
       plot = VE_A,
       width = 14, height = 8, units = "in", dpi = 300)

#===============================================================================================

#vaccine impact per 100,000 older adults vaccinated from specific age cohort (ntotal)
VE_impact_validated <- dplyr::select(pop_country_df, country, agey, ntotal) %>% 
    dplyr::rename(Vac.age = agey) %>% 
    dplyr::inner_join(VE_impact_by_age, by = c("country", "Vac.age")) %>%
    mutate(Impact = Impact/ntotal*scale) %>%
    nest(data = c(sim, Impact)) %>%
    mutate(Q = map(data, ~quantile(.x$Impact, probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q)

#impact per 100,000 older adults vaccinated from specific age cohort (ntotal)
VE_B <- make_grid_plot(x = VE_impact_validated, 
                       ylab = "Impact (cases averted per 100,000 vaccinees)")  

ggsave(filename = "output/Fig3_vaccine_impact_vaccinee.png", 
       plot = VE_B,
       width = 14, height = 8, units = "in", dpi = 300)



#===============================================================================================
# Sam's
#===============================================================================================

impact_per_case <- ipd_mc %>%
    mutate(cases = map(.x = mc, .f = ~group_by(.x, sim) %>%
                           # make it per vaccinee
                           crossing(Vac.age = seq(55, 85, by = 5)) %>%
                           filter(agey >= Vac.age) %>%
                           # end per vaccinee
                           group_by(Vac.age, sim) %>%
                           summarise(cases = sum(fit), .groups = 'drop'))) %>%
    select(-data, -model, -mc) %>%
    unnest(cases) %>%
    inner_join(VE_impact_by_age) %>%
    mutate(rel_impact = Impact/cases) %>%
    group_by_at(.vars = vars(-c(sim, cases, Impact, rel_impact))) %>%
    nest %>%
    mutate(Q = map(.x = data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q)

plot_impact_per_case <- 
  make_grid_plot(x    = impact_per_case,  percent = TRUE,
                 ylab = "Vaccine impact (proportion of cases averted among vaccinees)")

ggsave(filename = "output/vaccine_impact_per_case.png", 
       plot = plot_impact_per_case,
       width = 14, height = 8, units = "in", dpi = 300)

#===============================================================================================

#vaccine impact per 100,000 total population 55+y (pop)
impact_per_vaccinee <- 
    
    # get the population over 55, 60, etc. as potential vaccinees
    pop_country_df %>%
    crossing(Vac.age = seq(55, 85, by = 5)) %>%
    filter(agey >= Vac.age) %>%
    group_by(country, Vac.age) %>%
    summarise(pop = sum(ntotal)) %>%
    
    # merge with Impact data (averted cases, absolute)
    inner_join(VE_impact_by_age) %>%
    
    # relative impact is per 100,000 total population 55+y
    mutate(rel_impact = Impact/pop*scale) %>%
    group_by_at(.vars = vars(-c(sim, pop, Impact, rel_impact))) %>%
    nest %>%
    mutate(Q = map(.x = data, ~quantile(.x$rel_impact, probs = c(0.025, 0.5, 0.975)))) %>%
    unnest_wider(Q)

plot_impact_per_vaccinee <- 
    make_grid_plot(x    = impact_per_vaccinee, 
                   ylab = "Vaccine impact (Cases averted per 100,000 vaccinees)")

ggsave(filename = "output/vaccine_impact_per_vaccinee.png", 
       plot = plot_impact_per_vaccinee,
       width = 14, height = 8, units = "in", dpi = 300)
