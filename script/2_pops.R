# written by Samuel Clifford & Deus Thindwa
# optimal age targeting for pneumoccocal vaccines against IPD in older adults
# 31/02/2022

#loaad datasets into memory
pop_en <- read_csv("data/total_pop_en.csv")
pop_mw <- read_csv("data/total_pop_mw.csv")
pop_sa <- read_csv("data/total_pop_sa.csv")
pop_br <- read_csv("data/total_pop_br.csv")

#ggplot comparing % populations in England/Wales versus UN SDG regions
pop_country <- c("Brazil" = "blue",
                 "England" = "red", 
                 "Malawi" = "Black",
                 "South Africa" = "orange")

pop_totals <- list(`England` = 16492465, # mid-2017 in England
                   `Malawi`        = 66589, # mid-2018 pop in Blantyre
                   `South Africa`  = 6907974, # projected from 2011 to mid-2016 using Thembisa HIV model 
                   `Brazil`        = 27734412) %>% # projected from 2011 to mid-2016 per world popn growth stats
  map_df(.id = "country", ~data.frame(N = .x))

#set smooth or unsmooth conditions
pop_use_totals <- FALSE
pop_smooth    <- TRUE

#unsmoothed population values
pop_country_df <- list(`England` = pop_en,
                       `Malawi`        = pop_mw,
                       `South Africa`  = pop_sa,
                       `Brazil`        = pop_br) %>%
  bind_rows(.id = "country") %>%
  group_by(country) %>%
  inner_join(pop_totals) %>%
  mutate(p = ntotal/(pop_use_totals*N + (1-pop_use_totals)*sum(ntotal)))

#smoothed population values
if (pop_smooth){
  
  pop_country_df %<>%
    split(.$country) %>%
    map(~mutate(.x, N_ = sum(ntotal),
                ntotal = c(head(.x$ntotal, 2),
                           zoo::rollmean(.x$ntotal, 5),
                           tail(.x$ntotal, 2)))) %>%
    map_df(~mutate(.x,
                   p = ntotal/sum(ntotal)))
  
}

pop_cases <- dplyr::left_join(ipd_curves, pop_country_df, by = c("agey", "country")) %>% 
  dplyr::mutate(cases  = `50%`/scale*ntotal,
                lcases = `2.5%`/scale*ntotal, 
                ucases = `97.5%`/scale*ntotal, Vac.age = agey) %>% 
  filter(serogroup != "All") 
