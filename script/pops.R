
#loaad datasets into memory
pop.ew <- read_csv(here("data", "EW_total_pop.csv"))
pop.un <- read_csv(here("data", "UN_total_pop.csv"))

#smooth curve for UN regional populations to get annual estimates
ssa <- subset(pop.un, region=="Sub-Saharan Africa")
naw <- subset(pop.un, region=="Northen Africa & Western Asia")
csa <- subset(pop.un, region=="Central & Southern Asia")
esa <- subset(pop.un, region=="Eastern & South-Eastern Asia")
lac <- subset(pop.un, region=="Latin America & the Caribbean")
anz <- subset(pop.un, region=="Australia & New Zealand")
oce <- subset(pop.un, region=="Ocenian")
ena <- subset(pop.un, region=="Europe & Northern America")





#fit for different UN pops
unpops <- list(ssa, naw, csa, esa, lac, anz, oce, ena)
pred <- list()
plot <- list()

for(i in 1:8){
theta0 <- min(unpops[[i]]$ntotal,na.rm=TRUE)*0.5  
model0 <- lm(log(ntotal-theta0) ~ sage, data=unpops[[i]])  
alpha0 <- (coef(model0)[1])
beta0  <- (coef(model0)[2])

model <- nls(data=unpops[[i]], trace=TRUE, algorithm="port", log(ntotal) ~ alpha*exp(beta*sage), 
                 nls.control(maxiter=200), start=list(beta=-(-beta0), alpha=alpha0))

pred[[i]] <- data.frame(sage = seq(55, 90, by = 1)) %>% mutate(ntotal = exp(predict(model, newdata = .)))
plot[[i]] <- ggplot(data=unpops[[i]]) + geom_point(aes(x=sage, y=ntotal)) + geom_line(data = pred[[i]], aes(x=sage, y=ntotal))
}

#fit for Sub-Saharan Africa
theta0 <- min(ssa$ntotal,na.rm=TRUE)*0.5  
model0 <- lm(log(ntotal-theta0) ~ sage, data=ssa)  
alpha0 <- (coef(model0)[1])
beta0  <- (coef(model0)[2])

ssa_model <- nls(data=ssa, trace=TRUE, algorithm="port", log(ntotal) ~ alpha*exp(beta*sage), 
                 nls.control(maxiter=200), start=list(beta=-(-beta0), alpha=alpha0))

ssa_pred <- data.frame(sage = seq(55, 90, by = 1)) %>% mutate(ntotal = exp(predict(ssa_model, newdata = .)))
ggplot(data=ssa) + geom_point(aes(x=sage, y=ntotal)) + geom_line(data = ssa_pred, aes(x=sage, y=ntotal))

#ggplot comparing % populations in England/Wales versus UN SDG regions
countries <- c("England/Wales"=rgb(205, 17, 39, maxColorValue = 255),
               "Malawi"="black")

pop.totals <- list(`England/Wales` = 56286961 + 3152879, # mid-2019
                   `Malawi`        = 18628747) %>%
  map_df(.id = "Country", ~data.frame(N = .x))

use.pop.totals <- FALSE
smooth.pops    <- TRUE

countries_df <- list(`England/Wales` = pop.ew,
                     `Malawi`        = pop.mw) %>%
  bind_rows(.id = "Country") %>%
  group_by(Country) %>%
  inner_join(pop.totals) %>%
  mutate(p = ntotal/(use.pop.totals*N + (1-use.pop.totals)*sum(ntotal)))

if (smooth.pops){
  
  countries_df %<>%
    split(.$Country) %>%
    map(~mutate(.x, N_ = sum(ntotal),
                ntotal = c(head(.x$ntotal, 2),
                           zoo::rollmean(.x$ntotal, 5),
                           tail(.x$ntotal, 2)))) %>%
    map_df(~mutate(.x,
                   p = ntotal/sum(ntotal)))
  
}

countries_plot <- ggplot(data = countries_df,
       aes(x = agey,
           y = p)) +
  geom_col(aes(fill = Country),
           position = position_dodge()) +
  scale_x_continuous(breaks  = seq(50, 90, by = 10),
                     labels  = function(x){gsub(pattern = "90",
                                                replacement = "90+",
                                                x = x)}) + 
  scale_y_continuous(labels  = scales::percent,
                     limits = c(0, NA)) + 
  theme_bw() + 
  labs(title="", x="Age (years)", 
       y=paste0("Age as share of\n",
               ifelse(use.pop.totals,"total national","55+" ),
               " population"),
       color="Countries") +
  scale_fill_manual(values=countries) +
  theme(axis.text=element_text(face="bold", size=10, color="black"),
        legend.position = "bottom")

ggsave(filename = here("output","countries.png"), 
       plot = countries_plot,
       width = 7, height = 3.5, units = "in", dpi = 300)
