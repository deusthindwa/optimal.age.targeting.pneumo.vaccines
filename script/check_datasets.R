# check proportions of IPD cases across years

rsa <- readr::read_csv(here("data", "rsa.csv"))
bra <- readr::read_csv(here("data", "bra.csv"))

A <- bra  %>% 
  mutate(year = as.factor(year), age = if_else(age >= 90, 90, age)) %>%
  ggplot(aes(x = age, y = pcv13, color = year)) +
  geom_line() +
  labs(title = "BRA", x = "Age,y", y = "Cases by PCV13 ST") +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  ylim(0,16) +
  theme_bw() +
  theme(legend.position = "none")

B <- bra  %>% 
  mutate(year = as.factor(year), age = if_else(age >= 90, 90, age)) %>%
  ggplot(aes(x = age, y = ppv23, color = year)) +
  geom_line() +
  labs(title = "BRA", x = "Age,y", y = "Cases by PPV23 ST") +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  ylim(0,16) +
  theme_bw() +
  theme(legend.position = "right")

C <- rsa  %>% 
  filter(pcv13 == 1 & hiv != "Positive") %>% 
  group_by(year, age) %>% tally() %>%
  mutate(year = as.factor(year), age = if_else(age >= 90, 90, age)) %>%
  ggplot(aes(x = age, y = n, color = year)) +
  geom_line() +
  labs(title = "RSA", x = "Age,y", y = "Cases by PCV13 ST") +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  ylim(0,16) +
  theme_bw() +
  theme(legend.position = "none")

D <- rsa  %>% 
  filter(ppv23a == 1 & hiv != "Positive") %>% 
  group_by(year, age) %>% tally() %>%
  mutate(year = as.factor(year), age = if_else(age >= 90, 90, age)) %>%
  ggplot(aes(x = age, y = n, color = year)) +
  geom_line() +
  labs(title = "RSA", x = "Age,y", y = "Cases by PPV23 ST") +
  scale_x_continuous(breaks = seq(55, 90, 5)) +
  ylim(0,16) +
  theme_bw() +
  theme(legend.position = "right")

(A | B) / (C | D)


