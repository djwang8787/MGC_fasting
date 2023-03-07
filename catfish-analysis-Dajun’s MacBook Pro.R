pacman::p_load(readr, dplyr, lubridate, lme4, MuMIn, season, gridExtra, Kendall)

# To print/save figures for manuscript use
make_tall_fig = function(object){
  filename = paste0(as.character(substitute(object)), ".png")
  png.test = png(filename = as.character(filename), width = 19.05, height = 22.23, units = 'cm', res = 300)
  gridExtra::grid.arrange(object, newpage = TRUE)
  dev.off()
}

make_wide_fig = function(object){
  filename = paste0(as.character(substitute(object)), ".png")
  png.test = png(filename = as.character(filename), width = 22.23, height = 11.10, units = 'cm', res = 300)
  gridExtra::grid.arrange(object, newpage = TRUE)
  dev.off()
}

#### Rain data 2022 ####
rf.2022 = read.csv("Rainfall/Rainfall-2022/combined-2022-data.csv")

rf.2022 %>%
  group_by(Month) %>%
  summarise(mean(Daily.Rainfall.Total..mm.)) ## Need to compile this with the 2022 dataset; Lianne will be compiling.

#### 1. Activity budget; 2019--2021 ####
data = read.csv("catfish.csv")
data = data %>%
  dplyr::select(Date, Eating, Not.eating, Normally, Erratically, Dark, Bright, Rainfall)

data = data %>%
  dplyr::mutate(Date = as.Date(Date, format = "%m/%d/%Y"),
                Year = as.factor(year(Date)),
                Week = as.factor(week(Date)),
                Month = as.factor(month(Date)),
                Eat.prop = Eating/6,
                Mek.Monsoon = ifelse(Month %in% c(6:11), 'Southwest', "Northeast"),
                SG.Monsoon = ifelse(Month %in% c(6:9), 'Southwest',
                                     ifelse(Month %in% c(12:3), "Northeast", "Transition")))

# Correlation
data %>%
  filter(Year == "2021") %>%
  dplyr::summarise(fast.color = cor(Eating,Dark),
                   fast.swim = cor(Eating,Normally))

correl.data = data %>%
  filter(Year == "2021")

ggplot(correl.data, aes(x = 1:12)) +
  geom_line(aes(y = Normally/6), color = "blue", size = 1., linetype = 2) +
  geom_line(aes(y = Dark/6), color = "red", size = 1., linetype = 3) +
  geom_line(aes(y = Eat.prop), color = "black", size = 1., linetype = 1) +
  geom_point(aes(y = Eat.prop), color = "black", size = 3.5) +
  geom_point(aes(y = Eat.prop), color = "white", size = 3.0) +
  scale_x_discrete(name ="Month", limits=c(1:12)) +
  ylab("Odds of catfish entering fasting cycle") +
  theme_minimal()

MannKendall(correl.data$Eating)
MannKendall(correl.data$Normally)
MannKendall(correl.data$Dark)

fig3 = ggpubr::ggarrange(
ggplot(correl.data, aes(x = 1:12, y = Normally/6)) +
  geom_bar(stat = "identity") +
  coord_polar() +
  scale_x_continuous(breaks = c(1:12),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May",
                                "Jun", "Jul", "Aug", "Sep", "Oct",
                                "Nov", "Dec")) +
  theme_minimal() +
  xlab(" ") +
  ylab("Proportion (%)")
,
ggplot(correl.data, aes(x = 1:12, y = Dark/6)) +
  geom_bar(stat = "identity") +
  coord_polar() +
  scale_x_continuous(breaks = c(1:12),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May",
                                "Jun", "Jul", "Aug", "Sep", "Oct",
                                "Nov", "Dec")) +
  theme_minimal() +
  xlab("Months") +
  ylab(" ") +
  theme(axis.text.y = element_blank())
,
ggplot(correl.data, aes(x = 1:12)) +
  geom_bar(aes(y = Eating/6), stat = "identity") +
  coord_polar() +
  scale_x_continuous(breaks = c(1:12),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May",
                                "Jun", "Jul", "Aug", "Sep", "Oct",
                                "Nov", "Dec")) +
  theme_minimal() +
  xlab(" ") +
  ylab(" ") +
  theme(axis.text.y = element_blank()),
nrow = 1,
labels = c("A", "B", "C")
)

make_wide_fig(fig3)

# To wait for 2022 data to complete; remind Lianne.

#### 2. Periodicity x feeding response data ####
# the projected periodicity / rhythmicity for years 2018--2021

data.monthglm = data %>%
  dplyr::mutate(weight = n(),
                month = as.numeric(Month),
                year = as.numeric(Year)) %>%
  dplyr::select(Eat.prop, month, year)


mmodel = monthglm(formula=Eat.prop~1,
                  data=data.monthglm,
                    family=gaussian(),
                  offsetmonth=TRUE,
                  refmonth=11)
summary(mmodel)
plot(mmodel, ylim=c(0,1))


data.monthglm2 = data %>%
  dplyr::mutate(weight = n(),
                month = as.numeric(Month),
                year = ifelse(Year == "2020", as.numeric("2020"),
                              ifelse(Year == "2021", as.numeric("2021"),
                                     ifelse(Year == "2018", as.numeric("2018"), as.numeric("2019"))))) %>%
  dplyr::select(Eat.prop, month, year, Rainfall, Mek.Monsoon, SG.Monsoon)


mmodel1 = monthglm(formula=Eat.prop*6~1, data=data.monthglm2,
                  family=poisson(),  offsetmonth=TRUE, refmonth=11) #AIC 158.5
mm1.aic = 158.5
mm1.aicc = mm1.aic - (2*1*(1+1)/(39-1-1))
mm1.ll = (2*(1) - mm1.aicc)/2
mm1.rl = exp(-0.5 * 0) # relative likelihood: 1
mm1.rl / (mm1.rl + exp(-0.5 * 0.2316081) + exp(-0.5 * 0.6899081)) # weight: 0.3847778



mmodel2 = monthglm(formula=Eat.prop*6~year, data=data.monthglm2,
                   family=poisson(),  offsetmonth=TRUE, refmonth=11) #AIC 159.8
159.8 - (2*4*(4+1)/(39-4-1)) # AIcc 158.6235
(2*(4) - 159.8)/2 # LL -75.9
158.6235 - mm1.aicc # delta; 0.2316081
exp(-0.5 * 0.2316081) # relative LL, 0.8906497
exp(-0.5 * 0.2316081) / (1 + exp(-0.5 * 0.2316081) + exp(-0.5 * 0.6899081)) # weight: 0.3427022


mmodel3 = monthglm(formula=Eat.prop*6~year + Rainfall, data=data.monthglm2,
                   family=poisson(),  offsetmonth=TRUE, refmonth=11) #AIC 159.0818
160.9 - (2*5*(5+1)/(39-5-1)) # AIcc 159.0818
(2*(5) - 159.0818)/2 # LL -74.5409
159.0818 - mm1.aicc #delta; 0.6899081
exp(-0.5 * 0.6899081) # relative LL, 0.7082529
exp(-0.5 * 0.6899081) / (1 + exp(-0.5 * 0.2316081) + exp(-0.5 * 0.6899081)) # weight: 0.27252



df1 = data.frame(unclass(summary(mmodel1)))

to.plot = data.frame(
  odds = exp(df1$month.ests.mean) / (1+exp(df1$month.ests.mean)),
  upper = exp(df1$month.ests.upper) / (1+exp(df1$month.ests.upper)),
  lower = exp(df1$month.ests.lower) / (1+exp(df1$month.ests.lower))
)

to.plot$upper.ci = to.plot$upper - to.plot$odds
to.plot$lower.ci = to.plot$odds - to.plot$lower

to.plot = to.plot %>% add_row(odds = 0.05614, upper = 0, lower = 0,
                              upper.ci = 0, lower.ci = 0, .before = 11)

to.plot$a = data.monthglm2$Eat.prop[4:15]
to.plot$b = data.monthglm2$Eat.prop[16:27]
to.plot$c = data.monthglm2$Eat.prop[28:39]


fig1 = ggplot(to.plot, aes(x = 1:12)) +
  geom_errorbar(aes(
    ymin = odds - lower.ci,
    ymax = odds + upper.ci),
    width = 0.1) +
  geom_line(aes(y = odds), color = "black", size = 1.5) +
  geom_point(aes(y = odds), color = "black", size = 3.5) +
  geom_point(aes(y = odds), color = "white", size = 3.0) +
  scale_x_discrete(name ="Month", limits=c(1:12)) +
  ylab("Odds of catfish entering fasting cycle") +
  theme_minimal()

# make_wide_fig(fig1)

fig2 = ggplot(to.plot, aes(x = 1:12)) +
  geom_errorbar(aes(
    ymin = odds - lower.ci,
    ymax = odds + upper.ci),
    width = 0.1) +
  geom_point(aes(y = odds), color = "black", size = 3.5) +
  geom_point(aes(y = odds), color = "white", size = 3.0) +
  geom_line(aes(y = a), color = "red", size = 1.0, linetype = 4) +
  geom_line(aes(y = b), color = "blue", size = 1.0, linetype = 2) +
  geom_line(aes(y = c), color = "green", size = 1.5, linetype = 3) +
  scale_x_discrete(name ="Month", limits=c(1:12)) +
  ylab("Odds of catfish entering fasting cycle") +
  theme_minimal()
# make_wide_fig(fig2)

data.monthglm2 %>%
  filter(!year == "2018") %>%
  ggplot(aes(x=1:36)) +
  geom_line(aes(y = Eat.prop*100)) +
  geom_point(aes(y = Eat.prop*100)) +
  scale_x_discrete(name ="Month", limits=c(1:36)) +
  ylab("Proportion of fasting catfish (%)") +
  theme_minimal()

#### Mann-Kendall test ####
acf(data.monthglm2$Eat.prop)
pacf(data.monthglm2$Eat.prop)

MannKendall(data.monthglm2$Rainfall)

MK1 = data.monthglm2 %>%
  filter(year == "2020")
MannKendall(MK1$Rainfall)

#### Defunct ####
#### 2.1a Periodicity x feeding response data ####
# # This assumes multicosinor models which examines across 12 months
# # not suitable as it underfits the 12 months in comparison to monthglm()
#
# data.rm = data %>%
#   group_by(Year) %>%
#   dplyr::mutate(weight = n(),
#                 Month.integer = as.numeric(Month))
#
# rm1 = glmer(Eat.prop ~ 1 + (1|Year),  #null
#             data = data.rm,
#             weights = weight,
#             family = binomial(link='logit'))
#
# rm2 = glmer(Eat.prop ~ cos(2*pi*Month.integer/12) + (1|Year),
#             data = data.rm, # uni-modal; one peak
#             weights = weight,
#             family = binomial(link='logit'))
#
# rm3 = glmer(Eat.prop ~ cos(2*pi*Month.integer/6) + (1|Year),
#             data = data.rm, # bi-modal; two peaks
#             weights = weight,
#             family = binomial(link='logit'))
#
# rm4 = glmer(Eat.prop ~  cos(2*pi*Month.integer/4) + (1|Year),
#             data = data.rm, # tri-modal; three peaks; or cathemeral
#             weights = weight,
#             family = binomial(link='logit'))
#
#
# rm.output = model.sel(rm1, rm2, rm3, rm4)
# rm.output # uni-modal has the best fit.
# hour.sel = rep(1:12, times = 1)
#
# rm2.cr = rm2@beta[1] +
#   rm2@beta[2]*cos(2*pi*hour.sel/12)
# plot(rm2.cr)
#

#### 1.1a fine-scale 2022 data ####
data2022 = read.csv("Mekong_catfish_binary_data_2022.csv")

data2022 = data2022[1:1944,1:8] # Take away all remarks

data2022$Date = as.POSIXct(data2022$Date, format = "%d/%m/%y")
head(data2022)

data2022_cleaned = data2022 %>%
  mutate(month = lubridate::month(Date),
         day = lubridate::day(Date)) %>%
  filter(Feeding_Time == "AM") %>%
  mutate(Feed = ifelse(Pellet_Cake %in% c(1), 1,
                       ifelse(Fish_Prawn %in% c(1), 1, 0))) %>%
  dplyr::select(-Rubbing) %>%
  na.omit()


data2022_cleaned %>%
  group_by(month, day) %>%
  mutate(total.feed = sum(Feed),
         total.swim = sum(Swimming_Eratically),
         total) %>%
  head()


data2022[is.na(data2022)] = 0


##### 1.2 Fine-scale ####
data2 = read.csv("mgc-dataset-2021.csv")
data2$Date = data2[,1]
data2 = data2 %>%
  select(Date, Eating, Not.eating, Normally, Erratically, Dark, Bright, Rainfall)

data2 = data2 %>%
  dplyr::mutate(Eat.prop = Eating/6,
                Swim.prop = Normally/6,
                Colour = Dark/6,
                Date = as.Date(Date, format = "%m/%d/%Y"),
                Year = as.factor(year(Date)),
                Week = as.factor(week(Date)),
                Month = as.factor(month(Date)),
                Day = yday(Date),
                Mek.Monsoon = ifelse(Month %in% c(6:11), 'Southwest', "Northeast"),
                SG.Monsoon = ifelse(Month %in% c(6:9), 'Southwest',
                                    ifelse(Month %in% c(12:3), "Northeast", "Transition")))

# generalised linear models
m0 = glm(Eat.prop ~ 1,
         data = data2,
         family = binomial) # Null model

m1 = glm(Eat.prop ~ SG.Monsoon,
         data = data2,
         family = binomial)

m2 = glm(Eat.prop ~ Mek.Monsoon,
         data = data2,
         family = binomial)

m3 = glm(Eat.prop ~  Rainfall,
         data = data2,
         family = binomial)

m4 = glm(Eat.prop ~ SG.Monsoon + Rainfall,
         data = data2,
         family = binomial)

m5 = glm(Eat.prop ~ Mek.Monsoon + Rainfall,
         data = data2,
         family = binomial)

model.sel(m0, m1, m2, m3, m4, m5)

