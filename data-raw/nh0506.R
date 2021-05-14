## Code to prepare `nh0506Homocysteine` dataset

# These data files are not included and come from
# `https://wwwn.cdc.gov/nchs/nhanes/ContinuousNhanes/Default.aspx?BeginYear=2005`
DEMO <- foreign::read.xport("../data_example3/DEMO_D.XPT")
HCY <- foreign::read.xport("../data_example3/HCY_D.XPT")
SMQ <- foreign::read.xport("../data_example3/SMQ_D.XPT")
BMX <- foreign::read.xport("../data_example3/BMX_D.XPT")
COT <- foreign::read.xport("../data_example3/COT_D.XPT")

d <- merge(DEMO, HCY, by="SEQN", all.x=TRUE)
d <- merge(d, SMQ, by="SEQN", all.x=TRUE)
d <- merge(d, COT, by="SEQN", all.x=TRUE)
d <- merge(d, BMX, by="SEQN", all.x=TRUE)

rm(DEMO, HCY, SMQ, COT, BMX)

SEQN <- d$SEQN
age <- d$RIDAGEYR
race <- factor(d$RIDRETH1, levels = 1:5,
               labels = c("Mexican American",
                          "Other Hispanic",
                          "Non-Hispanic White",
                          "Non-Hispanic Black",
                          "Other Race - Including Multi-Racial"))
sex <- factor(d$RIAGENDR, levels = 1:2, labels = c("Male", "Female"))
education <- factor(d$DMDEDUC2, levels = c(1:5, 7, 9, NA),
                    labels = c("< Grade 9", "9-11th grade",
                               "High school grad/GED",
                               "Some college or AA degree",
                               "College graduate or above",
                               rep("Unknown education", 3)),
                    exclude = NULL)
povertyr <- d$INDFMPIR
homocysteine <- d$LBXHCY
bmi <- d$BMXBMI
cotinine <- d$LBXCOT

cigs100life <- d$SMQ020 == 1
cigs100life[d$SMQ020 > 3] <- NA
smokenow <- d$SMQ040 < 2.5
smokenow[!cigs100life] <- FALSE
cigsdays30 <- d$SMD641
cigsdays30[cigsdays30 > 32] <- NA
cigsdays30[!smokenow] <- 0
cigsperday30 <- d$SMD650
cigsperday30[cigsperday30 > 100] <- NA
cigsperday30[!smokenow] <- 0

dailysmoker <- cigs100life & (cigsdays30==30) & smokenow & (cigsperday30 >= 10)
neversmoker <- (!cigs100life) & (!smokenow)
z <- dailysmoker
z[(!neversmoker) & (!dailysmoker)] <- NA
z <- as.numeric(z)

ds <- data.frame(SEQN, z, sex, age, race, education, povertyr, bmi,
                 cigsperday30, cotinine, homocysteine)
nh0506 <- ds[age >= 20 & !is.na(z) & !is.na(homocysteine) &
                           !is.na(cotinine), ]
# As there is only one person remaining with unknown education at this step, remove them
nh0506 <- nh0506[nh0506$education != "Unknown education", ]
nh0506$education <- factor(nh0506$education, levels = levels(nh0506$education),
                           exclude = "Unknown education")

rm(SEQN, sex, age, education, povertyr, homocysteine,
   cotinine,  cigs100life, smokenow, cigsdays30, cigsperday30,  dailysmoker,
   neversmoker, z, race, bmi, d, ds)

usethis::use_data(nh0506, overwrite = TRUE)
