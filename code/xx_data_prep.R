# EDSD 2023-2024
# Analysis of mortality disturbances course
# Instructor: Enrique Acosta (CED)
# Preparing environment

# installing missing packages ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
libs <- c("devtools",
          "tidyverse", 
          "lubridate",
          "remotes",
          "readr",
          "viridisLite",
          "viridis", 
          "mgcv",
          "HMDHFDplus",
          "ISOweek",
          "countrycode",
          "patchwork",
          "isoband",
          "RColorBrewer")

for (i in libs){
  if (!(i %in% rownames(installed.packages()))) install.packages(i)
}
options(timeout = 600)
if (!("MortalitySmooth" %in% rownames(installed.packages()))) remotes::install_github("kikeacosta/MortalitySmooth")
if (!("wpp2022" %in% rownames(installed.packages()))) remotes::install_github("PPgp/wpp2022", force = TRUE)

# Loading required packages 
lapply(libs, require, character.only = T)
library("MortalitySmooth")
library("wpp2022")

# avoiding scientific notation
options(scipen=999)
# let's keep the same seed for reproducibility of results
set.seed(2019) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# population data from WPP 2022 ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# loading population data from WPP 2022 estimates

if (!file.exists("data_input/wpp2022_pop.rds")) {
  
  data(pop1dt)
  data(popproj1dt)
  
  pop_hist <- 
    pop1dt %>% 
    as_tibble() %>% 
    # filter(name =='Spain') %>% 
    select(name, year, pop, popM, popF)
  
  pop_proj <- 
    popproj1dt %>% 
    as_tibble() %>% 
    # filter(name =='Spain') %>% 
    select(name, year, pop, popM, popF)
  
  pop <- 
    bind_rows(pop_hist, pop_proj) %>% 
    arrange(name, year) %>% 
    mutate(code = countrycode(name, origin = "country.name",
                              destination = "iso3c")) %>% 
    drop_na(code)
  
  # write_rds(pop, "data_input/total_annual_population_all_countries.rds")
  write_rds(pop, "data_input/wpp2022_pop.rds")
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deaths and exposures from the HMD ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (!file.exists("data_input/hmd_dts_pop.rds")) {
  
  # HMD user and password
  # usethis::edit_r_environ()
  # hmd_us="acosta@demogr.mpg.de"
  # hmd_pw="Secreto_1"
  
  # getting HMD username and password from the R environment
  hmd_us <- Sys.getenv("hmd_us")
  hmd_pw <- Sys.getenv("hmd_pw")
  
  
  cat("Downloading deaths and exposures from HMD\n")
  cds_hmd <- getHMDcountries() %>% pull(CNTRY)
  
  hmd <- tibble()
  for(ct in cds_hmd){
    cat(paste0(ct, "\n"))
    chunk_d <- 
      readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
      as_tibble() %>%
      mutate(Code = ct)
    
    hmd <- 
      hmd %>%
      bind_rows(chunk_d)
  }
  
  hmd
  
  hmd_e <- tibble()
  for(ct in cds_hmd){
    cat(paste0(ct, "\n"))
    chunk_e <- 
      readHMDweb(ct, "Exposures_1x1", hmd_us, hmd_pw) %>%
      as_tibble() %>%
      mutate(Code = ct)
    
    hmd_e <- 
      hmd_e %>%
      bind_rows(chunk_e)
  }
  
  hmd2 <- 
    hmd %>% 
    rename_all(tolower) %>% 
    select(-openinterval) %>% 
    gather(female, male, total, key = sex, value = dts)
  
  hmd_e2 <- 
    hmd_e %>% 
    rename_all(tolower) %>% 
    select(-openinterval) %>% 
    gather(female, male, total, key = sex, value = pop)
  
  hmd_all <- 
    hmd2 %>% 
    left_join(hmd_e2) 
  
  write_rds(hmd_all, "data_input/hmd_dts_pop.rds",
            compress = "xz")
  
}

# other HMD data?
# New Zealand Maories and non Maories 
# Germany East and West

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deaths and exposures from the full HMD ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!file.exists("data_input/hmd_dts_pop_vcomplement.rds")) {
  
  dts_files <- 
    list.files("data_input/Deaths_1x1") %>% 
    as_tibble() %>% 
    filter(str_detect(value, "GBR_SCO|DEUTE|DEUTW|NZL_MA|NZL_NM")) %>% 
    pull(value)
  
  pop_files <- 
    list.files("data_input/Exposures_1x1") %>% 
    as_tibble() %>% 
    filter(str_detect(value, "GBR_SCO|DEUTE|DEUTW|NZL_MA|NZL_NM")) %>% 
    pull(value)
  
  # loading deaths from all countries in HMD
  db_d <- tibble()
  for(i in 1:length(dts_files)){
    
    txt_file <- dts_files[i]
    
    print(txt_file)
    
    temp <- 
      read_tsv(paste0("data_input/Deaths_1x1/", txt_file),
               skip = 1) %>%
      rename(var =1) %>% 
      mutate(var = str_replace_all(var, "\\.", "")) %>% 
      separate(1, c("year", "age", "female", "male", "total")) %>% 
      mutate(code = str_replace(txt_file, ".Deaths_1x1.txt", ""))
    
    db_d <- db_d %>% 
      bind_rows(temp)
  }
  
  # loading exposures from all countries in HMD
  db_p <- tibble()
  for(i in 1:length(pop_files)){
    
    txt_file <- pop_files[i]
    
    print(txt_file)
    
    temp <- 
      read_tsv(paste0("data_input/Exposures_1x1/", txt_file),
               skip = 1) %>%
      rename(var =1) %>% 
      mutate(var = str_replace_all(var, "\\.", "")) %>% 
      separate(1, c("year", "age", "female", "male", "total")) %>% 
      mutate(code = str_replace(txt_file, ".Exposures_1x1.txt", ""))
    
    db_p <- db_p %>% 
      bind_rows(temp)
  }
  
  db_d2 <- 
    db_d %>% 
    gather(male, female, total, key = sex, value = dts) %>% 
    left_join(db_p %>% 
                gather(male, female, total, key = sex, value = pop)) %>% 
    mutate(dts = dts %>% as.double(),
           pop = pop %>% as.double(),
           age = age %>% as.integer(),
           year = year %>% as.integer())
  
  write_rds(db_d2, "data_input/hmd_dts_pop_vcomplement.rds")
}

if (!file.exists("data_input/hmd_dts_pop_v2.rds")) {
  
  hmd1 <- read_rds("data_input/hmd_dts_pop.rds")
  hmd2 <- read_rds("data_input/hmd_dts_pop_vcomplement.rds")
  
  out <-  
    bind_rows(hmd1, hmd2)
  
  write_rds(out, "data_input/hmd_dts_pop_v2.rds",
            compress = "xz")
  
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# weekly mortality from the STMF ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (!file.exists("data_input/STMFinput.zip")) {
  download.file("https://www.mortality.org/File/GetDocument/Public/STMF/Inputs/STMFinput.zip",
                "data_input/STMFinput.zip")
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Age-specific fertility rates from the HFD ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (!file.exists("data_input/hfd_asfr.rds")) {
  
  hmd_us <- Sys.getenv("hmd_us")
  hmd_pw <- Sys.getenv("hmd_pw")
  
  # some countries having issues when trying to download... 
  # excluding them for now
  cts_issues <- c("BGR", "CAN", "POL", "KOR", "RUS", "CHE")
  cts_issues <- c()
  
  cat("Downloading ASFR from HFD\n")
  
  cds_hfd <- 
    getHFDcountries() %>% 
    filter(!CNTRY %in% cts_issues) %>% 
    pull(CNTRY)
  
  hfd <- tibble()
  for(ct in cds_hfd[c(1:3, 6:43)]){
    cat(paste0(ct, "\n"))
    
    chunk_f <- 
      readHFDweb(ct, "asfrRR", hmd_us, hmd_pw) %>%
      as_tibble() %>%
      mutate(Code = ct)
    
    hfd <- 
      hfd %>%
      bind_rows(chunk_f) %>% 
      unique()
  }
  hfd
  
  hfd2 <- 
    hfd %>% 
    rename_all(tolower) %>% 
    select(-openinterval)
  
  write_rds(hfd2, "data_input/hfd_asfr.rds")
  
}

