################################################################################
###                   Formatting & Tidying METADATA                       #####
###############################################################################
## - Alejandro, I've copied this from mine - you might have to update to match your
## metadata sheet (labsheet)

rm(list=ls()) # Clear memory

## Load Libraries ----

install.packages('remotes')
library('remotes')
options(timeout=9999999)
remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(tidyverse)
library(here)


#set study name
name <- "Baitcomp_All" ##maybe change this and add _allspecies or something 


#metadata (labsheet)
metadata <- read_metadata(here::here("./data/raw/em export"), #update path
                          method = "BRUVs") %>%
  dplyr::select(-c(campaignid, comment, maxn_checker, system, observer_count,
                   observer_length)) %>%
  dplyr::mutate(date_time = mdy_hm(date_time, tz = "GMT")) %>% 
  dplyr::mutate(date_time = with_tz(date_time, tzone = "Asia/Singapore"))%>% ##updating date & time
  dplyr::mutate(date_time = format(date_time, "%Y/%m/%dT%H:%M:%S")) %>%
  dplyr::mutate(date = substr(date_time, 1, 10))%>%
  dplyr::mutate(time = substr(date_time, 12, 19))%>%
  dplyr::mutate(time_of_day = as.POSIXct(time, format = "%H:%M:%S"))%>%
  dplyr::mutate(time_sec = as.numeric(format(time_of_day, "%H")) * 3600 +
                  as.numeric(format(time_of_day, "%M")) * 60 +
                  as.numeric(format(time_of_day, "%S")))%>%
  dplyr::mutate(time_hr = as.numeric(time_sec/3600))%>%
  filter(!is.na(opcode)data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC) %>% 
  dplyr::mutate(depth_m = as.integer(depth_m))%>% ## TODO
  glimpse()
## PUT 1 IN THE LAST LINE OF THE CONSOLE (Just the first time)
sum(metadata$successful_count == "Yes")


## Checking to see if depth is same in field metadata as final metadata
field_meta <- read.csv("./data/raw/metadata_field.csv")%>%
  clean_names()%>%
  dplyr::rename(depth_m=depth_m_)%>%
  dplyr::select(sample, depth_m, notes, x, y)%>%
  dplyr::mutate(opcode = sprintf("%03d", sample))%>%
  glimpse()

depth_mismatches <- field_meta %>%
  dplyr::select(opcode, depth_m) %>%
  anti_join(
    metadata %>% select(opcode, depth_m),
    by = c("opcode", "depth_m"))


# saving metadata as RDS
saveRDS(metadata, file = here::here(paste0("./data/tidy/",
                                           name, "_Metadata.rds")))
