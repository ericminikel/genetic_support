library(tidyverse)
library(writexl)
library(readxl)

target_indication_highest_status_after_2000 <- read_xlsx("./../input_sources/target_indication_highest_status_after_2000.xlsx")
target_indication_current_status_after_2000 <- read_xlsx("./../input_sources/target_indication_current_status_after_2000.xlsx")
target_indication_single_target_highest_status_after_2000 <- read_xlsx("./../input_sources/target_indication_single_target_highest_status_after_2000.xlsx")
target_indication_single_target_current_status_after_2000 <- read_xlsx("./../input_sources/target_indication_single_target_current_status_after_2000.xlsx")
target_indication_with_orphan_drugs <- read_xlsx("./../input_sources/target_indication_with_orphan_drugs.xlsx")

target_indication_highest_status_after_2000 %>% 
  left_join(target_indication_with_orphan_drugs, by = c("gene", "mesh_id")) %>% 
  replace_na(
    list(
      is_orphan_drug = FALSE
    )
  ) %>% 
  left_join(target_indication_current_status_after_2000, by = c("gene", "mesh_id")) %>% 
  left_join(target_indication_single_target_highest_status_after_2000, by = c("gene", "mesh_id")) %>% 
  left_join(target_indication_single_target_current_status_after_2000, by = c("gene", "mesh_id")) -> target_indication_pair_with_all_info

## export
target_indication_pair_with_all_info %>% 
  write_xlsx("./../input_sources/target_indication_pair_with_all_info.xlsx")

## all together without year filter
target_indication_highest_status_all_years <- read_xlsx("./../input_sources/target_indication_highest_status_after_2000.xlsx")
target_indication_current_status_all_years <- read_xlsx("./../input_sources/target_indication_current_status_after_2000.xlsx")
target_indication_single_target_highest_status_all_years <- read_xlsx("./../input_sources/target_indication_single_target_highest_status_after_2000.xlsx")
target_indication_single_target_current_status_all_years <- read_xlsx("./../input_sources/target_indication_single_target_current_status_all_years.xlsx")
target_indication_with_orphan_drugs <- read_xlsx("./../input_sources/target_indication_with_orphan_drugs.xlsx")

## export
target_indication_highest_status_all_years %>% 
  left_join(target_indication_with_orphan_drugs_all_years, by = c("gene", "mesh_id")) %>% 
  replace_na(
    list(
      is_orphan_drug = FALSE
    )
  ) %>% 
  left_join(target_indication_current_status_all_years, by = c("gene", "mesh_id")) %>% 
  left_join(target_indication_single_target_highest_status_all_years, by = c("gene", "mesh_id")) %>% 
  left_join(target_indication_single_target_current_status_all_years, by = c("gene", "mesh_id"))  -> target_indication_pair_with_all_info_all_years

## export
target_indication_pair_with_all_info_all_years %>% 
  write_xlsx("./../input_sources/target_indication_pair_with_all_info.xlsx")

