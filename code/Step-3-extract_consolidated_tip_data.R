# Extract EPI_ISL IDs and Complete Tip Names from Consolidated PARNAS Outputs
# Works from main folder, processes consolidated_parnas_outputs directory
# Creates clean Excel files named {Group}_tip_data.xlsx

# Required packages
required_packages <- c("readr", "dplyr", "writexl", "stringr")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cran.r-project.org")
}

# Load packages
library(readr)
library(dplyr)
library(writexl)
library(stringr)

cat("======================================================================\n")
cat("EPI_ISL AND TIP NAME EXTRACTOR FOR CONSOLIDATED PARNAS OUTPUTS\n")
cat("======================================================================\n\n")

# Configuration
CONSOLIDATED_DIR <- "consolidated_parnas_outputs"
OUTPUT_DIR <- "epi_isl_extractions"

# Function to find optimal subtree files in consolidated directory
find_consolidated_subtrees <- function() {
  cat("Searching for optimal subtree files in", CONSOLIDATED_DIR, "...\n")
  
  # Check if consolidated directory exists
  if (!dir.exists(CONSOLIDATED_DIR)) {
    cat("Error: Consolidated directory not found:", CONSOLIDATED_DIR, "\n")
    cat("Make sure you've run the complete pipeline first.\n")
    return(character(0))
  }
  
  # Look for optimal subtree files
  subtree_files <- list.files(
    path = CONSOLIDATED_DIR, 
    pattern = "*_optimal_subtree\\.tre$", 
    full.names = TRUE
  )
  
  if (length(subtree_files) == 0) {
    cat("Warning: No optimal subtree files found in", CONSOLIDATED_DIR, "\n")
    cat("Expected files like: H3_Am_optimal_subtree.tre\n")
    
    # Fallback: look for any .tre files
    all_tre_files <- list.files(
      path = CONSOLIDATED_DIR,
      pattern = "\\.tre$",
      full.names = TRUE
    )
    
    if (length(all_tre_files) > 0) {
      cat("Found", length(all_tre_files), "tree files instead:\n")
      for (f in all_tre_files) {
        cat("  ", basename(f), "\n")
      }
      return(all_tre_files)
    } else {
      return(character(0))
    }
  }
  
  return(subtree_files)
}

# Function to extract tip names from NEXUS format
extract_tip_names_nexus <- function(tree_content) {
  # Look for TAXLABELS section in NEXUS format
  if (str_detect(tree_content, "TAXLABELS")) {
    # Extract everything between TAXLABELS and ;
    taxlabels_section <- str_extract(tree_content, "TAXLABELS[\\s\\S]*?;")
    
    if (!is.na(taxlabels_section)) {
      # Extract quoted names from TAXLABELS
      quoted_names <- str_extract_all(taxlabels_section, "'([^']+)'")[[1]]
      # Remove quotes
      tip_names <- str_remove_all(quoted_names, "'")
      return(tip_names)
    }
  }
  
  return(character(0))
}

# Function to extract tip names from Newick format (in TREES section or standalone)
extract_tip_names_newick <- function(tree_content) {
  # Extract from TREES section if NEXUS, otherwise from whole content
  if (str_detect(tree_content, "BEGIN TREES")) {
    trees_section <- str_extract(tree_content, "BEGIN TREES[\\s\\S]*?END;")
    if (!is.na(trees_section)) {
      tree_content <- trees_section
    }
  }
  
  # Extract quoted names from tree string
  quoted_pattern <- "'([^']+)'"
  quoted_tips <- str_extract_all(tree_content, quoted_pattern)[[1]]
  quoted_tips <- str_remove_all(quoted_tips, "'")
  
  # Also extract unquoted names (for cases without quotes)
  # Pattern to match names before :, ,, or ) that aren't just numbers
  unquoted_pattern <- "([A-Za-z][^:,()\\s']*?)(?=[:,)])"
  unquoted_tips <- str_extract_all(tree_content, unquoted_pattern)[[1]]
  
  # Combine and filter
  all_tips <- c(quoted_tips, unquoted_tips)
  
  # Filter to keep only realistic tip names
  filtered_tips <- all_tips[
    (str_detect(all_tips, "/") | str_detect(all_tips, "EPI_ISL_")) &
    nchar(all_tips) > 10 &
    !str_detect(all_tips, "^[0-9.-]+$")
  ]
  
  return(unique(filtered_tips))
}

# Function to extract tip names from a tree file (handles both NEXUS and Newick)
extract_tip_names_from_tree <- function(tree_file) {
  # Read the tree file
  tree_content <- readr::read_file(tree_file)
  tree_content <- str_trim(tree_content)
  
  tip_names <- character(0)
  
  # Check if it's a NEXUS file
  if (str_detect(tree_content, "#NEXUS|BEGIN TAXA|TAXLABELS")) {
    cat("    Format: NEXUS detected\n")
    
    # Try TAXLABELS first (most reliable for NEXUS)
    tip_names <- extract_tip_names_nexus(tree_content)
    
    # If no TAXLABELS found, extract from tree string
    if (length(tip_names) == 0) {
      tip_names <- extract_tip_names_newick(tree_content)
    }
  } else {
    cat("    Format: Newick detected\n")
    tip_names <- extract_tip_names_newick(tree_content)
  }
  
  return(tip_names)
}

# Function to extract EPI_ISL IDs from tip names
extract_epi_isl_ids <- function(tip_names) {
  # Extract all EPI_ISL_XXXXXXXX patterns
  epi_pattern <- "EPI_ISL_[0-9]+"
  all_epi_ids <- c()
  
  for (tip in tip_names) {
    epi_ids <- str_extract_all(tip, epi_pattern)[[1]]
    all_epi_ids <- c(all_epi_ids, epi_ids)
  }
  
  # Return unique sorted EPI_ISL IDs
  sort(unique(all_epi_ids))
}

# Function to parse tip name components
parse_tip_name <- function(tip_name) {
  # Expected format: A/mallard/USA/023213-008/2025|2025-01-01|EPI_ISL_20176624
  
  # Split by pipe character
  parts <- str_split(tip_name, "\\|")[[1]]
  
  result <- list(
    full_name = tip_name,
    virus_name = ifelse(length(parts) >= 1, parts[1], ""),
    date = ifelse(length(parts) >= 2, parts[2], ""),
    epi_isl = ifelse(length(parts) >= 3, parts[3], "")
  )
  
  # Also try to extract EPI_ISL with regex if not in third part
  if (result$epi_isl == "" || !str_detect(result$epi_isl, "EPI_ISL_")) {
    epi_match <- str_extract(tip_name, "EPI_ISL_[0-9]+")
    result$epi_isl <- ifelse(is.na(epi_match), "", epi_match)
  }
  
  return(result)
}

# Function to extract group name from consolidated filename
get_group_name_from_consolidated <- function(tree_file) {
  filename <- basename(tree_file)
  
  # Remove consolidated filename patterns
  # From: H3_Am_optimal_subtree.tre → H3_Am
  group_name <- filename
  
  # Remove specific consolidated suffixes
  suffixes <- c("_optimal_subtree.tre", "_optimal_colored.tre", "_optimal_clusters.tab",
                "_optimal.tre", "_subtree.tre", ".tre", ".nex", ".nexus")
  
  for (suffix in suffixes) {
    if (str_ends(group_name, suffix)) {
      group_name <- str_remove(group_name, paste0(suffix, "$"))
      break
    }
  }
  
  return(group_name)
}

# Main execution
cat("Starting extraction from consolidated PARNAS outputs...\n\n")

# Step 1: Find consolidated subtree files
tree_files <- find_consolidated_subtrees()

if (length(tree_files) == 0) {
  stop("No tree files found! Make sure the complete pipeline has been run.")
}

cat("Found", length(tree_files), "tree files:\n")
for (f in tree_files) {
  cat("  ", basename(f), "\n")
}
cat("\n")

# Step 2: Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
  cat("Created output directory:", OUTPUT_DIR, "\n")
} else {
  cat("Output directory exists:", OUTPUT_DIR, "\n")
}
cat("\n")

# Step 3: Process each tree file
cat("Processing tree files...\n")

successful_extractions <- 0
total_tips <- 0
total_epi_ids <- 0
all_data <- list()

for (tree_file in tree_files) {
  group_name <- get_group_name_from_consolidated(tree_file)
  cat("Processing:", group_name, "\n")
  cat("  Input: ", basename(tree_file), "\n")
  
  tryCatch({
    # Extract tip names
    tip_names <- extract_tip_names_from_tree(tree_file)
    
    if (length(tip_names) == 0) {
      cat("  ✗ No tip names found!\n\n")
      next
    }
    
    # Extract EPI_ISL IDs
    epi_isl_ids <- extract_epi_isl_ids(tip_names)
    
    cat("  ✓ Extracted", length(tip_names), "tip names\n")
    cat("  ✓ Extracted", length(epi_isl_ids), "EPI_ISL IDs\n")
    
    # Create detailed data frame
    data_rows <- data.frame(
      Tip_Name = character(0),
      Virus_Name = character(0),
      Date = character(0),
      EPI_ISL_ID = character(0),
      Group = character(0),
      Source_File = character(0),
      stringsAsFactors = FALSE
    )
    
    for (tip in sort(tip_names)) {
      # Parse tip name components
      parsed <- parse_tip_name(tip)
      
      data_rows <- rbind(data_rows, data.frame(
        Tip_Name = parsed$full_name,
        Virus_Name = parsed$virus_name,
        Date = parsed$date,
        EPI_ISL_ID = parsed$epi_isl,
        Group = group_name,
        Source_File = basename(tree_file),
        stringsAsFactors = FALSE
      ))
    }
    
    # Create Excel file with clean naming: {Group}_tip_data.xlsx
    excel_file <- file.path(OUTPUT_DIR, paste0(group_name, "_tip_data.xlsx"))
    
    # Prepare sheets
    sheets <- list()
    
    # Sheet 1: Complete detailed data
    sheets$Complete_Data <- data_rows
    
    # Sheet 2: Just EPI_ISL IDs
    if (length(epi_isl_ids) > 0) {
      sheets$EPI_ISL_Only <- data.frame(
        EPI_ISL_ID = epi_isl_ids,
        Group = group_name,
        stringsAsFactors = FALSE
      )
    }
    
    # Sheet 3: Just tip names
    sheets$Tip_Names_Only <- data.frame(
      Tip_Name = sort(tip_names),
      Group = group_name,
      stringsAsFactors = FALSE
    )
    
    # Sheet 4: Summary statistics
    sheets$Summary <- data.frame(
      Metric = c("Total Tips", "Unique EPI_ISL IDs", "Tips with EPI_ISL", "Tips with Date"),
      Count = c(
        length(tip_names),
        length(epi_isl_ids),
        sum(data_rows$EPI_ISL_ID != ""),
        sum(data_rows$Date != "")
      ),
      stringsAsFactors = FALSE
    )
    
    # Write Excel file
    write_xlsx(sheets, excel_file)
    cat("  ✓ Excel file saved:", basename(excel_file), "\n")
    
    # Add to combined data
    all_data[[length(all_data) + 1]] <- data_rows
    
    successful_extractions <- successful_extractions + 1
    total_tips <- total_tips + length(tip_names)
    total_epi_ids <- total_epi_ids + length(epi_isl_ids)
    
    # Show example of tip name
    if (length(tip_names) > 0) {
      cat("  Example tip name:", tip_names[1], "\n")
    }
    
  }, error = function(e) {
    cat("  ✗ Error processing file:", e$message, "\n")
  })
  
  cat("\n")
}

# Step 4: Create combined Excel file
if (length(all_data) > 0) {
  cat("Creating combined Excel file...\n")
  
  # Combine all data
  combined_df <- bind_rows(all_data)
  combined_file <- file.path(OUTPUT_DIR, "all_groups_combined.xlsx")
  
  # Prepare combined sheets
  combined_sheets <- list()
  
  # Sheet 1: All data
  combined_sheets$All_Data <- combined_df
  
  # Sheet 2: Summary by group
  summary_df <- combined_df %>%
    group_by(Group) %>%
    summarise(
      Total_Tips = n(),
      EPI_ISL_Count = sum(EPI_ISL_ID != ""),
      With_Date = sum(Date != ""),
      .groups = 'drop'
    ) %>%
    arrange(desc(Total_Tips))
  
  combined_sheets$Summary_by_Group <- summary_df
  
  # Sheet 3: All unique EPI_ISL IDs
  all_epi_ids <- combined_df %>%
    filter(EPI_ISL_ID != "") %>%
    pull(EPI_ISL_ID) %>%
    unique() %>%
    sort()
  
  if (length(all_epi_ids) > 0) {
    combined_sheets$All_EPI_ISL_IDs <- data.frame(
      EPI_ISL_ID = all_epi_ids,
      stringsAsFactors = FALSE
    )
  }
  
  # Sheet 4: All unique virus names (without EPI_ISL suffix)
  all_virus_names <- combined_df %>%
    filter(Virus_Name != "") %>%
    pull(Virus_Name) %>%
    unique() %>%
    sort()
  
  if (length(all_virus_names) > 0) {
    combined_sheets$All_Virus_Names <- data.frame(
      Virus_Name = all_virus_names,
      stringsAsFactors = FALSE
    )
  }
  
  # Write combined file
  write_xlsx(combined_sheets, combined_file)
  cat("✓ Combined file saved:", basename(combined_file), "\n")
}


# List output files
cat("Output files created:\n")
excel_files <- list.files(OUTPUT_DIR, pattern = "\\.xlsx$", full.names = FALSE)
for (excel_file in sort(excel_files)) {
  full_path <- file.path(OUTPUT_DIR, excel_file)
  size_kb <- round(file.info(full_path)$size / 1024, 1)
  cat("  ", excel_file, " (", size_kb, " KB)\n", sep = "")
}

cat("\nScript completed successfully!\n")
