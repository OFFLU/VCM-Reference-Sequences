#!/bin/bash

# COMPLETE PHYLOGENETIC ANALYSIS PIPELINE
# Runs RAxML → Splicer → PARNAS → Consolidation (ignoring CVV files)
# Processes regular .fasta files only, skips *_with_cvvs.fasta

echo "======================================================================="
echo "COMPLETE PHYLOGENETIC ANALYSIS PIPELINE"
echo "======================================================================="
echo "Pipeline: RAxML → Splicer → PARNAS (10-15 reps) → Consolidation"
echo "Ignoring: *_with_cvvs.fasta files"
echo ""

# Automatically discover folders containing FASTA files
# This will process any folder that contains a .fasta file (excluding CVV files)

echo "Discovering folders with FASTA files..."

# Function to find folders with appropriate FASTA files
discover_folders() {
    local found_folders=()
    
    # Look for folders containing .fasta files (but not *_with_cvvs.fasta)
    for dir in */; do
        # Remove trailing slash
        folder="${dir%/}"
        
        # Skip if not a directory
        [ ! -d "$folder" ] && continue
        
        # Look for FASTA files in this folder
        # Check for folder-named FASTA file first (e.g., H3_Am.fasta in H3_Am/)
        if [ -f "$folder/${folder}.fasta" ]; then
            echo "  Found: $folder (contains ${folder}.fasta)"
            found_folders+=("$folder")
        # Check for reference.fasta in reference folder
        elif [ "$folder" = "reference" ] && [ -f "$folder/reference.fasta" ]; then
            echo "  Found: $folder (contains reference.fasta)"
            found_folders+=("$folder")
        # Check for any .fasta file that's not a CVV file
        else
            fasta_files=($(find "$folder" -maxdepth 1 -name "*.fasta" ! -name "*_with_cvvs.fasta" 2>/dev/null))
            if [ ${#fasta_files[@]} -gt 0 ]; then
                echo "  Found: $folder (contains ${#fasta_files[@]} FASTA files)"
                found_folders+=("$folder")
            fi
        fi
    done
    
    # Return the array
    printf '%s\n' "${found_folders[@]}"
}

# Discover folders automatically
echo ""
# Create temporary file for folder list (more compatible)
temp_folder_list=$(mktemp)
discover_folders > "$temp_folder_list"

# Read folders into array (compatible with older bash)
folders=()
while IFS= read -r folder; do
    [ -n "$folder" ] && folders+=("$folder")
done < "$temp_folder_list"

# Clean up
rm -f "$temp_folder_list"

if [ ${#folders[@]} -eq 0 ]; then
    echo "Error: No folders with FASTA files found!"
    echo "Expected folder structure:"
    echo "  folder_name/"
    echo "  ├── folder_name.fasta  (or reference.fasta for reference folder)"
    echo "  └── ..."
    echo ""
    echo "Make sure you're running this script in the directory containing your data folders."
    exit 1
fi

echo ""
echo "Will process ${#folders[@]} folders:"
for folder in "${folders[@]}"; do
    echo "  - $folder"
done

# PARNAS configuration
TARGET_REPRESENTATIVES=(10 12 15)
decomp_folder="mydecomposition"  # Only process regular decomposition

# Record overall start time
overall_start=$(date)
echo "Pipeline started: $overall_start"
echo ""

# =============================================================================
# PHASE 1: RAxML and Splicer Analysis
# =============================================================================

echo "PHASE 1: RAxML AND SPLICER ANALYSIS"
echo "======================================================================="
echo ""

# Function to process each folder for RAxML and Splicer
process_raxml_splicer() {
    local folder=$1
    echo "Processing folder: $folder"
    
    # Check if folder exists
    if [ ! -d "$folder" ]; then
        echo "  Warning: Folder $folder does not exist, skipping..."
        return 1
    fi
    
    # Change to the folder
    cd "$folder" || {
        echo "  Error: Could not change to directory $folder"
        return 1
    }
    
    # Look for the main FASTA file (folder-named file)
    fasta_file="${folder}.fasta"
    
    # If not found, try common alternatives
    if [ ! -f "$fasta_file" ]; then
        # Try reference.fasta for reference folder
        if [ "$folder" = "reference" ] && [ -f "reference.fasta" ]; then
            fasta_file="reference.fasta"
        else
            echo "  Warning: No suitable FASTA file found in $folder"
            echo "  Expected: ${folder}.fasta or reference.fasta"
            cd ..
            return 1
        fi
    fi
    
    echo "  Using FASTA file: $fasta_file"
    
    # Step 1: RAxML-NG
    echo "  Step 1: Running RAxML-NG..."
    if raxml-ng --fast --msa "$fasta_file" --model DNA; then
        echo "  ✓ RAxML-NG completed successfully"
        
        # Check if the tree file was created
        tree_file="${fasta_file}.raxml.bestTree"
        if [ -f "$tree_file" ]; then
            # Step 2: Splicer decomposition
            echo "  Step 2: Running Splicer decomposition..."
            if splicer decomp -t "$tree_file" -s "$fasta_file" -n "$decomp_folder"; then
                echo "  ✓ Splicer decomposition completed successfully"
                echo "  ✓ Created: $decomp_folder/"
            else
                echo "  ✗ Error: Splicer decomposition failed"
                cd ..
                return 1
            fi
        else
            echo "  ✗ Error: Tree file not found, skipping Splicer step"
            cd ..
            return 1
        fi
    else
        echo "  ✗ Error: RAxML-NG failed"
        cd ..
        return 1
    fi
    
    # Return to parent directory
    cd ..
    echo "  ✓ Finished processing $folder"
    echo ""
    
    return 0
}

# Process each folder for RAxML and Splicer
raxml_splicer_success=0
for folder in "${folders[@]}"; do
    if process_raxml_splicer "$folder"; then
        raxml_splicer_success=$((raxml_splicer_success + 1))
    fi
done

echo "Phase 1 Summary: $raxml_splicer_success/${#folders[@]} folders processed successfully"
echo ""

# =============================================================================
# PHASE 2: PARNAS Analysis (10-15 representatives)
# =============================================================================

echo "PHASE 2: PARNAS ANALYSIS"
echo "======================================================================="
echo "Target representative counts: ${TARGET_REPRESENTATIVES[*]}"
echo "Selection strategy: Diversity-based optimal selection"
echo ""

# Function to run PARNAS analysis on a scaffold tree with multiple representative counts
run_parnas_analysis() {
    local folder=$1
    local scaffold_tree="$folder/$decomp_folder/scaffold_tree.plain.tre"
    
    echo "Processing PARNAS for: $folder"
    echo "  Scaffold tree: $scaffold_tree"
    
    # Check if scaffold tree exists
    if [ ! -f "$scaffold_tree" ]; then
        echo "  ✗ Warning: Scaffold tree not found, skipping..."
        return 1
    fi
    
    # Create output directory for PARNAS results
    output_dir="$folder/$decomp_folder/parnas_output"
    mkdir -p "$output_dir"
    
    echo "  → Running PARNAS diversity analysis..."
    
    # First, run diversity analysis to understand the dataset
    diversity_file="$output_dir/diversity_analysis.csv"
    if ! ~/.local/bin/parnas -t "$scaffold_tree" -n 20 --diversity "$diversity_file" 2>/dev/null; then
        echo "  ✗ Error: PARNAS diversity analysis failed"
        return 1
    fi
    
    # Parse diversity file to determine optimal representative count
    if [ -f "$diversity_file" ]; then
        # Find diversity achieved by different representative counts
        diversity_10=$(awk -F',' '$1 == 10 { print $2 }' "$diversity_file" 2>/dev/null || echo "0")
        diversity_12=$(awk -F',' '$1 == 12 { print $2 }' "$diversity_file" 2>/dev/null || echo "0")
        diversity_15=$(awk -F',' '$1 == 15 { print $2 }' "$diversity_file" 2>/dev/null || echo "0")
        
        # Decision logic: if 10 reps achieve >75% diversity, use 10; 
        # if 15 reps achieve >85%, use 15; otherwise use 12
        if (( $(echo "$diversity_10 > 75" | bc -l 2>/dev/null || echo "0") )); then
            optimal_n=10
            rationale="10 representatives achieve ${diversity_10}% diversity (>75% threshold)"
        elif (( $(echo "$diversity_15 > 85" | bc -l 2>/dev/null || echo "0") )); then
            optimal_n=15
            rationale="15 representatives needed for good coverage (${diversity_15}% diversity)"
        else
            optimal_n=12
            rationale="12 representatives as balanced choice (${diversity_12}% diversity)"
        fi
        
        echo "  → Dataset assessment: $rationale"
        echo "  → Selected optimal: $optimal_n representatives"
    else
        echo "  → Warning: Could not assess diversity, defaulting to 12 representatives"
        optimal_n=12
        rationale="Default selection (diversity analysis failed)"
    fi
    
    # Create summary file for this folder
    summary_file="$output_dir/representative_selection_summary.txt"
    echo "PARNAS Representative Selection Summary for $folder" > "$summary_file"
    echo "Generated on: $(date)" >> "$summary_file"
    echo "======================================================" >> "$summary_file"
    echo "Selected representatives: $optimal_n" >> "$summary_file"
    echo "Rationale: $rationale" >> "$summary_file"
    echo "" >> "$summary_file"
    
    # Test all target numbers and provide comparison
    echo "Comparison of different representative counts:" >> "$summary_file"
    printf "%-6s %-12s %-25s %-25s %-30s\n" "Count" "Diversity%" "Colored Tree" "Clusters" "Subtree" >> "$summary_file"
    echo "---------------------------------------------------------------------------------" >> "$summary_file"
    
    success_count=0
    for n in "${TARGET_REPRESENTATIVES[@]}"; do
        echo "  → Testing $n representatives..."
        
        # Set output files for this count
        colored_tree="$output_dir/parnas_n${n}_representatives.tre"
        clusters_file="$output_dir/parnas_n${n}_clusters.tab"
        subtree_file="$output_dir/parnas_n${n}_subtree.tre"
        
        # Run PARNAS with this number
        if ~/.local/bin/parnas -t "$scaffold_tree" \
            -n "$n" \
            --color "$colored_tree" \
            --clusters "$clusters_file" \
            --subtree "$subtree_file" 2>/dev/null; then
            
            # Get diversity for this n from our analysis
            diversity_n=$(awk -F',' -v target="$n" '$1 == target { print $2 }' "$diversity_file" 2>/dev/null || echo "unknown")
            
            echo "    ✓ Success - ${diversity_n}% diversity coverage"
            
            # Add to summary
            marker=""
            if [ "$n" = "$optimal_n" ]; then
                marker=" ← SELECTED"
            fi
            
            printf "%-6s %-12s %-25s %-25s %-30s%s\n" \
                "$n" \
                "${diversity_n}%" \
                "parnas_n${n}_representatives.tre" \
                "parnas_n${n}_clusters.tab" \
                "parnas_n${n}_subtree.tre" \
                "$marker" >> "$summary_file"
            
            success_count=$((success_count + 1))
        else
            echo "    ✗ Failed"
            printf "%-6s %-12s %-25s %-25s %-30s\n" \
                "$n" \
                "FAILED" \
                "N/A" \
                "N/A" \
                "N/A" >> "$summary_file"
        fi
    done
    
    # Add usage recommendations to summary
    echo "" >> "$summary_file"
    echo "USAGE RECOMMENDATIONS:" >> "$summary_file"
    echo "- Use the SELECTED files (n=$optimal_n) for downstream analyses" >> "$summary_file"
    echo "- parnas_n${optimal_n}_subtree.tre contains exactly $optimal_n representative sequences" >> "$summary_file"
    echo "- parnas_n${optimal_n}_clusters.tab shows which sequences each representative covers" >> "$summary_file"
    echo "- Open parnas_n${optimal_n}_representatives.tre in FigTree to visualize coverage" >> "$summary_file"
    
    echo "  ✓ Representative selection completed: $success_count/${#TARGET_REPRESENTATIVES[@]} counts successful"
    echo "  ✓ Recommended subtree: parnas_n${optimal_n}_subtree.tre ($optimal_n representatives)"
    echo ""
    
    return 0
}

# Process each folder for PARNAS
parnas_success=0
for folder in "${folders[@]}"; do
    # Check if decomposition exists
    if [ -d "$folder/$decomp_folder" ]; then
        if run_parnas_analysis "$folder"; then
            parnas_success=$((parnas_success + 1))
        fi
    else
        echo "Processing PARNAS for: $folder"
        echo "  ✗ No decomposition folder found: $folder/$decomp_folder"
        echo ""
    fi
done

echo "Phase 2 Summary: $parnas_success folders processed successfully for PARNAS"
echo ""

# =============================================================================
# PHASE 3: Consolidate PARNAS Outputs
# =============================================================================

echo "PHASE 3: CONSOLIDATING PARNAS OUTPUTS"
echo "======================================================================="
echo ""

# Create consolidated output directory
consolidated_dir="consolidated_parnas_outputs"
mkdir -p "$consolidated_dir"

echo "Created consolidated directory: $consolidated_dir"
echo ""

# Function to select optimal files for a group
consolidate_group_outputs() {
    local folder=$1
    local parnas_dir="$folder/$decomp_folder/parnas_output"
    
    echo "Consolidating: $folder"
    
    # Check if PARNAS output directory exists
    if [ ! -d "$parnas_dir" ]; then
        echo "  ✗ No PARNAS output directory found"
        return 1
    fi
    
    # Check if diversity analysis file exists
    diversity_file="$parnas_dir/diversity_analysis.csv"
    if [ ! -f "$diversity_file" ]; then
        echo "  ✗ No diversity analysis file found"
        return 1
    fi
    
    # Find the optimal count (same logic as PARNAS selection)
    diversity_10=$(awk -F',' '$1 == 10 { print $2 }' "$diversity_file" 2>/dev/null || echo "0")
    diversity_12=$(awk -F',' '$1 == 12 { print $2 }' "$diversity_file" 2>/dev/null || echo "0")
    diversity_15=$(awk -F',' '$1 == 15 { print $2 }' "$diversity_file" 2>/dev/null || echo "0")
    
    if (( $(echo "$diversity_10 > 75" | bc -l 2>/dev/null || echo "0") )); then
        optimal_n=10
        optimal_diversity="$diversity_10"
    elif (( $(echo "$diversity_15 > 85" | bc -l 2>/dev/null || echo "0") )); then
        optimal_n=15
        optimal_diversity="$diversity_15"
    else
        optimal_n=12
        optimal_diversity="$diversity_12"
    fi
    
    echo "  → Selected: $optimal_n representatives (${optimal_diversity}% diversity)"
    
    # Define source files
    source_colored="$parnas_dir/parnas_n${optimal_n}_representatives.tre"
    source_clusters="$parnas_dir/parnas_n${optimal_n}_clusters.tab"
    source_subtree="$parnas_dir/parnas_n${optimal_n}_subtree.tre"
    
    # Define destination files
    dest_colored="$consolidated_dir/${folder}_optimal_colored.tre"
    dest_clusters="$consolidated_dir/${folder}_optimal_clusters.tab"
    dest_subtree="$consolidated_dir/${folder}_optimal_subtree.tre"
    
    # Copy files with error checking
    files_copied=0
    
    if [ -f "$source_colored" ]; then
        cp "$source_colored" "$dest_colored"
        echo "  ✓ Copied colored tree: ${folder}_optimal_colored.tre"
        files_copied=$((files_copied + 1))
    else
        echo "  ✗ Colored tree not found: $source_colored"
    fi
    
    if [ -f "$source_clusters" ]; then
        cp "$source_clusters" "$dest_clusters"
        echo "  ✓ Copied clusters: ${folder}_optimal_clusters.tab"
        files_copied=$((files_copied + 1))
    else
        echo "  ✗ Clusters file not found: $source_clusters"
    fi
    
    if [ -f "$source_subtree" ]; then
        cp "$source_subtree" "$dest_subtree"
        echo "  ✓ Copied subtree: ${folder}_optimal_subtree.tre"
        files_copied=$((files_copied + 1))
    else
        echo "  ✗ Subtree not found: $source_subtree"
    fi
    
    echo "  → Files copied: $files_copied/3"
    echo ""
    
    return 0
}

# Consolidate outputs for each processed folder
consolidation_success=0
for folder in "${folders[@]}"; do
    if [ -d "$folder/$decomp_folder/parnas_output" ]; then
        if consolidate_group_outputs "$folder"; then
            consolidation_success=$((consolidation_success + 1))
        fi
    fi
done

echo "Phase 3 Summary: $consolidation_success folders consolidated successfully"
echo ""

# Create summary file for consolidated outputs
summary_file="$consolidated_dir/consolidation_summary.txt"
echo "PARNAS Consolidation Summary" > "$summary_file"
echo "Generated on: $(date)" >> "$summary_file"
echo "==============================" >> "$summary_file"
echo "" >> "$summary_file"
echo "Consolidated files from $consolidation_success folders:" >> "$summary_file"
echo "" >> "$summary_file"

# List consolidated files
consolidated_files=($(ls "$consolidated_dir"/*.tre "$consolidated_dir"/*.tab 2>/dev/null || true))
if [ ${#consolidated_files[@]} -gt 0 ]; then
    for file in "${consolidated_files[@]}"; do
        filename=$(basename "$file")
        size_kb=$(du -k "$file" | cut -f1)
        echo "$filename (${size_kb} KB)" >> "$summary_file"
    done
fi

echo "" >> "$summary_file"
echo "File types:" >> "$summary_file"
echo "- *_optimal_colored.tre    - Colored trees showing representative coverage" >> "$summary_file"
echo "- *_optimal_clusters.tab   - Cluster assignments for representatives" >> "$summary_file"
echo "- *_optimal_subtree.tre    - Subtrees with selected representatives only" >> "$summary_file"

# =============================================================================
# FINAL SUMMARY
# =============================================================================

# Record overall end time
overall_end=$(date)

echo "======================================================================="
echo "COMPLETE PIPELINE FINISHED SUCCESSFULLY!"
echo "======================================================================="
echo "Started:  $overall_start"
echo "Finished: $overall_end"
echo ""
echo "PHASE SUMMARY:"
echo "- Phase 1 (RAxML/Splicer): $raxml_splicer_success/${#folders[@]} folders processed"
echo "- Phase 2 (PARNAS): $parnas_success folders processed"
echo "- Phase 3 (Consolidation): $consolidation_success folders consolidated"
echo ""
echo "OUTPUT STRUCTURE:"
echo "1. Individual folder results:"
echo "   [folder]/mydecomposition/parnas_output/ - Full PARNAS analysis"
echo ""
echo "2. Consolidated results:"
echo "   $consolidated_dir/ - Optimal files from each group"
echo ""
echo "CONSOLIDATED FILES:"
if [ -d "$consolidated_dir" ]; then
    consolidated_count=$(ls "$consolidated_dir"/*.tre 2>/dev/null | wc -l)
    echo "- $consolidated_count optimal subtree files ready for analysis"
    echo "- Each represents 10-15 sequences covering optimal diversity"
    echo "- Use consolidated subtrees for downstream phylogenetic analysis"
    echo ""
    echo "Example files in $consolidated_dir/:"
    ls "$consolidated_dir" | head -6 | sed 's/^/  /'
    total_files=$(ls "$consolidated_dir" 2>/dev/null | wc -l)
    if [ "$total_files" -gt 6 ]; then
        echo "  ... and $((total_files - 6)) more files"
    fi
else
    echo "- No consolidated files created"
fi

echo ""
echo "NEXT STEPS:"
echo "1. Review consolidated subtrees in $consolidated_dir/"
echo "2. Use *_optimal_subtree.tre files for downstream analysis"
echo "3. Open *_optimal_colored.tre in FigTree to visualize representatives"
echo "4. Check diversity scores in individual parnas_output directories"
echo ""
echo "Pipeline completed successfully!"
