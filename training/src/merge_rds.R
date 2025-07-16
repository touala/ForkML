#!/usr/bin/env Rscript

library(docopt)
"Usage:
    merge_rds.R -o=FILE FILE...

Mandatory arguments:
    -o FILE, --output FILE      Output file
    
" -> doc
args = docopt(doc)
file_out=args$output
files_in=args$FILE


# Initialize all_data
all_data=NULL

# For all files
for (file in files_in) {
	# Read data
    data = readRDS(file)
    # Merge tables
    if (is.data.frame(data)) {
    	all_data = rbind(all_data, data)
    }
    # Merge lists
    else if (is.list(data)) {
    	all_data = c(all_data, data)
    }
}

# Save data
saveRDS(all_data, file_out)
