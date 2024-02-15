
# Define the input TSV file
input_file=$1

# Define the output TSV file
output_file=$2

# Read the selected row indices from the text file
selected_row_indices_file=$3

# Read row indices from row_indices.txt into an array
mapfile -t indices < $selected_row_indices_file

# Create an awk script dynamically based on the indices
awk_script='BEGIN { FS = "\t"; OFS = "\t" } '

for index in "${indices[@]}"; do
    awk_script="${awk_script}NR == $index { print; } "
done

# Use awk to select rows based on the dynamically generated script
awk "$awk_script" $input_file > $output_file