# Create empty DataFrames for storing matching rows
matching_rows_df = pd.DataFrame(columns=['Product', 'Count', 'Target Gene', 'Target Locus'])
matching_rows_reverse_df = pd.DataFrame(columns=['Product', 'Count', 'Target Gene', 'Target Locus'])

# Iterate over each row in result_df
for _, row in result_df.iterrows():
    key = (row['Product'], row['Length'])
    if key in primer_dict:
        matching_row = {
            'Product': row['Product'],
            'Count': row['Count'],
            'Target Gene': primer_dict[key]['target_gene'],
            'Target Locus': primer_dict[key]['target_locus']
        }
        matching_rows_df = matching_rows_df.append(matching_row, ignore_index=True)

# Iterate over each row in result_reverse_df
for _, row in result_reverse_df.iterrows():
    key = (row['Product'], row['Length'])
    if key in reverse_primer_dict:
        matching_row = {
            'Product': row['Product'],
            'Count': row['Count'],
            'Target Gene': reverse_primer_dict[key]['target_gene'],
            'Target Locus': reverse_primer_dict[key]['target_locus']
        }
        matching_rows_reverse_df = matching_rows_reverse_df.append(matching_row, ignore_index=True)

# Now you can output these DataFrames to CSV files or use them for further processing
matching_rows_df.to_csv('matching_rows.csv', index=False)
matching_rows_reverse_df.to_csv('matching_rows_reverse.csv', index=False)
