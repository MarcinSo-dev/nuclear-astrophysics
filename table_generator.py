import numpy as np

def csv_to_latex_table(csv_file, output_file, delimiter=',', decimal_mark='.'): 
    # Load the CSV file into a numpy array
    with open(csv_file, 'r') as f:
        data = np.genfromtxt(f, delimiter=delimiter, dtype=str)

    # Replace the decimal mark if needed
    if decimal_mark != '.':
        data = np.char.replace(data, decimal_mark, '.')

    # Separate headers and data
    headers = data[0]  # First row as headers
    rows = data[1:]    # Remaining rows as data

    # Open the output file to write the LaTeX code
    with open(output_file, 'w') as f:
        # Start the LaTeX code
        f.write("\\begin{longtable}{" + "".join(["c"] * len(headers)) + "}\hline\n")
        
        # Write the header row
        header = " & ".join([f"\\textbf{{{col}}}" for col in headers]) + " \\\\ \\\hline\n"
        f.write(header)
        f.write("\\endfirsthead\n")
        f.write(header)  # Repeat header for subsequent pages
        f.write("\\endhead\n\\hline\n\\endfoot\n\\hline\\endlastfoot\n")

        # Write the data rows with alternating colors
        for i, row in enumerate(rows):
            row_data = " & ".join(row) + " \\\\ \n"
            if i % 2 == 0:
                f.write("\\rowcolor[HTML]{C0C0C0} ")  # Gray background for even rows
            f.write(row_data)
            f.write("\\hline\n")

        # End the longtable environment
        f.write("\\end{longtable}")

# Example usage
csv_file = "example.csv"  # Path to your CSV file
output_file = "table.tex"  # Output LaTeX file
csv_to_latex_table(csv_file, output_file, delimiter=';', decimal_mark=',')

print(f"LaTeX table code has been written to {output_file}")
