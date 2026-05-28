#!/usr/bin/env python3

# RM html to tsv
# by Diego De Panis, 2025
# This script is part of the GAME pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

"""
Convert RepeatMasker HTML landscape output to TSV format for R plotting.
"""

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Convert RepeatMasker HTML landscape to TSV for R plotting"
    )
    parser.add_argument(
        "-html", "--html", 
        required=True, 
        help="Input HTML file from createRepeatLandscape.pl"
    )
    parser.add_argument(
        "-out", "--out", 
        required=True, 
        help="Output TSV file"
    )
    
    args = parser.parse_args()
    
    try:
        with open(args.html, 'r') as infile, open(args.out, 'w') as outfile:
            first = True
            data = False
            
            for line in infile:
                line = line.rstrip("\n").lstrip(" ")
                
                # Parse column headers
                if line.startswith("data.addColumn"):
                    val = line.rstrip("');\n").split(",")[1]
                    if first:
                        outfile.write(val.lstrip(" '"))
                        first = False
                    else:
                        outfile.write("\t" + val.lstrip(" '"))
                
                # Start of data section
                elif line.startswith("data.addRows"):
                    data = True
                    outfile.write("\n")
                
                # End of data section
                elif line.startswith("]);"):
                    data = False
                
                # Parse data rows
                elif line.startswith("[") and data:
                    values = line.lstrip("['").rstrip("],\n").split(",")
                    values[0] = values[0].rstrip("'")
                    tab_line = "\t".join(values)
                    outfile.write(tab_line + "\n")
        
        print(f"Successfully converted {args.html} to {args.out}")
        
    except FileNotFoundError:
        print(f"Error: Could not find input file {args.html}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
