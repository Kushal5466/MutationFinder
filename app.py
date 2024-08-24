from flask import Flask, jsonify, request
import subprocess
import os

# Initialize the Flask application
app = Flask(__name__)

# Function to run the variant caller using the specified BAM file and reference genome
def run_variant_caller(bam_file, ref_genome):
    # Name of the output VCF file where the variant caller will store the results
    vcf_file = "output_variants.vcf"
    
    # Command to run the variant caller (e.g., FreeBayes) with the reference genome and BAM file as inputs
    # The results will be written to the VCF file
    command = f"freebayes -f {ref_genome} {bam_file} > {vcf_file}"
    
    # Execute the command in the shell
    subprocess.run(command, shell=True)
    
    # Return the name of the generated VCF file
    return vcf_file

# Define a route for running the variant caller
@app.route('/run-variant-caller', methods=['POST'])
def run_variant_caller_endpoint():
    # Extract the BAM file and reference genome file names from the POST request's form data
    bam_file = request.form['bam_file']
    ref_genome = request.form['ref_genome']
    
    # Construct the full file paths by joining the file names with the 'Dataset' directory
    # This ensures the files are correctly located within the server's file system
    bam_file_path = os.path.join('Dataset', bam_file)
    ref_genome_path = os.path.join('Dataset', ref_genome)
    
    # Call the function to run the variant caller with the specified file paths
    vcf_file = run_variant_caller(bam_file_path, ref_genome_path)
    
    # Return a JSON response with a success message and the name of the generated VCF file
    return jsonify({"message": "SNPs identified", "vcf_file": vcf_file})

# Define a route to retrieve SNPs from the generated VCF file
@app.route('/get-snps', methods=['GET'])
def get_snps():
    # Get the VCF file name from the query parameters, default to 'output_variants.vcf' if not provided
    vcf_file = request.args.get('vcf_file', 'output_variants.vcf')
    snps = []

    # Open the VCF file and read it line by line
    with open(vcf_file, 'r') as file:
        for line in file:
            # Skip header lines that start with '#'
            if line.startswith("#"):
                continue
            
            # Split the line into fields based on tab separation
            fields = line.strip().split("\t")
            chrom = fields[0]       # Chromosome
            pos = int(fields[1])     # Position on the chromosome
            ref = fields[3]          # Reference allele
            alt = fields[4]          # Alternate allele
            
            # Append the SNP information to the list as a dictionary
            snps.append({
                "chrom": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "type": "snp"         # Specify the type of variant as SNP
            })

    # Return the list of SNPs as a JSON response
    return jsonify(snps)

# Run the Flask application in debug mode when the script is executed directly
if __name__ == "__main__":
    app.run(debug=True)
