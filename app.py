import subprocess
import os
from flask import Flask, jsonify, request

# Initialize the Flask application
app = Flask(__name__)

# Function to run the variant caller using the specified BAM file and reference genome
def run_variant_caller(bam_file, ref_genome):
    # Name of the output VCF file where the variant caller will store the results
    vcf_file = "output_variants.vcf"
    
    # Command to run the variant caller (e.g., FreeBayes) with the reference genome and BAM file as inputs
    # The results will be written to the VCF file
    command = f"freebayes -f {ref_genome} {bam_file} > {vcf_file}"

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

# Run the Flask application in debug mode when the script is executed directly
if __name__ == "__main__":
    app.run(debug=True)
