from flask import Flask, jsonify, request
from flask_cors import CORS
import subprocess
import os

# Initialize the Flask application
app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Function to run the variant caller using the specified BAM file and reference genome
def run_variant_caller(bam_file, ref_genome):
    # Name of the output VCF file
    vcf_file = "output_variants.vcf"
    
    # Command to run Freebayes for variant calling
    command = f"freebayes -f {ref_genome} {bam_file} > {vcf_file}"
    
    # Execute the command in the shell
    subprocess.run(command, shell=True)
    
    # Return the name of the generated VCF file
    return vcf_file

# Define a route for running the variant caller
@app.route('/run-variant-caller', methods=['POST'])
def run_variant_caller_endpoint():
    # Retrieve BAM file and reference genome from the POST request
    bam_file = request.form['bam_file']
    ref_genome = request.form['ref_genome']
    
    # Construct full paths to the files
    bam_file_path = os.path.join('Dataset', bam_file)
    ref_genome_path = os.path.join('Dataset', ref_genome)
    
    # Run the variant caller and get the output VCF file
    vcf_file = run_variant_caller(bam_file_path, ref_genome_path)
    
    # Return a JSON response indicating success and the VCF file name
    return jsonify({"message": "SNPs identified", "vcf_file": vcf_file})

# Function to classify QUAL values into High, Moderate, or Low
def classify_qual(qual_str):
    # Convert the string to a float for comparison
    qual = float(qual_str)
    
    if qual > 100:
        return "High"
    elif 20 <= qual <= 100:
        return "Moderate"
    else:
        return "Low"

# Define a route to retrieve SNPs from the generated VCF file
@app.route('/get-snps', methods=['GET'])
def get_snps():
    # Get the VCF file name from the request, or use a default name
    vcf_file = request.args.get('vcf_file', 'output_variants.vcf')
    
    # Initialize an empty list to store SNP data
    snps = []

    # Check if the VCF file exists before attempting to open it
    if not os.path.exists(vcf_file):
        return jsonify({"error": "VCF file not found"}), 404

    # Open the VCF file and parse the SNP data
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):  # Skip header lines
                continue
            
            # Split the line into fields
            fields = line.strip().split("\t")
            
            # Extract relevant fields
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            annotations = fields[7]
            
            # Append the SNP data to the list
            snps.append({
                "chrom": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "type": "snp",
                "qual": qual,
                "qual_label": classify_qual(qual),
                "annotations": annotations
            })

    # Return the SNP data as a JSON response
    return jsonify(snps)

# A simple test route to check if the server is running
@app.route('/test')
def test_route():
    return jsonify({"message": "This route is accessible"}), 200

# Main entry point of the application
if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=8000)
