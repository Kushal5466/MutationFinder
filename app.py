from flask import Flask, jsonify, request
from flask_cors import CORS
import subprocess
import os
import matplotlib.pyplot as plt
import io
import base64

# Initialize the Flask application
app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Function to run the variant caller using the specified BAM file and reference genome
def run_variant_caller(bam_file, ref_genome):
    vcf_file = "output_variants.vcf"
    command = f"freebayes -f {ref_genome} {bam_file} > {vcf_file}"
    subprocess.run(command, shell=True)
    return vcf_file

# Define a route for running the variant caller
@app.route('/run-variant-caller', methods=['POST'])
def run_variant_caller_endpoint():
    bam_file = request.form['bam_file']
    ref_genome = request.form['ref_genome']
    bam_file_path = os.path.join('Dataset', bam_file)
    ref_genome_path = os.path.join('Dataset', ref_genome)
    vcf_file = run_variant_caller(bam_file_path, ref_genome_path)
    return jsonify({"message": "SNPs identified", "vcf_file": vcf_file})

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
    print("test")
    vcf_file = request.args.get('vcf_file', 'output_variants.vcf')
    snps = []

    # Check if the VCF file exists before attempting to open it
    if not os.path.exists(vcf_file):
        return jsonify({"error": "VCF file not found"}), 404

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            snps.append({
                "chrom": chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "type": "snp",
                "qual":qual,
                "qual_label": classify_qual(qual)
            })

    snps = snps[:100]
    # variance_sequence = snps[:10]
    # print(len(snps))

    # # Extract positions, QUAL values, and labels
    # positions = [variant["position"] for variant in variance_sequence]
    # qual_values = [float(variant["qual"]) for variant in variance_sequence]
    # qual_labels = [variant["qual_label"] for variant in variance_sequence]

    # # Color mapping based on QUAL label
    # colors = {'Low': 'red', 'Moderate': 'orange', 'High': 'green'}
    # color_values = [colors[label] for label in qual_labels]

    # # Create a line plot (wave-like plot)
    # plt.figure(figsize=(10, 6))
    # plt.plot(positions, qual_values, color='gray', linestyle='--', marker='o', markersize=8)

    # # Color the points based on their QUAL label
    # for pos, qual, color in zip(positions, qual_values, color_values):
    #     plt.scatter(pos, qual, color=color, s=100)  # s is the size of the scatter points

    # plt.title("Wave Plot of QUAL Values Across Positions with Color Coding")
    # plt.xlabel("Position")
    # plt.ylabel("QUAL Value")
    # plt.grid(True)
    # plt.savefig("qual_wave_plot_colored.png")  # Save the plot as a PNG file


    # ## Extract positions and QUAL values and convert QUAL to floats
    # # positions = [variant["position"] for variant in variance_sequence]
    # # qual_values = [float(variant["qual"]) for variant in variance_sequence]

    # # # Create a line plot (wave-like plot)
    # # plt.figure(figsize=(10, 6))
    # # plt.plot(positions, qual_values, marker='o')
    # # plt.title("Wave Plot of QUAL Values Across Positions")
    # # plt.xlabel("Position")
    # # plt.ylabel("QUAL Value")
    # # plt.grid(True)

    # # # Save plot as a PNG file
    # # plt.savefig("qual_wave_plot.png")

    # # Save plot to a bytes buffer for base64 encoding
    # buf = io.BytesIO()
    # plt.savefig(buf, format='png')
    # buf.seek(0)

    # # Encode the bytes to base64
    # plot_base64 = base64.b64encode(buf.read()).decode('utf-8')
    # buf.close()

    # # # Extract QUAL values and convert them to floats
    # # qual_values = [float(variant["qual"]) for variant in variance_sequence]

    # # # Create a plot
    # # plt.figure(figsize=(8, 6))
    # # plt.boxplot(qual_values, vert=False, patch_artist=True)
    # # plt.title("Box Plot of QUAL Values")
    # # plt.xlabel("QUAL Value")

    # # # Save plot to a bytes buffer
    # # buf = io.BytesIO()
    # # plt.savefig(buf, format='png')
    # # buf.seek(0)

    # # # Encode the bytes to base64
    # # plot_base64 = base64.b64encode(buf.read()).decode('utf-8')
    # # buf.close()

    # # Return the base64 string as a JSON response
    # return jsonify({'plot': plot_base64})
    return jsonify(snps)

@app.route('/test')
def test_route():
    return jsonify({"message": "This route is accessible"}), 200


if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=8000)

