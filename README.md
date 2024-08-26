# MutationFinder

web application to identify and visualize single nucleotide polymorphisms (SNPs) in a sample whole genome sequencing dataset.

Before Starting of the Project Let's dicuss some of the tools that I have used to complete this project. I have list steps for mac

#Install Homebrew if you don't have it.
1./bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install BWA. BWA is used for aligning short reads to the reference genome:

2.brew install bwa

# Samtools is required for manipulating and processing the SAM/BAM files:

3.brew install samtools

# Download the Reference Genome

# Navigate to your working directory

cd /path/to/your/working/directory

# Download the reference genome

curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/571/405/GCA_002571405.2_ASM257140v2/GCA_002571405.2_ASM257140v2_genomic.fna.gz

# Download the Short Read Data

# You can use the SRA toolkit to download the short read data, but an easier method is to use the fastq-dump command from the SRA toolkit:

# Install the SRA toolkit

brew install sratoolkit

# Download the short read data

fastq-dump --split-files --gzip SRR5892450

# Align the short reads to the reference genome using BWA:

# Index the reference genome

bwa index GCA_002571405.2_ASM257140v2_genomic.fna

# Align the short reads to the reference genome

bwa mem GCA_002571405.2_ASM257140v2_genomic.fna SRR5892450_1.fastq.gz SRR5892450_2.fastq.gz > alignment.sam

# Convert the SAM file to BAM format and sort it using Samtools:

# Convert SAM to BAM

samtools view -S -b alignment.sam > alignment.bam

# Sort the BAM file

samtools sort alignment.bam -o alignment_sorted.bam

# Index the BAM file

samtools index alignment_sorted.bam

Mutation-finder/
│
├── backend/
│ ├── app.py # Flask application entry point contains routes.
│ ├── dataset/
│ │ ├──GCA_002571405.2_ASM257140v2_genomic.fna # Reference genome file
│ │ ├── alignment_reads_sorted.bam # Aligned BAM file (output from alignment step)
│ │ ├── output_variants.vcf # VCF file containing identified SNPs (output from variant caller)
│ │
│ ├── requirements.txt # Python dependencies for the Flask backend
│ └── README.md # Instructions for setting up and running the backend

To run the Backend Server Follow these Steps.

1.First Move to into the MutationFinder Directory. 2. Install all the dependencies from requirements.txt 3. To start the application you can use the command python3 app.py 4. Check whether you have all dependencies are installed. 5. Check the API responses in Postman or Swagger to ensure they are returning the expected data.

frontend
│
├── snp-visualization/
│ ├── react-app/ # React application directory (if using React)
│ │ ├── src/
│ │ │ ├── App.js # Main React component
| | | |-- components
│ │ │ | ├── GeonmeViewer.js # Component to visualize SNPs
│ │ │ | ├── ReferenceSequenceViewer.js # Component to display the reference sequence
│ │ │ ├── SNPSummaryTable.js # Component to show SNPs in a table format
│ │ | |-- VarianceBoxPlot.js #Displays the QUAL value and position
│ │ | |-- ZoomControls.js # Zoom in and Zoom out options.
│ │ ├── public/
│ │ ├── package.json # React dependencies
│ │ └── README.md # Instructions for setting up and running the
