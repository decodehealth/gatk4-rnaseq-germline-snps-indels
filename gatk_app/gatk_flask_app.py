from flask import Flask, request
import subprocess

app = Flask(__name__)

# Create variables for fastq_r1_path, fastq_r2_path, genome_dir_path, num_threads,
# out_SAM_type, twopass_mode, outfile_name_prefix, genome_SA_sparse_D,
# genome_Chr_Bin_N_bits, and log_path. Create log file to write to, created a string
# for invoking STAR, invoked command for STAR, and ran the command in a shell, logging
# stdout and stderr and returning a pid number
  
@app.route("/run_gatk_mark_dupes", methods=["POST"])
def run_gatk_mark_dupes():
    rjson = request.get_json()
    ## get vals
    path_to_input_bam = rjson["input_bam_path"]
    path_to_deduped_bam = rjson["path_to_deduped_bam"]
    str_validation_stringency = rjson["validation_stringency_str"]
    path_to_metrics_file = rjson["metrics_file_path"]
    path_to_log = rjson["log_path"]

    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar MarkDuplicates "
      f"--INPUT {path_to_input_bam}	--OUTPUT {path_to_deduped_bam} --CREATE_INDEX true "
      f"--VALIDATION_STRINGENCY {str_validation_stringency} --METRICS_FILE {path_to_metrics_file}"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)

@app.route("/run_gatk_splitncigar", methods=["POST"])
def run_gatk_splitncigar():
    rjson = request.get_json()
    ## get vals
    path_to_dict_fa = rjson["path_to_dict_fa"]
    path_to_deduped_bam = rjson["path_to_deduped_bam"]
    path_to_split_bam = rjson["path_to_split_bam"]
    path_to_log = rjson["log_path"]
    
    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar SplitNCigarReads -R {path_to_dict_fa} -I {path_to_deduped_bam} -O {path_to_split_bam}"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)

@app.route("/run_gatk_labelread", methods=["POST"])
def run_gatk_labelread():
    rjson = request.get_json()
    ## get vals
    path_to_labeled_bam = rjson["path_to_labeled_bam"]
    path_to_split_bam = rjson["path_to_split_bam"]
    sample_name = rjson["sample_name"]
    path_to_log = rjson["log_path"]
    
    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/picard.jar AddOrReplaceReadGroups I= {path_to_split_bam} O= {path_to_labeled_bam} RGPL=illumina RGPU=unit1 RGLB=lib1 RGSM={sample_name} RGID=1	"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)

@app.route("/run_gatk_baserecal", methods=["POST"])
def run_gatk_baserecal():
    rjson = request.get_json()
    ## get vals
    path_to_dict_fa = rjson["path_to_dict_fa"]
    path_to_bam = rjson["path_to_bam"]
    path_to_csv = rjson["path_to_csv"]
    path_to_log = rjson["log_path"]
    
    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar BaseRecalibrator -R {path_to_dict_fa} -I {path_to_bam} --use-original-qualities -O {path_to_csv} -known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_clinically_associated.vcf -known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_phenotype_associated.vcf -known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_somatic_incl_consequences.vcf -known-sites /mnt/iquityazurefileshare1/GTFs/SNPs/homo_sapiens_structural_variations.vcf"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)

@app.route("/run_gatk_bqsr", methods=["POST"])
def run_gatk_bqsr():
    rjson = request.get_json()
    ## get vals
    path_to_dict_fa = rjson["path_to_dict_fa"]
    path_to_labeled_bam = rjson["path_to_labeled_bam"]
    path_to_recaled_bam = rjson["path_to_recaled_bam"]
    path_to_recaled_csv = rjson["path_to_recaled_csv"]
    path_to_log = rjson["log_path"]
      
    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar ApplyBQSR --add-output-sam-program-record  -R {path_to_dict_fa} -I {path_to_labeled_bam} --use-original-qualities -O {path_to_recaled_bam} --bqsr-recal-file {path_to_recaled_csv}"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)

@app.route("/run_gatk_analyzecov", methods=["POST"])
def run_gatk_analyzecov():
    rjson = request.get_json()
    ## get vals
    path_to_labeled_csv = rjson["path_to_labeled_csv"]
    path_to_recaled_csv = rjson["path_to_recaled_csv"]
    path_to_cov_pdf = rjson["path_to_cov_pdf"]
    path_to_log = rjson["log_path"]
    
    log = open(path_to_log, "w")
    
    cmd = ("conda activate gatk")
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    log.flush()
 
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar AnalyzeCovariates -before {path_to_labeled_csv} -after {path_to_recaled_csv} -plots {path_to_cov_pdf}"
    )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    log.flush()
    return str(proc.pid + 1)  

  
@app.route("/run_gatk_halotypecaller", methods=["POST"])
def run_gatk_halotypecaller():
    rjson = request.get_json()
    ## get vals
    path_to_vcf = rjson["path_to_vcf"]
    path_to_dict_fa = rjson["path_to_dict_fa"]
    path_to_recaled_bam = rjson["path_to_recaled_bam"]
    path_to_log = rjson["log_path"]
    
    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar HaplotypeCaller -R {path_to_dict_fa} -I {path_to_recaled_bam} -L /mnt/iquityazurefileshare1/GTFs/Homo_sapiens.GRCh38.104.gtf.exons.interval_list -O {path_to_vcf} -dont-use-soft-clipped-bases false --standard-min-confidence-threshold-for-calling 30"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)

## optional --dbsnp flag
##	--dbsnp /mnt/iquityazurefileshare1/IQuity/GATK_Testing/gatk-workflows/inputs/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf	


@app.route("/run_gatk_selectvariants", methods=["POST"])
def run_gatk_selectvariants():
    rjson = request.get_json()
    ## get vals
    path_to_vcf = rjson["path_to_vcf"]
    path_to_dict_fa = rjson["path_to_dict_fa"]
    path_to_selected_vcf = rjson["path_to_selected_vcf"]
    str_selected_variant_type = rjson["select_variant_type"] #'SNP' or 'indel'
    path_to_log = rjson["log_path"]
    
    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar SelectVariants -R {path_to_dict_fa} -V {path_to_vcf} -select-type {str_selected_variant_type} -O {path_to_selected_vcf}"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)

@app.route("/run_gatk_variantfiltration", methods=["POST"])
def run_gatk_variantfiltration():
    rjson = request.get_json()
    ## get vals
    path_to_vcf = rjson["path_to_vcf"]
    path_to_dict_fa = rjson["path_to_dict_fa"]
    path_to_filtered_vcf = rjson["path_to_filtered_vcf"]
    path_to_log = rjson["log_path"]
    
    log = open(path_to_log, "w")
    
    cmd = (
      f"/usr/bin/java -jar /usr/lib/gatk/gatk-package-4.2.6.1-local.jar VariantFiltration --R {path_to_dict_fa} --V {path_to_vcf} --window 35 --cluster 3 --filter-name 'FS' --filter 'FS>30.0' --filter-name 'QD' --filter 'QD<2.0' -O {path_to_filtered_vcf}"
      )
    
    proc = subprocess.Popen(cmd, shell=True, stdout=log, stderr=log) 
    
    log.flush()
    return str(proc.pid + 1)
  
  
