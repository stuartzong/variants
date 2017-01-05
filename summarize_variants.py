#! /usr/bin/env python
"""
This pipeline does the following:
-> Combine DNA and RNA variants from mpileup, mutseq, strelka pipeline etc,
-> Summarize variants by gene and variant level,
-> Rerun samtools mpileup at default setting for all variant positions,
-> Parse mpileup results to get base count info,
-> Combining base count info for tumour-normal pairs if applicable,
-> Adjust tumour base count according to tumour content,
-> Perform Fisher Exact test to determine p value,
-> Annotate variant in summary with base counts and cosmic,
-> Indicate the caller(s) for each variant.
"""


import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from pprint import pprint
from itertools import islice
import operator
import ConfigParser
from collections import defaultdict
import headers as HEADER


from jinja2 import Environment, FileSystemLoader
import logging
import colorlog

logger = colorlog.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter('%(log_color)s%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(handler)


def get_files(bam_vcf_files, input_headers):
    """ 
    Dictionary: holding all file paths
    {patient ->
             {status ->
                     {file_identifier -> file_path}}}  
    """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        headers = set(records.fieldnames)
        input_headers = set(input_headers)
        # check if all mandatory columns prsent in input file
        if input_headers.issubset(headers):
            logger.info("input file has all mandatory columns! Continue...")
            for line in records:
                patient = line['patient']
                status = line['status']
                if patient not in patient_files:
                    patient_files[patient] = {}
                if status not in patient_files[patient]:
                    patient_files[patient][status] = line
                else:
                    logger.error('Duplicate status entries!'\
                                 'One tissue sequenced multiple times?')
                    sys.exit()
        else:
            logger.critical("Input file doesn't contain all mandatory headers!\n %s"
                            % list(input_headers))
            sys.exit()
        return patient_files


def example_display(n, iterable):
    """ 
    Return first n items of the iterable as a list! 
    """
    return list(islice(iterable, n))


def group_readable(file_path):
    m = os.stat(file_path).st_mode
    #other_execute  = bool(m & 0001)
    #other_write = bool(m & 0002)
    #other_read  = bool(m & 0004)
    group_read = bool(m & stat.S_IRGRP)
    return group_read


def check_files(files):
    missing_files = []
    for file in files:
        short_name = file.split("/")[-1]
        if os.path.exists(file):
            readable = group_readable(file)
            if readable:
                logger.info("%s: OK" % short_name)
            else:
                missing_files.append(file)
    if missing_files:
        logger.error("ERROR: The following files missing")
        for file in missing_files:
            logger.info("Missing %s" % file)
        sys.exit()


def delete_files_by_extension(extension):
    files = [ f for f in os.listdir(".") if f.endswith(tuple(extension))]
    for f in files:
        os.remove(f)


def delete_files(files):
    for f in files:
        os.remove(f)


def run_fisher_exact_test(Rscript_path, r_script):
    p = subprocess.Popen('%s/Rscript %s' % (Rscript_path, r_script),
                         shell=True, stdout=subprocess.PIPE)
    output,  err = p.communicate()
    #logger.info(output)


def count_ref_alt_bases(snp_list, base_count_file, af_file):
   """ Generate allele frequency dict to append af info to summary file """
   snps = dict()
   print "snp_list is %s" % snp_list
   for snp in snp_list:
       sl = snp.rstrip().split(':')
       chr, pos, ref, alt = sl
       key = ":".join([chr, pos, ref])
       if key not in snps:
           snps[key] = [alt]
       else:
           snps[key].append(alt)
   # parse mpileup output to make af dict
   afs = dict()
   logger.info("Parsing mpileup file: %s!" % base_count_file)
   with open (base_count_file) as fh:
       for line in fh:
           sl = line.strip().split('\t')
           snp, cov = sl[:2]
           pos_ref = sl[0].split(":")
           chr, pos, ref = pos_ref
           a, c, g, t, n = sl[2:7]
           d = {"A": a, "C": c, "G": g, "T": t, "N": n}
           if ref in d:
               ref_count = d[ref]
               d.pop(ref, None)
               alts = snps[snp]
               for alt in alts:
                   alt_count = d[alt]
                   adj_cov = int(alt_count)+int(ref_count)
                   print "cov = %s, adj_cov = %s" % (cov, adj_cov)
                   if (adj_cov != 0):
                       af = str("{:.2f}".format(float(alt_count)/adj_cov))
                   else:
                       af = str("{:.2f}".format(float(alt_count)/int(cov)))
                   value = ":".join([alt, cov, ref_count, alt_count, af])
                   try:
                       afs[snp].append(value)
                   except KeyError:
                       afs[snp] = [value]
   # logger.info(afs)
   with open(af_file,  'wb') as opf:
       writer = csv.writer(opf,  delimiter='\t', lineterminator="\n")
       for snp in afs:
           sp_snp = snp.split(':')
           for af in afs[snp]:
               sp_af = af.split(':')
               content = sp_snp + sp_af
               writer.writerow(content)
   return afs


def parse_snv_vcf(vcf_file, caller, impacts, gene_pats, gene_variants,
              snv_pos, snv_list, patient_status):
    logger.info("Parsing %s, Variant caller is: %s! " %
                (vcf_file.split("/")[-1], caller))
    with open(vcf_file,  'r') as fh:
        ## Following line is useful for testing the script
        # fh = [next(fh) for x in xrange(350)]
        for line in fh:
            if not line.startswith("#"):
                sl = line.strip().split("\t")
                Info = sl[7]
                if (not Info.startswith("INDEL") and
                    "EFF=" in Info and
                    any(x in Info for x in impacts)):
                    splitInfo = Info.split(";")
                    # get gmaf value
                    try:
                       gmaf = [j for j in splitInfo if
                           j.startswith("GMAF=")][0].split("=")[1]
                    except:
                       gmaf = "gmaf_unknown"
                    transcripts = [i for i in splitInfo if
                        i.startswith("EFF=")][0].split("=")[1].split(",")
                    rs_cos_id = sl[2]
                    chr, start, ref = sl[0], sl[1], sl[3]
                    alts = sl[4].split(",")
                    snv_pos.append("\t".join([chr, start]))
                    for alt in alts:
                        snv = ":".join([chr, start, ref, alt])
                        snv_list.append(snv)
                        snp_id_tmp = rs_cos_id.split(";")[0]
                        if snp_id_tmp.startswith("rs"):
                           snp_id = snp_id_tmp
                        else:
                           snp_id = "novel_snp"
                        try:
                           cosmic_id = rs_cos_id.split(";")[1]
                        except:
                           if snp_id_tmp.startswith("COS"):
                               cosmic_id = snp_id_tmp
                           else:
                               cosmic_id = "not_in_cosmic64"
                        selected_transcripts = [k for k in transcripts
                                                if any(x in k for x in impacts)]
                        # pick the first transcript to get the gene name
                        trans_details = selected_transcripts[0]
                        gene = trans_details.split("|")[5]
                        if (gene == ""):
                            gene = "INTERGENIC"
                        # make gene -> patient dict to count how many patients
                        # have variant in this gene
                        try:
                           gene_pats[gene].append(patient_status)
                        except:
                           gene_pats[gene] = [patient_status]
                        details = ":".join([snp_id, cosmic_id, gmaf, trans_details, caller])
                        # put all relevant info into a triple nested dictionary:
                        # gene -> snv -> patient -> details
                        try:
                           gene_variants[gene][snv][patient_status].append(details)
                        except:
                           gene_variants[gene][snv][patient_status] = [details]
    return [gene_pats, gene_variants, snv_pos, snv_list]


def summarize_snvs(patient_files, sum_header_wn_tmp, hm_sum_wn_tmp, hm_impacts):
    """
    Generate triple nested dictionary for each
    Combine all the SNV positions from each pileline and patient status for pileup
    Positions = mpileup + strelka + primary + relaspe
    Get snp list for all patients
    Summarize all SNVs with high or moderate or low or modifier impact
    Annotate as if a SNV in mpileup or/and strelka
    pool all snv positions from various vcf files: single, paired, mutseq, strelka
    if a variant exists in one tissue, time point, check all samples from the same patient
    """
    tree = lambda: defaultdict(tree)
    hm_variants_dict = tree()
    hm_gene_dict = dict()
    patient_snv = dict()
    patient_callers = dict()
    snv_pos_files = []
    for patient in patient_files:
        logger.info("Parsing all vcf files related to: %s >>>>>>" % patient)
        snv_pos_file = ".".join([patient, "vcf.snp.pos"])
        snv_pos_files.append(snv_pos_file)
        snv_pos = []
        snv_list = []
        for status in patient_files[patient]:
            if ("normal" not in status):
                status_vcfs = patient_files[patient][status]
                logger.info("Parsing vcf files related to: %s, %s ......" %
                            (patient, status))
                patient_status = "_".join([patient, status])
                for key in status_vcfs:
                    caller = key.replace('_vcf', '').replace('_snv', '')
                    vcf = status_vcfs[key]
                    if 'vcf' in key and 'indel' not in key and vcf != 'NA':
                        # print(caller, vcf)
                        try:
                            patient_callers[patient_status].append(caller)
                        except KeyError:
                            patient_callers[patient_status] = [caller]
                        two_dicts = parse_snv_vcf(vcf, caller, hm_impacts,
                                              hm_gene_dict, hm_variants_dict,
                                              snv_pos, snv_list, patient_status)
                        hm_gene_dict = two_dicts[0]
                        hm_variants_dict = two_dicts[1]
                        snv_pos = two_dicts[2]
                        snv_list = two_dicts[3]
        snv_pos = list(set(snv_pos))
        snv_list = list(set(snv_list))
        patient_snv[patient] = snv_list
        with open(snv_pos_file, 'w') as pos_writer:
            for snv_coord in sorted(snv_pos):
                pos_writer.write(snv_coord)
                pos_writer.write("\n")
    #logger.info("Example patient snv dict are: \n")
    # n_items = example_display(2, patient_snv.iteritems())
    # logger.info(n_items)
    write_summary(hm_variants_dict, hm_gene_dict, hm_sum_wn_tmp,
                  sum_header_wn_tmp, patient_callers, patient_files)
    return [patient_snv, snv_pos_files]


def write_summary(variants_dict, gene_dict, summary_file,
                  sum_header_wn_tmp, patient_callers, patient_files):
    writer = open(summary_file, 'w')
    writer.write("\t".join(sum_header_wn_tmp))
    writer.write("\n")
    logger.info("Writing variants_dict into summary file: %s" % summary_file)
    for gene in variants_dict:
        short_patient = [i.split("_")[0] for i in gene_dict[gene]]
        num_patients_gene_level = len(list(set(short_patient)))
        numsnvs = len(variants_dict[gene].keys())
        for snv in variants_dict[gene]:
            patient_statuses_snv_level = list(set([i for i in
                                          variants_dict[gene][snv].keys()]))
            patients_snv_level = [i.split("_")[0] for
                                  i in patient_statuses_snv_level]
            num_patients_snv_level = len(patients_snv_level)
            snv_details = []
            # logger.info("patients_snv_level are: %s" % patients_snv_level)
            # all patient status combinations
            all_patient_statuses = []
            for pat in patients_snv_level:
                for sta in patient_files[pat]:
                    if not re.search(r'normal', sta, re.M|re.I):
                        all_patient_statuses.append( "_".join([pat, sta]))
            all_patient_statuses = list(set(all_patient_statuses))
            patient_statuses_not_called = list(set(all_patient_statuses)-
                                               set(patient_statuses_snv_level))
            # logger.info("all_patient_statuses is %s, %s, %s, %s" %
            #             (gene, all_patient_statuses, patient_statuses_snv_level,
            #              patient_statuses_not_called))
            for patient in variants_dict[gene][snv]:
                chr, start, ref, alt = snv.split(':')
                # anno_details: rs8473:COSM1128106:0.4748:NON_SYNONYMOUS_
                # CODING(MODERATE|MISSENSE|Aaa/Gaa|K2857E|2896|MKI67|protein_
                # coding|CODING|ENST00000368653|13|1):pileup
                anno_details_all = variants_dict[gene][snv][patient]
                callers = [i.split(":")[-1] for i in anno_details_all]
                anno_details = list(set(anno_details_all))[0].split(":")
                rs_id = anno_details[0]
                cosmic_id = anno_details[1]
                gmaf = anno_details[2]
                snp_eff = anno_details[3]
                # see if variant called by both mpileup and strelka etc
                callers_run = patient_callers[patient]
                (in_paired_mpileup,
                 in_mutseq, in_strelka,
                 in_DNA_single_mpileup,
                 in_RNA_single_mpileup) = determine_callers(callers_run, callers)
                content = "\t".join([gene, str(num_patients_gene_level),
                                     str(numsnvs), chr, start, ref, alt,
                                     str(num_patients_snv_level), patient,
                                     rs_id, gmaf, cosmic_id, snp_eff,
                                     in_DNA_single_mpileup, in_paired_mpileup,
                                     in_mutseq, in_strelka, in_RNA_single_mpileup])
                writer.write(content)
                writer.write("\n")
            # include variant for all statuses as long it is called in status
            for patient_status in patient_statuses_not_called:
                # patient = patient_status
                callers_run = patient_callers[patient_status]
                callers = []
                (in_paired_mpileup,
                 in_mutseq, in_strelka,
                 in_DNA_single_mpileup,
                 in_RNA_single_mpileup) = determine_callers(callers_run, callers)
                content = '\t'.join([gene, str(num_patients_gene_level),
                                     str(numsnvs), chr, start, ref, alt,
                                     str(num_patients_snv_level), patient_status,
                                     rs_id, gmaf, cosmic_id, snp_eff,
                                     in_DNA_single_mpileup, in_paired_mpileup,
                                     in_mutseq, in_strelka, in_RNA_single_mpileup])
                writer.write(content)
                writer.write("\n")


def determine_callers(callers_run, snv_callers):
    # callers_run are callers that have run for a patient_status bam
    # snv_callers are the callers reporting this snv
    in_paired_mpileup = "not_in_paired_mpileup"
    in_mutseq = "not_in_mutseq"
    in_strelka = "not_in_strelka"
    in_DNA_single_mpileup = "not_in_DNA_single_mpileup"
    in_RNA_single_mpileup = "not_in_RNA_single_mpileup"
    if 'strelka' in callers_run:
        if 'strelka' in snv_callers:
            in_strelka = "in_strelka"
    else:
        # print 'ooooooooooo', callers_run, snv_callers
        in_strelka = "strelka_not_run"
    if 'paired_mpileup' in callers_run:
        if 'paired_mpileup' in snv_callers:
            in_paired_mpileup = 'in_paired_mpileup'
    else:
        in_paired_mpileup = 'paired_mpileup_not_run'
    if 'DNA_single' in callers_run:
        if 'DNA_single' in  snv_callers:
            in_DNA_single_mpileup = 'in_DNA_single_mpileup'
    else:
        in_DNA_single_mpileup = 'DNA_single_mpileup_not_run'
    if 'RNA_single' in callers_run:
        if 'RNA_single' in  snv_callers:
            in_RNA_single_mpileup = 'in_RNA_single_mpileup'
    else:
        in_RNA_single_mpileup = 'RNA_single_mpileup_not_run'
    if 'mutseq' in callers_run:
        if 'mutseq' in snv_callers:
            in_mutseq = 'in_mutseq'
    else:
        in_mutseq = 'mutseq_not_run'
    return (in_paired_mpileup, in_mutseq, in_strelka,
            in_DNA_single_mpileup, in_RNA_single_mpileup)


def calculate_af (cov, refC, altC):
    # ref + alt <= cov because low qaulity bases filtered
    adj_cov = int(refC) + int(altC)
    if (adj_cov != 0):
        af = "{:.2f}".format(float(altC)/adj_cov)
    elif (int(cov) != 0):
        af = "{:.2f}".format(float(altC)/int(cov))
    else:
        af = 0.00
    return af


def summarize_indels(patient_files, sum_header_wn_tmp,
                     hm_sum_wn_tmp, hm_impacts):
    tree = lambda: defaultdict(tree)
    hm_variants_dict = tree()                                         
    hm_gene_dict = dict()
    hm_indel_dict = dict()
    patient_callers = dict()
    modifier_indel_dict = dict()
    for patient in patient_files:
        logger.info("Parsing vcf files related to %s >>>>>>" % patient)
        for status in patient_files[patient]:
            if ("normal" not in status):
                status_vcfs = patient_files[patient][status]
                logger.info("Parsing vcf related to %s %s ......" %
                            (patient, status))
                strelka_indel_vcf = status_vcfs['strelka_indel_vcf']
                DNA_single_vcf = status_vcfs['DNA_single_vcf']
                RNA_single_vcf = status_vcfs['RNA_single_vcf']
                paired_mpileup_vcf = status_vcfs['paired_mpileup_vcf']
                #mutseq does not have indel result
                # mutseq_snv_vcf = status_vcfs['mutseq_snv_vcf']
                patient_status = "_".join([patient, status])
                patient_status = "_".join([patient, status])
                for key in status_vcfs:
                    # do not parse mutseq or strelka snv vcfs
                    caller = key.replace('_vcf', '').replace('_indel', '')
                    vcf = status_vcfs[key]
                    if 'vcf' in key and 'strelka_indel' in key and vcf != 'NA':
                        try:
                            patient_callers[patient_status].append(caller)
                        except KeyError:
                            patient_callers[patient_status] = [caller]
                        two_dicts = parse_strelka_indel(vcf, caller, hm_impacts,
                                              hm_gene_dict, hm_variants_dict,
                                              patient_status)
                        hm_gene_dict, hm_variants_dict = two_dicts
                    elif ('vcf' in key and 'strelka' not in key and
                         'mutseq' not in key and 'snv' not in key and
                          vcf != 'NA'):
                        try:
                            patient_callers[patient_status].append(caller)
                        except KeyError:
                            patient_callers[patient_status] = [caller]
                        two_dicts = parse_mpileup_indel(vcf, caller, hm_impacts,
                                              hm_gene_dict, hm_variants_dict,
                                              patient_status)
                        hm_gene_dict, hm_variants_dict = two_dicts
    write_indel_summary(hm_variants_dict, hm_gene_dict, hm_sum_wn_tmp,
                        sum_header_wn_tmp, patient_callers, patient_files)


def parse_mpileup_indel(vcf_file, caller, impacts,
                        gene_patients, gene_variants, patient_status):
    logger.info("%s variant caller is: %s! " %  (vcf_file.split("/")[-1], caller))
    with open(vcf_file,  'r') as fh:
        # Following line is useful for testing the script
        # pileup header 153 lines
        # fh = [next(fh) for x in xrange(200)]
        for line in fh:
            if not line.startswith("#"):
                sl = line.strip().split("\t")
                Info = sl[7]
                if (Info.startswith("INDEL") and
                    "EFF=" in Info and
                    any(x in Info for x in impacts)):
                    #logger.info(Info, impacts)
                    splitInfo = Info.split(";")
                    cov =  [i for i in splitInfo if
                            i.startswith("DP=")][0].split("=")[1]
                    refC_altC =  [i for i in splitInfo if
                                  i.startswith("DP4=")][0].split("=")[1].split(",")
                    tmp_refC = refC_altC[0:2]
                    tmp_altC = refC_altC[2:]
                    refC = sum(int(i) for i in tmp_refC)
                    altC = sum(int(i) for i in tmp_altC)
                    af = calculate_af (cov, refC, altC)
                    mpileup_cov = [cov, str(refC), str(altC), str(af)]
                    strelka_cov = ["na", "na", "na", "na", "na", "na", "na", "na"]
                    cov_info = "\t".join(mpileup_cov + strelka_cov)
                    try:
                        gmaf = [j for j in splitInfo if
                                j.startswith("GMAF=")][0].split("=")[1]
                    except:
                       gmaf = "gmaf_unknown"
                    transcripts = [i for i in splitInfo if
                                   i.startswith("EFF=")][0].split("=")[1].split(",")
                    rs_cos_id = sl[2]
                    chromosome, start, ref = sl[0], sl[1], sl[3]
                    alts = sl[4].split(",")
                    for alt in alts:
                        indel = ":".join([chromosome, start, ref, alt])
                        snp_id_tmp = rs_cos_id.split(";")[0]
                        if snp_id_tmp.startswith("rs"):
                           snp_id = snp_id_tmp
                        else:
                           snp_id = "novel_snp"
                        try:
                           cosmic_id = rs_cos_id.split(";")[1]
                        except:
                           if snp_id_tmp.startswith("COS"):
                               cosmic_id = snp_id_tmp
                           else:
                               cosmic_id = "not_in_cosmic64"
                        selected_transcripts = [k for k in transcripts if
                                                any(x in k for x in impacts)]
                        # pick the first transcript to get the gene name
                        trans_details = selected_transcripts[0]
                        gene = trans_details.split("|")[5]
                        if (gene == ""):
                            gene = "INTERGENIC"
                        # make gene -> patient dict to count
                        # how many patients have variant in this gene
                        try:
                           gene_patients[gene].append(patient_status)
                        except:
                           gene_patients[gene] = [patient_status]
                        details = ":".join([snp_id, cosmic_id, gmaf,
                                            trans_details, cov_info, caller])
                        # put all relevant info into a triple nested dictionary:
                        # gene -> snv -> patient -> details
                        try:
                           gene_variants[gene][indel][patient_status].append(details)
                        except:
                           gene_variants[gene][indel][patient_status] = [details]
    return [gene_patients, gene_variants]


def parse_strelka_indel(vcf_file, caller, impacts, gene_pats,
                        gene_variants, patient_status):
    logger.info("%s variant caller is: %s!" %  (vcf_file.split("/")[-1], caller))
    with open(vcf_file,  'r') as fh:
        # Following line is useful for testing the script
        # strelka header lines 324
        # fh = [next(fh) for x in xrange(300)]
        for line in fh:
            if (not line.startswith("#")):
                sl = line.strip().split("\t")
                Info = sl[7]
                if ("EFF=" in Info) and any(x in Info for x in impacts):
                    # use tier1 counts, which is with mapping quality >40
                    normal_info = sl[9].split(":")
                    tumor_info = sl[10].split(":")
                    splitInfo = Info.split(";")
                    nor_cov =  normal_info[0]
                    nor_refC = normal_info[2].split(",")[0]
                    nor_indelC = normal_info[3].split(",")[0]
                    tum_cov =  tumor_info[0]
                    tum_refC = tumor_info[2].split(",")[0]
                    tum_indelC = tumor_info[3].split(",")[0]
                    #nor_adj_cov = int(nor_refC) + int(nor_indelC)
                    #tum_adj_cov = int(tum_refC) + int(tum_indelC)
                    nor_af = calculate_af (nor_cov, nor_refC, nor_indelC)
                    tum_af = calculate_af (tum_cov, tum_refC, tum_indelC)
                    mpileup_cov = ["na", "na", "na", "na"]
                    strelka_cov = [nor_cov, nor_refC, nor_indelC,
                                   str(nor_af), tum_cov, tum_refC,
                                   tum_indelC, str(tum_af)]
                    cov_info = "\t".join(mpileup_cov + strelka_cov)
                    try:
                        gmaf = [j for j in splitInfo if
                                j.startswith("GMAF=")][0].split("=")[1]
                    except:
                       gmaf = "gmaf_unknown"
                    transcripts = [i for i in splitInfo if
                        i.startswith("EFF=")][0].split("=")[1].split(",")
                    rs_cos_id = sl[2]
                    chr, start, ref = sl[0], sl[1], sl[3]
                    alts = sl[4].split(",")
                    for alt in alts:
                        indel = ":".join([chr, start, ref, alt])
                        snp_id_tmp = rs_cos_id.split(";")[0]
                        if snp_id_tmp.startswith("rs"):
                           snp_id = snp_id_tmp
                        else:
                           snp_id = "novel_snp"
                        try:
                           cosmic_id = rs_cos_id.split(";")[1]
                        except:
                           if snp_id_tmp.startswith("COS"):
                               cosmic_id = snp_id_tmp
                           else:
                               cosmic_id = "not_in_cosmic64"
                        selected_transcripts = [k for k in transcripts
                                                if any(x in k for x in impacts)]
                        trans_details = selected_transcripts[0]
                        gene = trans_details.split("|")[5]
                        if (gene == ""):
                            gene = "INTERGENIC"
                        try:
                           gene_pats[gene].append(patient_status)
                        except:
                           gene_pats[gene] = [patient_status]
                        details = ":".join([snp_id, cosmic_id, gmaf,
                                            trans_details, cov_info, caller])
                        try:
                           gene_variants[gene][indel][patient_status].append(details)
                        except:
                           gene_variants[gene][indel][patient_status] = [details]
    return [gene_pats, gene_variants]



def write_indel_summary(variants_dict, gene_dict, summary_file,
                        sum_header_wn_tmp, patient_callers, patient_files):
    writer = open(summary_file, 'w')
    writer.write("\t".join(sum_header_wn_tmp))
    writer.write("\n")
    logger.info("Writing variants_dict into %s" % summary_file)
    for gene in variants_dict:
        short_patient = [i.split("_")[0] for i in gene_dict[gene]]
        num_patients_gene_level = len(list(set(short_patient)))
        numSNVs = len(variants_dict[gene].keys())
        for snv in variants_dict[gene]:
            patients_snv_level = [i.split("_")[0] for i in
                                  variants_dict[gene][snv]]
            num_patients_SNV_level = len(list(set(patients_snv_level)))
            snv_details = []
            for patient in variants_dict[gene][snv]:
                # snv_details: 10:64927823:C:G
                chr, start, ref, alt = snv.split(":")
                # #anno_details: rs8473:COSM1128106:0.4748:NON_SYNONYMOUS_
                # CODING(MODERATE|MISSENSE|Aaa/Gaa|K2857E|2896|MKI67|protein_
                #        coding|CODING|ENST00000368653|13|1):pileup
                anno_details_all = variants_dict[gene][snv][patient]
                
                callers = [i.split(":")[-1] for i in anno_details_all]
                if ("strelka" in callers):
                    anno_details = [i for i in anno_details_all
                                    if "strelka" in i][0].split(":")
                else:
                    # callers = [i.split(":")[-1] for i in anno_details_all]
                    # pick highest ref+alt if called both in DNA and RNA
                    multi_cov = [i.split(':')[4] for i in anno_details_all]
                    # print('****************************') 
                    # print(multi_cov, anno_details_all)
                    multi_ref_alt= [int(i.split('\t')[1]) +
                                    int(i.split('\t')[2])
                                    for i in multi_cov]
                    index, value = max(enumerate(multi_ref_alt),
                                       key=operator.itemgetter(1))
                    anno_details = anno_details_all[index].split(':')
                (rs_id, cosmic_id, gmaf,
                 snp_eff, cov_info ) = anno_details[0:5]
                # see if variant called by both mpileup and strelka
                callers_run = patient_callers[patient]
                (in_paired_mpileup,
                 in_mutseq, in_strelka,
                 in_DNA_single_mpileup,
                 in_RNA_single_mpileup) = determine_callers(callers_run,
                                                            callers)
                content = "\t".join([gene, str(num_patients_gene_level),
                                     str(numSNVs), chr, start, ref, alt,
                                     str(num_patients_SNV_level), patient,
                                     rs_id, gmaf, cosmic_id, snp_eff,
                                     in_DNA_single_mpileup,
                                     in_paired_mpileup,
                                     in_strelka, in_RNA_single_mpileup,
                                     cov_info])

                writer.write(content)
                writer.write("\n")


def adjust_allele_frequency(type, com_af_files,patient_files):
    tc_key = "_".join([type, "tc"])
    for file in com_af_files:
        file = file.strip()
        sp_f = re.split('[.]',file)
        patient =sp_f[0]
        status = sp_f[1]
        pat_status = "_".join([patient, status])
        af_aj_file = file + ".adjusted"
        writer = csv.writer(open(af_aj_file, 'wb'),delimiter='\t', lineterminator="\n")
        with open (file) as handle:
            for line in handle:
                sl = line.rstrip().split('\t')
                n_totC = int(sl[4])
                n_refC = int(sl[5])
                n_altC = int(sl[6])
                t_totC = int(sl[9])
                t_refC = int(sl[10])
                t_altC = int(sl[11])
                n_af = float(sl[7])
                t_af = float(sl[12])
                # skip adjustment for second low frequency allele
                tlow_totC = t_refC + t_altC
                tc_tmp = patient_files[patient][status][tc_key]
                if ( n_totC > 0 and t_totC > 0 and
                     (tlow_totC > 0.7*t_totC) and
                     ("unknow" not in tc_tmp) ):
                    nlz_factor = float(float(t_totC)/n_totC)
                    tc = float(int(tc_tmp)/100.0)
                    n_nlz_refC = n_refC*nlz_factor
                    n_nlz_altC = n_altC*nlz_factor
                    t_aj_ref = float(t_refC-n_nlz_refC*(1-tc))/tc
                    t_aj_alt = float(t_altC-n_nlz_altC*(1-tc))/tc
                    if (t_aj_ref >= 0 and t_aj_alt >= 0):
                        t_aj_tot = t_aj_ref + t_aj_alt
                        if t_aj_tot != 0:
                            t_aj_af = t_aj_alt/t_aj_tot
                            if t_aj_af >= 1:
                                t_aj_af =1
                            elif t_aj_af<0:
                                t_aj_af = 0
                        else:
                            t_aj_tot, t_aj_ref = t_totC, t_refC
                            t_aj_alt, t_aj_af = t_altC, t_af
                    elif (t_aj_ref < 0 and t_aj_alt >= 0):
                        t_aj_alt = t_aj_alt - t_aj_ref
                        t_aj_ref = 0
                        t_aj_tot = t_aj_ref + t_aj_alt
                        if t_aj_tot != 0:
                            t_aj_af = t_aj_alt/t_aj_tot
                            if t_aj_af >= 1:
                                t_aj_af =1
                            elif t_aj_af < 0:
                                t_aj_af = 0
                        else:
                            t_aj_tot, t_aj_ref = t_totC, t_refC
                            t_aj_alt, t_aj_af = t_altC, t_af
                    elif (t_aj_ref >= 0 and t_aj_alt < 0):
                        t_aj_ref = t_aj_ref - t_aj_alt
                        t_aj_alt = 0
                        t_aj_tot = t_aj_ref + t_aj_alt
                        if t_aj_tot != 0:
                            t_aj_af = t_aj_alt/t_aj_tot
                            if t_aj_af >= 1:
                                t_aj_af = 1
                            elif t_aj_af < 0:
                                t_aj_af = 0
                        else:
                            t_aj_tot, t_aj_ref = t_totC, t_refC
                            t_aj_alt, t_aj_af = t_altC, t_af
                else:
                    t_aj_tot, t_aj_ref = t_totC, t_refC
                    t_aj_alt, t_aj_af = t_altC, t_af
                    try:
                        tc = float(int(tc_tmp)/100.0)
                    except:
                        tc = "tc_unknown"
                t_aj_tot = "{:.1f}".format( float(t_aj_tot) )
                t_aj_ref = "{:.1f}".format( float(t_aj_ref) )
                t_aj_alt = "{:.1f}".format( float(t_aj_alt) )
                t_aj_af = "{:.2f}".format( float(t_aj_af) )
                writer.writerow([sl[0],sl[1],sl[2],sl[3],
                                 sl[4],sl[5],sl[6],sl[7],sl[8],
                                 t_totC, t_refC, t_altC, t_af,
                                 t_aj_tot,t_aj_ref,t_aj_alt,t_aj_af,tc])
        logger.info("Finishing af adjustment for: %s" % (pat_status))


def make_patient_wn_af_dict(af_files):
   patient_af_dict = dict()
   for file in af_files:
       if (file != "NA"):
           patient = "_".join(re.split('[.]', file)[:2])
           with open(file) as f:
              next( f )
              for line in f:
                 sl = [i.strip() for i in line.rstrip().split('\t')]
                 chr, pos, ref, alt = sl[1], sl[2], sl[3], sl[4]
                 key = "_".join([patient, chr, pos, ref, alt])
                 n_cov, n_ref, n_alt, n_af = sl[5], sl[6], sl[7], sl[8]
                 t_cov, t_ref, t_alt, t_af = sl[10], sl[11], sl[12], sl[13]
                 n_af = "{:.2f}".format( float(n_af) )
                 t_af = "{:.2f}".format( float(t_af) )
                 ajt_cov, ajt_ref, ajt_alt, ajt_af = sl[14], sl[15], sl[16], sl[17]
                 ajt_cov = "{:.1f}".format( float(ajt_cov) )
                 ajt_ref = "{:.1f}".format( float(ajt_ref) )
                 ajt_alt = "{:.1f}".format( float(ajt_alt) )
                 ajt_af = "{:.2f}".format( float(ajt_af) )
                 p_value, tc = sl[19], sl[18]
                 value = ",".join([n_cov, n_ref, n_alt, n_af,
                                   t_cov, t_ref, t_alt, t_af,
                                   ajt_cov, ajt_ref, ajt_alt, ajt_af,
                                   p_value])
                 try:
                   patient_af_dict[key].append(value)
                 except KeyError:
                   patient_af_dict[key] = [value]
   return patient_af_dict


def make_patient_non_af_dict(af_files):
   patient_af_dict = dict()
   for file in af_files:
       if (file != "NA"):
           patient = "_".join(re.split('[.]', file)[:2])
           with open(file) as f:
              for line in f:
                 sl = [i.strip() for i in line.rstrip().split('\t')]
                 chr, pos, ref, alt = sl[0], sl[1], sl[2], sl[3]
                 key = "_".join([patient, chr, pos, ref, alt])
                 t_cov, t_ref, t_alt, t_af = sl[4], sl[5], sl[6], sl[7]
                 value = ",".join([t_cov, t_ref, t_alt, t_af])
                 try:
                   patient_af_dict[key].append(value)
                 except KeyError:
                   patient_af_dict[key] = [value]
   return patient_af_dict


def qsub_scripts(scripts):
    """ qsub scripts """
    wkdir = os.getcwd()
    for script in scripts:
        p = subprocess.Popen('ssh tachpc \"cd %s;  qsub %s\"' %
                             (wkdir, script),  shell=True,
                             stdout=subprocess.PIPE)
        output,  err = p.communicate()


def parse_pileup_output(patient_files, patient_snv_wn, patient_snv_non):
    af_files = dict()
    com_af_files = dict()
    af_dicts = dict()
    data_types = ["DNA", "RNA"]
    for type in data_types:
        af_files[type] = []
        com_af_files[type] = []
        af_dicts[type] = dict()
    for patient in patient_files:
        if ("normal" in patient_files[patient]):
            snp_list = patient_snv_wn[patient]
        else:
            snp_list = patient_snv_non[patient]
        # get base count  and af for all samples of a patient
        for status in patient_files[patient]:
            bam_types = dict()
            DNA_bam = patient_files[patient][status]['DNA_bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            bam_types["DNA"] = DNA_bam
            bam_types["RNA"] = RNA_bam
            for type in bam_types:
                if (bam_types[type] != "NA"):
                    #for type in data_types:
                    base_count_file  = ".".join([patient, status,
                                                 type, "pileup.AFcounts"])
                    af_file = ".".join([patient, status,
                                        type, "pileup.AFcounts.af"])
                    af_files[type].append( af_file )
                    # variant bp count dict for primary, relapse, and normal
                    af_dict  = count_ref_alt_bases(snp_list,
                                                   base_count_file,
                                                   af_file)
                    # '5:131892979:G': ['G:88:1:1:0.01']
                    try:
                        #logger.info(patient, status)
                        af_dicts[type][patient][status] = af_dict
                    except KeyError:
                        if patient not in af_dicts[type]:
                            af_dicts[type][patient] = {}
                        if status not in af_dicts[type][patient]:
                            af_dicts[type][patient][status] = af_dict
        nor_bams = dict()
        if ("normal" in patient_files[patient]):
            nor_bams["DNA"] = patient_files[patient]["normal"]['DNA_bam']
            nor_bams["RNA"] = patient_files[patient]["normal"]['RNA_bam']
        else:
            nor_bams["DNA"] = "NA"
            nor_bams["RNA"] = "NA"
        for type in nor_bams:
            if (nor_bams[type] != "NA"):
                n_af_dict = af_dicts[type][patient]["normal"]
                for status in patient_files[patient]:
                    # if tumor bam exist
                    bam_key = "_".join([type, "bam"])
                    t_bam = patient_files[patient][status][bam_key]
                    if (not status == "normal") and (t_bam != "NA"):
                        logger.info("merging %s %s %s with its normal!" %
                                    (patient, status, type))
                        com_af_file = ".".join([patient, status,
                                                type, "af.combined"])
                        com_af_files[type].append(com_af_file)
                        t_af_dict = af_dicts[type][patient][status]
                        with open(com_af_file,  'wb') as opf:
                            writer = csv.writer(opf,  delimiter='\t', lineterminator="\n")
                            #logger.info("snp_list is: ", snp_list)
                            for snp in snp_list:
                               sp_snp = snp.split(':')
                               chr, start, ref = sp_snp[0], sp_snp[1], sp_snp[2]
                               alt = sp_snp[3]
                               af_key = ":".join([chr, start, ref])
                               if ("no_matched_normal" not in n_af_dict):
                                   try:
                                       normal_afs = n_af_dict[af_key]
                                       normal_af = [i for i in normal_afs
                                                    if i[0] == alt][0].split(":")
                                   except KeyError:
                                       normal_af = ["N", "0", "0", "0", "0"]
                               else:
                                   alt_tag = "_".join(["no", type, "normal"])
                                   normal_af = [alt_tag, "0", "0", "0", "0"]
                               n_totC, n_refC, n_altC, n_af = normal_af[1:5]
                               try:
                                   tumour_afs = t_af_dict[af_key]
                                   tumour_af = [i for i in tumour_afs if
                                                i[0] == alt][0].split(":")
                               except KeyError:
                                   tumour_af = ["N", "0", "0", "0", "0"]
                               t_totC, t_refC, t_altC, t_af = tumour_af[1:5]
                               t_alt = tumour_af[0]
                               writer.writerow([chr, start, ref, alt,
                                                n_totC, n_refC, n_altC, n_af,
                                                t_alt,
                                                t_totC, t_refC, t_altC, t_af])
            else:
                n_af_dict = dict()
                n_af_dict["no_matched_normal"] = "True"
                logger.info("no %s %s normal bam, merging not required!" %
                            (patient, type))
    return com_af_files


def detect_cluster_job_status(completeion_file_list):
    completed = False
    for file in completeion_file_list:
        if (os.path.exists(file)):
            completed = True
        else:
            completed = False
            break
    return completed


def detect_cluster_jobs(complete_stamps):
    """ detect if job on cluster finised """
    job_status = False
    logger.info("Waiting for cluster jobs to finish!\n")
    while (not job_status):
        time.sleep(10)
        job_status = detect_cluster_job_status(complete_stamps)
    logger.info("All cluster jobs finished? %s\n" % job_status)


def make_pileup_scripts(patient_files, template_dir):
    """ genearate pileup scripts, which also parses pileup results! """
    pileup_scripts = []
    pileup_stamps = []
    for patient in patient_files:
        for status in patient_files[patient]:
            snv_positions = ".".join([patient, "vcf.snp.pos"])
            DNA_bam = patient_files[patient][status]['DNA_bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            if (DNA_bam != 'NA'):
                script = populate_template(patient, status, 'DNA',
                                     snv_positions, DNA_bam, template_dir)
                pileup_scripts.append(script)
                pileup_stamp = ".".join([patient, status, 'DNA',
                                         "pileup.complete"])
                pileup_stamps.append( pileup_stamp )
            if (RNA_bam != 'NA'):
                script = populate_template(patient, status, 'RNA',
                                     snv_positions, RNA_bam, template_dir)
                pileup_scripts.append(script)
                pileup_stamp = ".".join([patient, status, 'RNA',
                                         "pileup.complete"])
                pileup_stamps.append( pileup_stamp )

    return [pileup_scripts, pileup_stamps]


def populate_template(patient, status, data, mpileup_positions,
                bam_file, template_dir):
    name = ".".join([patient, status, data])
    pileup_script = ".".join([name, "pileup.sh"])
    mpileup_output = ".".join([name, "pileup"])
    mpileup_AFcounts  = ".".join([name, "pileup.AFcounts"])
    mpileup_stamp = ".".join([name, "pileup.complete"])
    # patient_status = "_".join([patient, status, data])
    log_file = "variant_summary.log"
    jinja2_env = Environment(loader=FileSystemLoader([template_dir]),
                             trim_blocks=True)
    template = jinja2_env.get_template('mpileup_template.sh')
    with open(pileup_script, 'wb') as opf:
        content = template.render(mpileup_output = mpileup_output,
                                  mpileup_AFcounts = mpileup_AFcounts,
                                  mpileup_stamp = mpileup_stamp,
                                  log_file = log_file,
                                  mpileup_positions = mpileup_positions,
                                  bam_file = bam_file,
                                  patient_status = name)
        opf.write(content)
        logger.info('templated {0}'.format(pileup_script))
    return pileup_script


def make_af_dict(patient_files):
    """
    Reading in af files for all patients and making af_dicts:
    {'DNA': [{'Axxxxx_unknown_10_100148176_A_G': ['2,0,2,1.00,2,0,2,1.00,2.0,0.0,2.0,1.00,1.000'],
              'Axxxxx_unknown_X_99657566_G_T': ['4,2,2,0.50,4,2,2,0.50,4.0,2.0,2.0,0.50,1.000'],
              'Axxxxx_unknown_X_99941705_T_C': ['5,0,5,1.00,5,0,5,1.00,5.0,0.0,5.0,1.00,1.000']},
             {}],
     'RNA': [{},
             {'Axxxxx_unknown_10_100148176_A_G': ['2,0,2,1.00'],
              'Axxxxx_unknown_X_99941705_T_C': ['5,0,5,1.00']}]}
    """
    AFcounts_afs_wn = dict()
    af_dicts = dict()
    types = ["DNA", "RNA"]
    for type in types:
        bam_key = "_".join([type, "bam"])
        affs_wn = []
        affs_non = []
        af_dicts[type] = []
        AFcounts_afs_wn[type] = []
        for patient in patient_files:
            # see if DNA or RNA normal BAMs available
            if ("normal" in patient_files[patient]):
                nor_bam = patient_files[patient]["normal"][bam_key]
            else:
                nor_bam = "NA"
            if (nor_bam != "NA" ):
                logger.info("%s has %s normal bam!" % (patient, type))
                for status in patient_files[patient]:
                    #bam = patient_files[patient][status][bam_key]
                    if (not status == "normal"): #and (bam != "NA"):
                        aff  = ".".join([patient, status,
                                         type, "af.combined.adjusted.fisher"])
                        affs_wn.append(aff)
                        AFcounts_af  = ".".join([patient, status,
                                                 type, "pileup.AFcounts.af"])
                        AFcounts_afs_wn[type].append( AFcounts_af )
            elif (nor_bam == "NA"):
                logger.info("%s does not have %s normal bam!" % (patient, type))
                for status in patient_files[patient]:
                    bam = patient_files[patient][status][bam_key]
                    #if (bam != "NA"):
                    if (not status == "normal") and (bam != "NA"):
                        aff  = ".".join([patient, status,
                                         type, "pileup.AFcounts.af"])
                        affs_non.append( aff )
        logger.info("%s no normal af files are: %s" % (type, affs_non))
        logger.info("%s with normal af files are: %s" % (type, affs_wn))
        af_non_dict = make_patient_non_af_dict( affs_non )
        #'PASWLN_Primary_3_50879176_A_T': ['18,14,4,0.22']
        af_wn_dict = make_patient_wn_af_dict( affs_wn )
        #'PASWAJ_Relapse_8_116635942_C_T': ['96,55,41,0.43,116,61,55,0.47,116.0,61.0,55.0,0.47,0.579,Unknown']
        af_dicts[type] = [af_wn_dict, af_non_dict]
    return [af_dicts, AFcounts_afs_wn]


def add_af_anno_for_unpaired(snv_sum, snv_sum_tmp,
                             DNA_af_wn_dict, RNA_af_wn_dict,
                             DNA_af_non_dict, RNA_af_non_dict,
                             patient_files, sum_header):
    with open(snv_sum,  'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator="\n")
        writer.writerow(sum_header)
        no_coverage = ["0", "0", "0", "0"]
        no_bam = ["na", "na", "na", "na"]
        with open(snv_sum_tmp,  'r') as snv_fh:
            records = csv.DictReader(snv_fh,  delimiter='\t')
            header = records.fieldnames
            for line in records:
                line_content = [line[i] for i in header ]
                sl = line['patient_ID'].split("_")
                patient = sl[0]
                status = '_'.join(sl[1:])
                DNA_bam = patient_files[patient][status]['DNA_bam']
                RNA_bam = patient_files[patient][status]['RNA_bam']
                DNA_tc = patient_files[patient][status]['DNA_tc']
                RNA_tc = patient_files[patient][status]['RNA_tc']
                # variant = "_".join([line[patient_id, chr, pos, ref, alt])
                variant = "_".join([line['patient_ID'],
                                    line['chromosome'],
                                    line['position'],
                                    line['ref_base'],
                                    line['alt_base']])
                # DNA af info
                try:
                    DNA_af = DNA_af_non_dict[variant][0].split(',')
                except KeyError:
                    if (DNA_bam == "NA"):
                        DNA_af = no_bam
                    else:
                        DNA_af = no_coverage
                # RNA af info
                try:
                    RNA_af = RNA_af_non_dict[variant][0].split(',')
                except KeyError:
                    if (RNA_bam == "NA"):
                        RNA_af = no_bam
                    else:
                        RNA_af = no_coverage
                final_content = line_content + DNA_af + [DNA_tc] + RNA_af + [RNA_tc]
                writer.writerow(final_content)


def add_af_anno_for_paired(snv_sum, snv_sum_tmp,
                           DNA_af_wn_dict, RNA_af_wn_dict,
                           DNA_af_non_dict, RNA_af_non_dict,
                           patient_files, sum_header):
    no_bam = ["na", "na", "na", "na", "na", "na", "na",
              "na", "na", "na", "na", "na", "na"]
    no_coverage = ["0", "0", "0", "0", "0", "0", "0",
                   "0", "0", "0", "0", "0", "0"]
    nas = ["na", "na", "na", "na"]
    zeros = ["0", "0", "0", "0"]
    with open(snv_sum,  'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator="\n")
        writer.writerow(sum_header)
        with open(snv_sum_tmp,  'r') as snv_fh:
            records = csv.DictReader(snv_fh,  delimiter='\t')
            header = records.fieldnames
            for line in records:
                line_content = [line[i] for i in header]
                sl = line['patient_ID'].split("_")
                patient = sl[0]
                status = '_'.join(sl[1:])
                n_DNA_bam = patient_files[patient]['normal']['DNA_bam']
                n_RNA_bam = patient_files[patient]['normal']['RNA_bam']
                DNA_bam = patient_files[patient][status]['DNA_bam']
                RNA_bam = patient_files[patient][status]['RNA_bam']
                DNA_tc = patient_files[patient][status]['DNA_tc']
                RNA_tc = patient_files[patient][status]['RNA_tc']
                variant = "_".join([line['patient_ID'],
                                    line['chromosome'],
                                    line['position'],
                                    line['ref_base'],
                                    line['alt_base']])
                # DNA af info
                if (n_DNA_bam != "NA"):
                    try:
                        DNA_af = DNA_af_wn_dict[variant][0].split(',')
                    except KeyError:
                        if (DNA_bam != "NA"):
                            DNA_af = no_coverage
                        else:
                            DNA_af = no_bam
                elif (n_DNA_bam == "NA"):
                    try:
                        DNA_af_tmp = DNA_af_non_dict[variant][0].split(',')
                        DNA_af = nas + DNA_af_tmp + nas + ["na"]
                    except KeyError:
                        if (DNA_bam != "NA"):
                            DNA_af = nas + zeros + nas + ["na"]
                        else:
                            DNA_af = no_bam
                # RNA af info
                if (n_RNA_bam != "NA"):
                    try:
                        RNA_af = RNA_af_wn_dict[variant][0].split(',')
                    except KeyError:
                        RNA_af = no_coverage
                elif (n_RNA_bam == "NA"):
                    try:
                        RNA_af_tmp = RNA_af_non_dict[variant][0].split(',')
                        RNA_af = nas + RNA_af_tmp + nas + ["na"]
                    except KeyError:
                        if (RNA_bam != "NA"):
                            RNA_af = nas + zeros + nas + ["na"]
                        else:
                            RNA_af = no_bam
                final_content = line_content + DNA_af + [DNA_tc] + RNA_af + [RNA_tc]
                writer.writerow(final_content)


def group_patients(patient_files):
    patient_files_wn = dict()
    patient_files_non = dict()
    for patient in patient_files:
        tmp_dict = patient_files[patient]
        # make sure all lower case in input file
        if "normal" in patient_files[patient]:
            patient_files_wn[patient] = tmp_dict
        else:
            patient_files_non[patient] = tmp_dict
    return [patient_files_wn, patient_files_non]


def check_file_permission(patient_files):
    strelka_snv_vcfs = []
    strelka_indel_vcfs = []
    DNA_bams = []
    paired_mpileup_vcfs = []
    mutseq_snv_vcfs = []
    DNA_single_vcfs = []
    RNA_single_vcfs = []
    RNA_bams = []
    for patient in patient_files:
        for status in patient_files[patient]:
            strelka_snv_vcf = patient_files[patient][status]['strelka_snv_vcf']
            strelka_indel_vcf = patient_files[patient][status]['strelka_indel_vcf']
            DNA_single_vcf = patient_files[patient][status]['DNA_single_vcf']
            paired_mpileup_vcf = patient_files[patient][status]['paired_mpileup_vcf']
            mutseq_snv_vcf = patient_files[patient][status]['mutseq_snv_vcf']
            DNA_bam = patient_files[patient][status]['DNA_bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            RNA_single_vcf = patient_files[patient][status]['RNA_single_vcf']
            if (strelka_snv_vcf != "NA"):
                strelka_snv_vcfs.append(strelka_snv_vcf)
            if (strelka_indel_vcf != "NA"):
                strelka_indel_vcfs.append(strelka_indel_vcf)
            if (DNA_single_vcf != "NA"):
                DNA_single_vcfs.append(DNA_single_vcf)
            if (RNA_single_vcf != "NA"):
                RNA_single_vcfs.append(RNA_single_vcf)
            if (paired_mpileup_vcf != "NA"):
                paired_mpileup_vcfs.append(paired_mpileup_vcf)
            if (mutseq_snv_vcf != "NA"):
                mutseq_snv_vcfs.append(mutseq_snv_vcf)
            if (DNA_bam != "NA"):
                DNA_bams.append(DNA_bam)
            if (RNA_bam != "NA"):
                RNA_bams.append(RNA_bam)
    logger.info("Checking strelka snv vcf fileis permissions!")
    check_files(strelka_snv_vcfs)
    logger.info("Checking strelka indel vcf fileis permissions!")
    check_files(strelka_indel_vcfs)
    logger.info("Checking paired_mpileup vcf file permissions!")
    check_files(paired_mpileup_vcfs)
    logger.info("Checking single_mpileup vcf file permissions!")
    check_files(DNA_single_vcfs)
    check_files(RNA_single_vcfs)
    logger.info("Checking bam file permissions!")
    check_files(DNA_bams)
    logger.info("Checking mutseq snv vcf file permissions!")
    check_files(mutseq_snv_vcfs)
    logger.info("Checking RNA_bam file permissions!")
    check_files(RNA_bams)



def __main__():
    parser = argparse.ArgumentParser(description='Summarize variants at both gene and variant level')
    parser.add_argument('-d','--data_type', required=True,
                        help='specify data type as wgs or spc')
    parser.add_argument('-i','--input_file', required=True,
                        help='specify input file, which tells vcf, bam paths')
    args = vars(parser.parse_args())
    logger.info("Gene summary scripts starts at: %s\n" % datetime.datetime.now())
    input_file = args['input_file']
    # impact_types as annotated by snpeff: "HIGH", "MODERATE", "LOW", "MODIFIER"
    impact_types = ["HIGH", "MODERATE"]
    # impact_types = ["HIGH", "MODERATE",  "LOW", "MODIFIER"]
    summary_prefix = '_'.join(impact_types)
    Rscript_path = "/gsc/software/linux-x86_64-centos6/R-3.1.1/bin"
    # r_script = "/home/szong/projects/development/varsummarizer/fisherExactTest.r"
    r_script = os.path.join(os.path.dirname(__file__), 'fisherExactTest.r')
    template_dir = os.path.dirname(__file__)
    # mpileup script template file
    # template_dir = '/home/szong/projects/development/varsummarizer/'
    # snv and indel summaries for patient with or without normal
    hm_snv_sum_wn = '_'.join([summary_prefix, "SNV_summary_with_normal.txt"])
    hm_snv_sum_wn_tmp = ".".join([hm_snv_sum_wn, "tmp"])
    hm_snv_sum_non = '_'.join([summary_prefix, "SNV_summary_no_normal.txt"])
    hm_snv_sum_non_tmp = ".".join([hm_snv_sum_non, "tmp"])
    hm_indel_sum_wn = '_'.join([summary_prefix, "INDEL_summary_with_normal.txt"])
    hm_indel_sum_non = '_'.join([summary_prefix, "INDEL_summary_no_normal.txt"])
    # print help(HEADER)
    logger.info("Removing completion stamps to initialize the pipeline!")
    extension = ['intersect.complete', 'pileup.complete']
    delete_files_by_extension(extension)
    logger.info("Variant input file is: %s" % (input_file))
    logger.info("Generating patient_files dictionary!")
    input_headers = HEADER.INPUT_FILE_HEADER
    patient_files = get_files(input_file, input_headers)
    pprint(patient_files)
    out_dict = group_patients(patient_files)
    patient_files_wn = out_dict[0]
    patient_files_non = out_dict[1]
    
    logger.info("Patients with matched normal are:\n")
    logger.info((patient_files_wn))
    
    logger.info("Patients without matched normal are:\n")
    logger.info((patient_files_non))
    
    logger.info("Checking bam and vcf file permissions!")
    check_file_permission(patient_files)
    
    if args['data_type'] == 'wgs':
        logger.info("Summarizing indels for tumours with matching normals!")
        summarize_indels(patient_files_wn,
                         HEADER.SUM_INDEL_HEADER,
                         hm_indel_sum_wn,
                         impact_types)
    
        logger.info("Summarizing indels for tumours without normals!")
        summarize_indels(patient_files_non,
                         HEADER.SUM_INDEL_HEADER,
                         hm_indel_sum_non,
                         impact_types)
    
        logger.info("Summarizing snvs for tumours with matching normals!")
        (patient_snv_wn,
         pos_files_wn) = summarize_snvs(patient_files_wn,
                                        HEADER.SUM_SNV_HEADER,
                                        hm_snv_sum_wn_tmp,
                                        impact_types)
    
        logger.info("Summarizing snvs for tumours without matching_normal!")
        (patient_snv_non,
         pos_files_non) = summarize_snvs(patient_files_non,
                                         HEADER.SUM_SNV_HEADER,
                                         hm_snv_sum_non_tmp,
                                         impact_types)
        logger.info("Generating mpileup scripts!")
        (pileup_scripts,
         pileup_complete_stamps) = make_pileup_scripts(patient_files, template_dir)

        logger.info("Qsub mpileup scripts!")
        qsub_scripts(pileup_scripts)

        logger.info("Detecting if pileup jobs on cluster finised!")
        detect_cluster_jobs(pileup_complete_stamps)
    
        # logger.info("Deleting intermediate files!")
        # delete_files(pileup_scripts)
        # delete_files (pos_files_non)
        # delete_files (pos_files_wn)
    
        logger.info("Parsing pileup output to get base counts!")
        com_af_files = parse_pileup_output(patient_files,
                                           patient_snv_wn,
                                           patient_snv_non)
    
        logger.info("Adjusting tumour base count by tumour content!")
        for type in com_af_files:
            adjust_allele_frequency(type,
                                    com_af_files[type],
                                    patient_files)

        logger.info("Performing Fisher Exact Test for tumour/normal pairs!")
        run_fisher_exact_test(Rscript_path, r_script)
    
        logger.info("Deleting intermediate files!" )
        extension = ['combined.adjusted', 'combined', '.pileup.log']
        delete_files_by_extension(extension)
    
        logger.info("Figuring out if DNA and/or RNA bams available!\n")
        logger.info("Reading in DNA and RNA af files and making af dictionary!")
        WITH_NORMAL = 0
        NO_NORMAL = 1
        DICT = 0
        COUNTS = 1
        af_out = make_af_dict(patient_files)
        DNA_af_wn_dict = af_out[DICT]["DNA"][WITH_NORMAL]
        DNA_af_non_dict = af_out[DICT]["DNA"][NO_NORMAL]
        RNA_af_wn_dict = af_out[DICT]["RNA"][WITH_NORMAL]
        RNA_af_non_dict = af_out[DICT]["RNA"][NO_NORMAL]
        DNA_AFcounts_afs_wn = af_out[COUNTS]["DNA"]
        RNA_AFcounts_afs_wn = af_out[COUNTS]["RNA"]
        # DNA_af_wn_dict = af_out[0]["DNA"][0]
        # DNA_af_non_dict = af_out[0]["DNA"][1]
        # RNA_af_wn_dict = af_out[0]["RNA"][0]
        # RNA_af_non_dict = af_out[0]["RNA"][1]
        # DNA_AFcounts_afs_wn = af_out[1]["DNA"]
        # RNA_AFcounts_afs_wn = af_out[1]["RNA"]
        
        logger.info("Deleting AFcounts.af files!")
        # delete_files(AFcounts_afs_wn)
        
        logger.info("Combining snv_non summary, af, and annotation results!")
        add_af_anno_for_unpaired(hm_snv_sum_non,
                                 hm_snv_sum_non_tmp,
                                 DNA_af_wn_dict,
                                 RNA_af_wn_dict,
                                 DNA_af_non_dict,
                                 RNA_af_non_dict,
                                 patient_files,
                                 HEADER.SUM_SNV_HEADER_NON)

        logger.info("Combining snv_wn summary, af, and annotation results!")
        add_af_anno_for_paired(hm_snv_sum_wn,
                               hm_snv_sum_wn_tmp,
                               DNA_af_wn_dict,
                               RNA_af_wn_dict,
                               DNA_af_non_dict,
                               RNA_af_non_dict,
                               patient_files,
                               HEADER.SUM_SNV_HEADER_WN)
        
        logger.info("Deleting intermediate files!")
        extension = ['adjusted.fisher', '.vcf', 'AFcounts.af']
        delete_files_by_extension(extension)
        
    # space hold for code development for specifc capture data
    elif args['data_type'] == 'spc':
        print "Summarizing specific capture variants!"

    else:
        print "Unrecognized data type, exiting!"
        sys.exit()
    logger.info("Summarization scripts finished on: %s" %
                datetime.datetime.now())

if __name__ == '__main__':
    __main__()

