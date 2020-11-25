import sys
import argparse
import pathogenprofiler as pp
from collections import Counter
import csv


acgt = ["A","C","G","T"]

def get_fst_cat(fst):
    if fst==1:
        return "1"
    elif fst<1 and fst >=0.99:
        return "0.99"
    elif fst<0.99 and fst>=0.95:
        return ">=0.95&<0.99"

def main(args):
    vcf = pp.vcf(args.vcf)
    if args.samples_blacklist:
        samples_blacklist = [x.strip() for x in open(args.samples_blacklist).readlines()]
    else:
        samples_blacklist = []

    in_samples = [x.strip().replace(" ","_") for x in open(args.in_samples).readlines()]
    in_samples = [s for s in in_samples if s in vcf.samples and s not in samples_blacklist]
    out_samples = [s for s in vcf.samples if s not in in_samples and s not in samples_blacklist]

    tmp_in_file = pp.get_random_file()
    with open(tmp_in_file,"w") as O:
        O.write("\n".join(in_samples)+"\n")

    tmp_out_file = pp.get_random_file()
    with open(tmp_out_file,"w") as O:
        O.write("\n".join(out_samples)+"\n")

    if pp.nofile(f"{args.prefix}.weir.fst") or args.redo:
        pp.run_cmd(f"vcftools --gzvcf {args.vcf} --weir-fst-pop {args.in_samples} --weir-fst-pop {tmp_out_file} --out {args.prefix}")


    high_fst_snps = []
    tmp_bed_file = pp.get_random_file(extension=".bed")
    with open(tmp_bed_file,"w") as O:
        for line in open(args.prefix+".weir.fst"):
            row = line.strip().split()
            if row[2]=="WEIR_AND_COCKERHAM_FST": continue
            if float(row[2]) >= args.fst_cutoff:
                high_fst_snps.append( {"pos": int(row[1]), "fst": float(row[2])} )
                O.write("Chromosome\t%s\t%s\n" % (int(row[1])-1,row[1]))

    vcf_alleles = {}
    for line in pp.cmd_out("bcftools view %s -R %s | bcftools csq -p a -f %s -g %s | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%BCSQ[\\t%%IUPACGT]\\n' " % (args.vcf,tmp_bed_file,args.ref,args.gff)):
        row = line.strip().split()
        new_row = {"pos":int(row[0]),"ref":row[1], "alt": row[2], "bcsq":row[3]}
        for i in range(4,len(row)):
            new_row[vcf.samples[i-4]] = row[i]
        vcf_alleles[int(row[0])] = new_row


    coll_snps = {}
    for line in open(args.coll_snps):
        row = line.strip().split()
        coll_snps[int(row[1])] = row[0]


    results = []

    for snp in high_fst_snps:

        allele = Counter([vcf_alleles[snp["pos"]][s] for s in in_samples if vcf_alleles[snp["pos"]][s] in acgt]).most_common(1)[0][0]

        result_row = {
            "lineage": args.prefix,
            "pos": snp["pos"],
            "In_coll": 1 if snp["pos"] in coll_snps else 0,
            "Fst": snp["fst"],
            "Fst_cat": get_fst_cat(snp["fst"]),
            "ref": vcf_alleles[snp["pos"]]["ref"],
            "alt": vcf_alleles[snp["pos"]]["alt"],
            "allele": allele
        }

        if "," in vcf_alleles[snp["pos"]]["alt"]: continue ### Remove multi alleleics

        csq_str = vcf_alleles[snp["pos"]]["bcsq"]
        csq_list = csq_str.split("|")
        print(csq_str)
        if csq_str[0]=="@":
            continue
        elif csq_str==".":
            result_row["change_type"] = "intergenic"
            result_row["gene"] = "intergenic"
            result_row["ensembl_ID"] = "intergenic"
            result_row["biotype"] = "intergenic"
            result_row["coding_strand"] = "NA"
            result_row["aa_pos"] = "NA"
            result_row["change"] = "NA"
        elif csq_list[0]=="start_lost":
            result_row["change_type"] = csq_list[0]
            result_row["gene"] = csq_list[1]
            result_row["ensembl_ID"] = csq_list[2]
            result_row["biotype"] = csq_list[3]
            result_row["coding_strand"] = csq_list[4]
            result_row["aa_pos"] = "NA"
            result_row["change"] = "NA"
        elif csq_list[0]=="non_coding":
            result_row["change_type"] = csq_list[0]
            result_row["gene"] = csq_list[1]
            result_row["ensembl_ID"] = "NA"
            result_row["biotype"] = csq_list[3]
            result_row["coding_strand"] = "NA"
            result_row["aa_pos"] = "NA"
            result_row["change"] = "NA"
        else:
            result_row["change_type"] = csq_list[0]
            result_row["gene"] = csq_list[1]
            result_row["ensembl_ID"] = csq_list[2]
            result_row["biotype"] = csq_list[3]
            result_row["coding_strand"] = csq_list[4]
            result_row["aa_pos"] = csq_list[5]
            result_row["change"] = csq_list[6]


        results.append(result_row)

    with open(args.prefix+".results.txt","w") as O:
        writer = csv.DictWriter(O,delimiter=args.sep,fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    pp.rm_files([tmp_out_file,tmp_bed_file])

parser = argparse.ArgumentParser(description='Run fst pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--prefix',help='Prefix for output files',required=True)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--ref',help='Reference fasta file',required=True)
parser.add_argument('--gff',help='Reference gff file',required=True)
parser.add_argument('--in-samples',help='Samples belonging to lienage',required=True)
parser.add_argument('--coll-snps',help='List of Coll et al SNPs',required=True)
parser.add_argument('--fst-cutoff',type=float,default=0.95,help='Cutoff for retaining SNPs')
parser.add_argument('--sep',type=str,default="\t",help='Column seperator for output file')
parser.add_argument('--samples-blacklist',type=str,help='Sample to ignore from analysis')
parser.add_argument('--redo',action="store_true",help='Overwrite existing files')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
