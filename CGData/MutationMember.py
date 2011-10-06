import re
import sys

import CGData
import CGData.GenomicMutation
import CGData.VCF


def get_sample_id(list_of_metadata, sample_id):
    sample_re = re.compile("SAMPLE.*ID="+sample_id)
    for meta in list_of_metadata:
        if sample_re.search(meta):
            id_match = re.search("TCGA-[0-9]{2}-[0-9]{4}-[0-9]{2}[A-Z]-[0-9]{2}[DRTWG]", meta)
            if id_match:
                return id_match.group()
            else:
                print "WARNING: no TCGA id found"
                return ""

SS_DICT = {
           "0": "NA",
           "1": "Germline",
           "2": "Somatic",
           "3": "NA",
           "4": "NA"
           }

class MutationMember(CGData.GenomicMutation.GenomicMutation,
                     CGData.CGGroupMemberSQL):
    def __init__(self):
        CGData.GenomicMutation.GenomicMutation.__init__(self)
        self.vcf = CGData.VCF.VCF()

    def read(self, handle):
        self.vcf.read(handle)
        for line in self.vcf.data:
            if ("NORMAL" not in line.genotype) or ("PRIMARY" not in line.genotype):
                print "expected NORMAl and PRIMARY"
                sys.exit(1)
        
            if line.filter != ["PASS"]:
                pass
                continue
        
            normal = line.genotype["NORMAL"]
            primary = line.genotype["PRIMARY"]
        
            normal_gt = map(int, normal["GT"].split("/"))
            normal_dp = int(normal["DP"])
            normal_ad = map(int, normal["AD"].split(",")) #includes depth of reference
            normal_fa = map(float, normal["FA"].split(",")) #frequency of all ALT combined
            normal_sample_id = get_sample_id(self.vcf.meta, "NORMAL")
            normal_ss = normal["SS"]
            primary_gt = map(int, primary["GT"].split("/"))
            primary_dp = int(primary["DP"])
            primary_ad = map(int, primary["AD"].split(","))
            primary_fa = map(float, primary["FA"].split(","))
            primary_sample_id = get_sample_id(self.vcf.meta, "PRIMARY")
            primary_ss = primary["SS"]
            
            normal_gt.sort()
            primary_gt.sort()
            genotype = [line.ref] + line.alt
            
            for ind, ad in enumerate(normal_ad):
                if ind == 0: #ignore reference
                    continue
                if ad == 0: #ignore alts with no reads
                    continue
                mutation_line = []
                mutation_line.append(normal_sample_id[:-8])
                mutation_line.append(line.chrom)
                mutation_line.append(line.pos)
                mutation_line.append(line.pos+len(line.ref)-1)
                mutation_line.append(line.ref)
                mutation_line.append(genotype[ind])
                mutation_line.append(normal_dp)
                mutation_line.append(1.0*ad/normal_dp)
                mutation_line.append("NA")
                mutation_line.append("NA")
                mutation_line.append(SS_DICT[normal_ss])
                mutation_line.append(normal_sample_id)
                mutation_line.append("GENE_NAME")
                mutation_line.append("VARIANT_CLASSIFICATION")
                self.data.append(CGData.GenomicMutation.GenomicMutationLine(*mutation_line))
        
            for ind, ad in enumerate(primary_ad):
                if ind == 0:
                    continue
                if ad == 0:
                    continue
                mutation_line = []
                mutation_line.append(primary_sample_id[:-8])
                mutation_line.append(line.chrom)
                mutation_line.append(line.pos)
                mutation_line.append(line.pos+len(line.ref)-1)
                mutation_line.append(line.ref)
                mutation_line.append(genotype[ind])
                mutation_line.append(primary_dp)
                mutation_line.append(1.0*ad/primary_dp)
                mutation_line.append("NA")
                mutation_line.append("NA")
                mutation_line.append(SS_DICT[primary_ss])
                mutation_line.append(primary_sample_id)
                mutation_line.append("GENE_NAME")
                mutation_line.append("VARIANT_CLASSIFICATION")
                self.data.append(CGData.GenomicMutation.GenomicMutationLine(*mutation_line))

    def gen_sql_head(self):
        yield "drop table if exists mutation_%s;\n" % self.attrs['group']
        yield "CREATE TABLE mutation_%s (\n" % self.attrs['group']
        yield "'id' INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,\n"
        yield "'samplename' VARCHAR(255) DEFAULT NULL,\n"
        yield "'chrom' VARCHAR(255) DEFAULT NULL,\n"
        yield "'startposition' INT(10) UNSIGNED DEFAULT NULL,\n"
        yield "'endposition' INT(10) UNSIGNED DEFAULT NULL,\n"
        yield "'ref' VARCHAR(255) DEFAULT NULL,\n"
        yield "'alt' VARCHAR(255) DEFAULT NULL,\n"
        yield "'coverage' INT(10) UNSIGNED DEFAULT NULL,\n"
        yield "'freq' DOUBLE DEFAULT NULL,\n"
        yield "'rnacoverage' INT(10) UNSIGNED DEFAULT NULL,\n"
        yield "'rnafreq' DOUBLE DEFAULT NULL,\n"
        yield "'gstype' ENUM( 'Germline', 'Somatic', 'NA' ) default NULL,\n"
        yield "'patientname' VARCHAR(255) DEFAULT NULL,\n"
        yield "'hugosymbol' VARCHAR(255) DEFAULT NULL,\n"
        yield ("'variant_classification' "
               "ENUM('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', "
               "'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', "
               "'Silent', 'Splice_Site', 'Translation_Start_Site', "
               "'Nonstop_Mutation', '3UTR', '3Flank', '5UTR', '5Flank', "
               "'IGR1' , 'Intron', 'RNA', 'Targeted_Region', 'Indel', "
               "'De_novo_Start_InFrame', 'De_novo_Start_OutOfFrame') "
               "DEFAULT NULL\n")
        yield ");\n"

    
    def gen_sql_lines(self, id_table):
        for line in self.data:
            yield "INSERT INTO mutation%s VALUES (%d, %s);\n"% (self.attrs['group'], id_table.get('mutation_id', line), line.sqlline())

    
    def gen_sql_tail(self):
        return []
