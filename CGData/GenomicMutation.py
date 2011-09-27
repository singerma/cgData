import CGData

class GenomicMutationLine:
    def __init__(self, samplename, chrom, start, end, ref, alt, coverage,
                 freq, rna_coverage, rna_freq, gstype, patient, hugo,
                 var_class):
        self.samplename = samplename
        self.chrom = chrom
        self.start = start
        self.end = end
        self.ref = ref
        self.alt = alt
        self.coverage = coverage
        self.freq = freq
        self.rna_coverage = rna_coverage
        self.rna_freq = rna_freq
        self.gstype = gstype
        self.patient = patient
        self.hugo = hugo
        self.var_class = var_class
    
    def sqlline(self):
        return ",".join(map(lambda x: """'%s'""" % str(x),
                            [self.samplename, self.chrom, self.start,
                             self.end, self.ref, self.alt, self.coverage,
                             self.freq, self.rna_coverage, self.rna_freq,
                             self.gstype, self.patient, self.hugo,
                             self.var_class]))

class GenomicMutation(CGData.CGDataSetObject, CGData.CGSQLObject):
    def __init__(self):
        CGData.CGDataSetObject.__init__(self)
        self.data = []
    
    def read(self, handle):
        for line in handle:
            splitline = line.strip().split("\t")
            self.data.append(GenomicMutationLine(*splitline))
    
    def __iter__(self):
        for line in self.data:
            yield line
    
    def init_schema(self):
        pass
    
    def build_ids(self, id_allocator):
        if self.light_mode:
            self.load()
        
        for line in self.data:
            id_allocator.alloc('mutation_id', line)
    
    def gen_sql(self, id_table):
        yield "drop table if exists mutation_%s;\n" % self.attrs['name']
        yield "CREATE TABLE mutation_%s (\n" % self.attrs['name']
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
        yield "'gstype' ENUM('Germline', 'Somatic', 'NA') default NULL,\n"
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
        
        for line in self:
            yield "INSERT INTO mutation%s VALUES (%d, %s);\n"% (self.attrs['name'], id_table.get('mutation_id', line), line.sqlline())

    
