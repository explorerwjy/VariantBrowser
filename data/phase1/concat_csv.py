# load multiple csv file with variants and write a 'VCF' file contains these variants
import csv

class Variant():
    def __init__(self,CHR,POS,ID,REF,ALT,QUAL,Filter,INFO,):
        self.CHR = CHR
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT
        self.QUAL = QUAL
        self.FIlter = Filter
        self.make_info(INFO) #self.INFO a dictionary
        self.formkey()
    def formkey(self):
        self.key = '|'.join([self.CHR,self.POS,self.REF,self.ALT])
    # Need to update later
    def make_info(self,INFO):
        res = INFO.split(';')
        self.INFO = {}
        for item in res:
            try:
                key,value = item.split('=')
                self.INFO[key]=value
            except:
                print item
                exit()

    def update_info(self,INFO):
        res = INFO.split(';')
        for item in res:
            try:
                key,value = item.split('=')
                if key in self.INFO:
                    if value not in self.INFO[key].split(','):
                        self.INFO[key] = self.INFO[key] + ',' + value
                else:
                    self.INFO[key] = value
            except:
                print item
                exit()

    def winfo(self):
        
        res = [k+'='+v for k,v in self.INFO.items()]
        return ';'.join(res) 
    def write(self):
        return '\t'.join([self.CHR,self.POS,self.ID,self.REF,self.ALT,self.QUAL,self.FIlter,self.winfo()])+'\n'

def loadcsv(fin,variants):
    fin = open(fin,'rb')
    header = fin.readline().strip().split(',')
    chr_idx = header.index("CHROM")
    pos_idx = header.index("POS")
    try:
        ID_idx = header.index("ID")
    except ValueError:
        ID_idx = -1
    ref_idx = header.index("REF")
    alt_idx = header.index('ALT')
    try:
        qual_idx = header.index('QUAL')
    except ValueError:
        qual_idx = -1
    try:
        filter_idx = header.index('FILTER')
    except ValueError:
        filter_idx = -1
    try:
        info_idx = header.index('INFO')
    except ValueError:
        indo_idx = -1
    reader = csv.reader(fin, delimiter=',')
    count1,count2,count3 = 0,0,0
    for row in reader:
        if qual_idx != -1:
            QUAL = row[qual_idx]
        else:
            QUAL = "30"
        if ID_idx != -1:
            ID = row[ID_idx]
        else:
            ID = "."
        if filter_idx != -1:
            FILTER = row[filter_idx]
        else:
            FILTER = "PASS"
        if info_idx != -1:
            INFO = row[info_idx]
        else:
            INFO = ''
        var = Variant(row[chr_idx],row[pos_idx],ID,row[ref_idx],row[alt_idx],QUAL,FILTER,INFO)
        count1 += 1
        if var.key not in variants:
            variants[var.key] = var
            count2 += 1
        else:
            count3 += 1
            variants[var.key].update_info(row[info_idx])
    print 'Total %d variants.\t%d variants added.\t%d variants updated.'%(count1,count2,count3)
    #return variants

def writeVCF(variants,out_file):
    fout = open(out_file,'wb')
    fout.write('##fileformat=VCFv4.1\n')
    fout.write('##INFO=<ID=Disease,Number=A,Type=String,Description="Diseases this variant might related">\n')
    fout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for var in variants.values():
        fout.write(var.write())
    return

def get_vcf_list(filelist):
    filelist = file(filelist).readlines()
    return [csv.strip() for csv in filelist]

def main():
    CSVs = get_vcf_list("csv.list")
    variants = {}
    for CSV in CSVs:
        print "Reading variants from %s"%CSV
        loadcsv(CSV,variants)

    print 'Total',len(variants.keys()),'variants loaded.'
    
    out_file = 'AllDeNovoVariants.vcf'
    writeVCF(variants,out_file)

if __name__=='__main__':
    main()
