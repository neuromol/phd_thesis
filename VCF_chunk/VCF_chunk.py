import os
import re
import pandas as pd
import sys
import argparse

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Phase and minimac3 - Josh Atkins c3114203 at uon.edu.au ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input vcf file',required='True', default='plink')
    parser.add_argument('-s', '--chunk_size', help='Chunk size',required='True', default='6000000')
    parser.add_argument('-f', '--over_lap', help='Flanking region ', default='500000')
    parser.add_argument('-o', '--output', help='output file DEFAULT: Results_label/', default='Results')
    parser.add_argument('-c', '--chr', help='Chromosome', default='1')
    parser.add_argument('-r', '--ref', help='Reference Panel directory' , default='/home/username/imputation/1000g_v5_impute/')
    parser.add_argument('-w', '--working', help='Working directory' , default='/home/username/imputation/dataset/')


    results = parser.parse_args(args)
    return (results.input , results.chunk_size , results.over_lap , results.output , results.chr , results.ref , results.working )


def count_comments(filename):
    comments = 0
    with open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                comments += 1
            else:
                break
    return comments

def header_line(f):
    t =open(f , 'r')
    for line in t:
        if line.startswith('#CHRO'):
            return(line)

def closest(df, col, val, direction):
    n = len(df[df[col] <= val])
    if(direction < 0):
        n -= 1
    if(n < 0 or n >= len(df)):
        print('err - value outside range')
        return None
    return df.ix[n, col]


#def range_file( title , startbp , endbp , Num_SNPs):
#    output = (output_vcf + "\t range: " + str(startbp) + "-" + str(endbp) + "\t  Num_SNPs= " + str(Num_SNPs) + "\n" )
#    f.write(output)
#


def create_header(file, header_count):
    m = open(file, 'r')
    temp = open("tempheader.txt", "w")
    for i in range(header_count):
        line=m.next().strip()
        temp.write(line + "\n")
    m.close()




#
def write_script(ref , chr , base_label , start_overlap , end_overlap , overlap , chunk , file, working):
    output_vcf = (base_label + "." +chunk + ".vcf" )
    chunk_sh = ((chunk) + ".sh")
    t = open(chunk_sh , 'w')



    t.write("""#PBS -l select=1:ncpus=8:mem=31G,walltime=12:00:00 \n""")
    t.write("\n")
    #t.write("source /etc/profile.d/modules.sh \n")
    #t.write("module load bcftools/1.4.1 \n")
    #t.write("module load samtools/1.4.1   \n")
    #t.write("module load eagle/2.3.2   \n")
    t.write("cd " + working +"\n")
    t.write("bgzip -c " + output_vcf + " > " + output_vcf + ".gz \n")
    t.write("tabix -p vcf " + output_vcf + ".gz \n")
    t.write("eagle --vcfRef "+ ref + "chr" + str(chr) + ".vcf.gz --vcfTarget " + working  + output_vcf + ".gz --allowRefAltSwap --bpStart " + str(start_overlap) + " --bpEnd " + str(end_overlap) + " --bpFlanking " + str
(overlap) + " --outPrefix "+ chunk +".phased --geneticMapFile=/Eagle_v2.3.2/tables/genetic_map_hg19_withX.txt.gz --numThreads 8 \n" )
    t.write("Minimac3-omp --refHaps "+ ref + "chr" + str(chr) +".vcf.gz --haps "+ chunk +".phased.vcf.gz --prefix "+ chunk +".imputed --chr "+ chr +" --start " + str(start_overlap) + " --end " + str(end_overlap) + " --
window " + str(overlap) + " --log --cpus 8 --noPhoneHome \n")
    t.write("plinknew --vcf "+ chunk +".imputed.dose.vcf.gz --maf 0.01 --make-bed --out " + chunk )


def main(file , chunk_size , over_lap , output , chr , ref , working , x):
    #os.chdir("chr"+chr)
    chunk_size = int(chunk_size)
    over_lap = int(over_lap)
    base_label = str(input[:-4])
    reformat = (header_line(file).split()  )
    comments = count_comments(file)
    start_line= (comments + 1)
    main_part = pd.read_csv(file, skiprows=comments, header=None ,sep="\t" )
    main_part.columns = reformat
    first_SNP = int(main_part['POS'].iloc[0])
    last_SNP = int(main_part['POS'].iloc[-1:])
    startbp = int(first_SNP)
    pluspos = int(first_SNP + chunk_size)
    endbp = int(startbp + chunk_size)
    f = open('ranges.txt','w')
    counter = 1
    startbp = first_SNP

    while startbp <= last_SNP :

        start_overlap = (int(startbp) - int(over_lap))
        end_overlap = (int(endbp) + int(over_lap))
        section = main_part[(main_part['POS'] <= end_overlap ) & (main_part['POS'] >= start_overlap )]
        Num_SNPs = len(section.axes[0])
        output_vcf = (base_label + ".chunk_" + str(counter) + ".vcf" )
        chunk = ("chunk_"  + str(counter) )
        Num_SNPs=len(section.axes[0])
        write_script( ref , chr , base_label , start_overlap , end_overlap , over_lap , chunk , file , working)
        #range_file()
        v = open(output_vcf, 'w')
        temp = open("tempheader.txt", "r")
        for line in temp:

            v.write(line)

        temp.close
        v.close
        section.to_csv(output_vcf, header=True , mode='a' , sep="\t" , index = False)
        print (Num_SNPs)
        counter += 1
        startbp = startbp + chunk_size
        endbp = endbp + chunk_size

if __name__ == '__main__':

    input , chunk_size , over_lap , output , chr , ref , working = check_arg(sys.argv[1:])
    x = count_comments(input)
    create_header(input,x)
    main(input , chunk_size , over_lap , output , chr , ref , working , x)
