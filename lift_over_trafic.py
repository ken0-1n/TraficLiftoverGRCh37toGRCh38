from __future__ import print_function 
import os, sys
import subprocess

def get_info_val(infos, key):
    ret = ""
    l_info = infos.split(';')
    for info in l_info:
        if info.startswith(key):
            ret = info.split("=")[1]
            break
    return ret

def replace_info_val(infos, key, new_val):
    val = ""
    l_info = infos.split(';')
    for info in l_info:
        if info.startswith(key):
            val = info.split("=")[1]
            return(infos.replace(key+val,key+new_val))
    return infos

def make_hash_chr(in_bed, h):
    with open(in_bed, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            key = F[3]
            val = F[0].replace("chr","") +"\t"+ F[2]
            h[key] = val
    return h
    
def make_hash_bkpb(in_bed, h, flg):
    with open(in_bed, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            F = line.split('\t')
            key = F[3]
            val = F[2] if flg else "Liftover_Fail"
            h[key] = val
    return h
    
def make_hash_src_tdc(in_bed, h, flg):
    with open(in_bed, 'r') as hin:
        for line in hin:
            if line.startswith("#"): continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrom = F[0].replace("chr","")
            start = F[1]
            end = F[2]
            key = F[3]
            reverse_flg = F[4]
            if int(reverse_flg) == 1:
                start = F[2]
                end = F[1]
            val = chrom +"_"+ start +"_"+ end if flg else "Liftover_Fail"
            h[key] = val
    return h
    
def check_position_format(l_src_val):
    reverse_flg = "0"
    start = l_src_val[1]
    end = l_src_val[2]
    if (l_src_val[2] < l_src_val[1]):
        start = l_src_val[2]
        end = l_src_val[1]
        reverse_flg = "1"
    return start, end, reverse_flg


def make_liftover_input(invcf, out_pref):
    
    hout_chr  = open(out_pref+ ".chr.bed", 'w')
    hout_bkpb = open(out_pref+ ".bkpb.bed", 'w')
    hout_src  = open(out_pref+ ".src.bed", 'w')
    hout_tdc  = open(out_pref+ ".tdc.bed", 'w')
    
    with open(in_vcf, 'r') as hin:
        for line in hin:
            if line.startswith("#"):
                # header = line.rstrip('\n')
                # print(header, file=hOUT)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrom = F[0]
            pos = F[1]
            info = F[7]
            key = chrom+":"+pos
            
            # lift_over chromosome position
            print("chr"+chrom+"\t"+str(int(pos)-1)+"\t"+pos+"\t"+key, file=hout_chr)
            
            # lift_over bkpb
            bkpb_val = get_info_val(info, "BKPB=")
            if bkpb_val != "":
                print("chr"+chrom+"\t"+str(int(bkpb_val)-1)+"\t"+bkpb_val+"\t"+key, file=hout_bkpb)
    	
    	    # lift_over src
            src_val = get_info_val(info, "SRC=")
            if src_val != "":
                l_src_val = src_val.split('_')
                start, end, reverse_flg = check_position_format(l_src_val)
                print("chr"+l_src_val[0]+"\t"+start+"\t"+end+"\t"+key+"\t"+reverse_flg, file=hout_src)
    
            # lift_over tdc
            tdc_val = get_info_val(info, "TDC=")
            if tdc_val != "":
                l_tdc_val = tdc_val.split('_')
                start, end, reverse_flg = check_position_format(l_tdc_val)
                print("chr"+l_tdc_val[0]+"\t"+start+"\t"+end+"\t"+key+"\t"+reverse_flg, file=hout_tdc)
    	
    hout_chr.close()
    hout_bkpb.close()
    hout_src.close()
    hout_tdc.close()


def make_hash_from_liftover_output():
    h_chr = {}
    h_chr = make_hash_chr(out_pref+ ".chr_GRCh38.bed", h_chr)
    
    h_bkpb = {}
    h_bkpb = make_hash_bkpb(out_pref+ ".bkpb_GRCh38.bed", h_bkpb, True)
    h_bkpb = make_hash_bkpb(out_pref+ ".bkpb_unmapped.bed", h_bkpb, False)
    
    h_src = {}
    h_src = make_hash_src_tdc(out_pref+ ".src_GRCh38.bed", h_src, True)
    h_src = make_hash_src_tdc(out_pref+ ".src_unmapped.bed", h_src, False)
    
    h_tdc = {}
    h_tdc = make_hash_src_tdc(out_pref+ ".tdc_GRCh38.bed", h_tdc, True)
    h_tdc = make_hash_src_tdc(out_pref+ ".tdc_unmapped.bed", h_tdc, False)

    return h_chr, h_bkpb, h_src, h_tdc
    

def write_result(out_vcf, h_chr, h_bkpb, h_src, h_tdc):
    
    hout  = open(out_vcf, 'w')

    with open(in_vcf, 'r') as hin:
        for line in hin:
            if line.startswith("#"):
                header = line.rstrip('\n')
                if header.startswith("##contig="):continue
                print(header.replace("##reference=hg19","##reference=hg38"), file=hout)
                continue
            line = line.rstrip('\n')
            F = line.split('\t')
            chrom = F[0]
            pos = F[1]
            infos = F[7]
            key = chrom+":"+pos
            
            # liftover chromosome position
            chrom, pos = h_chr[key].split("\t")
            
            # liftover bkpb
            if key in h_bkpb:
                infos = replace_info_val(infos, "BKPB=", h_bkpb[key])
    
            # liftover src
            if key in h_src:
                infos = replace_info_val(infos, "SRC=", h_src[key])
    
            # liftover tdc
            if key in h_tdc:
                infos = replace_info_val(infos, "TDC=", h_tdc[key])
    
            print("\t".join([chrom, pos]+ F[2:7] + [infos, F[8], F[9]]), file=hout)



### Main ###
in_vcf = sys.argv[1]
out_vcf = sys.argv[2]
map_chain = sys.argv[3]

out_pref, ext = os.path.splitext(out_vcf)

make_liftover_input(in_vcf, out_pref)

# liftover chromosome position
subprocess.check_call(["liftOver", out_pref+ ".chr.bed", map_chain, out_pref+ ".chr_GRCh38.bed", out_pref+ ".chr_unmapped.bed"])
        
# liftover bkpb
subprocess.check_call(["liftOver", out_pref+ ".bkpb.bed", map_chain, out_pref+ ".bkpb_GRCh38.bed", out_pref+ ".bkpb_unmapped.bed"])

# liftover src 
subprocess.check_call(["liftOver", out_pref+ ".src.bed", map_chain, out_pref+ ".src_GRCh38.bed", out_pref+ ".src_unmapped.bed"])

# liftover tdc 
subprocess.check_call(["liftOver", out_pref+ ".tdc.bed", map_chain, out_pref+ ".tdc_GRCh38.bed", out_pref+ ".tdc_unmapped.bed"])
        
h_chr, h_bkpb, h_src, h_tdc = make_hash_from_liftover_output()

write_result(out_vcf, h_chr, h_bkpb, h_src, h_tdc)



