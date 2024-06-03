from Bio import Entrez
from Bio import SeqIO
from Bio import SeqUtils
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import Align
from os import listdir
from tqdm import tqdm
from time import sleep
import time
import datetime
import os
import subprocess
import re
import random
import pandas as pd
import numpy as np
import primer3
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Primers_dict={}
Primers_dict2={}
Primers_dict_comb={}
Primers_dict_in_silico={}
Blastn_result_list=[]
Blastn_result_list_alig=[]

class RHDV_primer_lit:

    def _init_(self):
        pass

    def parse_primers_RHDV_test(self,database_file):
            print ("Global module to read primers for RHDV")
            print ("READING FOR: RHDV_primers.csv")
            RHDV_primers_data=open("CONSERVATION_SCORES\Test.csv","r")
            count_primers=0
            for line in RHDV_primers_data:
                #print (line)
                count_primers+=1
                primer_data=line.split(";")
                primer_id=primer_data[0]
                primer_sequence=primer_data[1]
                if primer_sequence.startswith("N") or primer_sequence.startswith("K") or primer_sequence.startswith("S") or primer_sequence.startswith("W") or primer_sequence.startswith("M") or primer_sequence.startswith("R") or primer_sequence.startswith("Y"):
                    primer_sequence=primer_sequence[1:]
                if primer_sequence.endswith("N") or primer_sequence.endswith("K") or primer_sequence.endswith("S") or primer_sequence.endswith("W") or primer_sequence.endswith("M") or primer_sequence.endswith("R") or primer_sequence.endswith("Y"):
                    primer_sequence=primer_sequence[:-1]
                primer_reference=primer_data[2]
                primer_reference2=primer_reference.replace(",","|")
                primer_reference_html=primer_reference2.split("|")
                primer_reference_table=""
                for ref in primer_reference_html:
                    ref="<a href="+ref+">"+ref+"</a><br>"
                    primer_reference_table +=ref
                primer_length=len(primer_sequence)
                #primer_length=primer_data[3]
                #primer_length2=primer_length.replace("\n","")
                Primers_dict[primer_id]=[primer_id,primer_sequence,primer_reference_table,primer_length] 
                print(primer_sequence) 
            #for primers in Primers_dict:
                #print (Primers_dict[primers])      
            print (Primers_dict)
            pass

    def parse_primers_RHDV(self,database_file):
            print ("Global module to read primers for RHDV")
            print ("READING FOR: RHDV_primers.csv")
            RHDV_primers_data=open("CONSERVATION_SCORES\Virus\RHDV\/finalprimersrhdv03072023.csv","r")
            count_primers=0
            for line in RHDV_primers_data:
                #print (line)
                count_primers+=1
                primer_data=line.split(";")
                primer_id=primer_data[0]
                primer_id2=primer_id.replace("ï»¿","")
                primer_sequence=primer_data[1]
                if primer_sequence.startswith("N") or primer_sequence.startswith("K") or primer_sequence.startswith("S") or primer_sequence.startswith("W") or primer_sequence.startswith("M") or primer_sequence.startswith("R") or primer_sequence.startswith("Y"):
                    primer_sequence=primer_sequence[1:]
                if primer_sequence.endswith("N") or primer_sequence.endswith("K") or primer_sequence.endswith("S") or primer_sequence.endswith("W") or primer_sequence.endswith("M") or primer_sequence.endswith("R") or primer_sequence.endswith("Y"):
                    primer_sequence=primer_sequence[:-1]
                primer_reference=primer_data[2]
                primer_reference2=primer_reference.replace(",","|")
                primer_reference3=primer_reference2.replace("'","")
                primer_reference4=primer_reference3.replace("[","")
                primer_reference5=primer_reference4.replace("]","")
                primer_reference6=primer_reference5.replace(" ","")
                primer_reference_html=primer_reference6.split("|")
                primer_reference_table=""
                for ref in primer_reference_html:
                    ref="<a href="+ref+">"+ref+"</a><br>"
                    primer_reference_table +=ref
                primer_length=len(primer_sequence)
                #primer_length=primer_data[3]
                #primer_length2=primer_length.replace("\n","")
                Primers_dict[primer_id]=[primer_id2, primer_sequence, primer_reference_table, primer_length] 
            #print(primer)
            #for primers in Primers_dict:
                #print (Primers_dict[primers])      
            #print (Primers_dict)
            pass
    
    def parse_primers_final(self,database_file):
            print ("Global module to read primers for RHDV")
            print ("READING FOR: RHDV_primers.csv")
            RHDV_final_data=open("CONSERVATION_SCORES\/final_AROLIT.csv","r")
            count_primers=0
            for line in RHDV_final_data:
                #print (line)
                count_primers+=1
                primer_data=line.split(";")
                primer_id=primer_data[0]
                primer_sequence=primer_data[1]
                primer_reference=primer_data[2]
                primer_length=primer_data[3]
                primer_ref_start=primer_data[4]
                primer_ref_end=primer_data[5]
                #primer_ref_seq=primer_data[9]
                gc_content=primer_data[6]
                melting_temp=primer_data[7]
                primer_type=primer_data[8]
                scores=primer_data[9]
                scores2=scores.replace("\n","")
                scores3=scores2.replace(" ","")
                Primers_dict2[primer_id]=[primer_id,primer_sequence,primer_reference,primer_length,primer_ref_start,primer_ref_end,gc_content,melting_temp,primer_type,scores3]
            #for primers in Primers_dict:
                #print (Primers_dict[primers])      
            #print (Primers_dict2)
            pass
    
    def parse_primers_in_silico(self,database_file,alignment_length):
        print ("Global module to read primers for RHDV")
        print ("READING FOR: RHDV_primers.csv")
        RHDV_in_silico_data=open("CONSERVATION_SCORES\Virus\RHDV\CSVsInSilicoPrimers1704.csv","r")
        count_primers=0
        for line in RHDV_in_silico_data:
            #print (line)
            count_primers+=1
            primer_data=line.split(";")
            primer_id=primer_data[0]
            primer_sequence=primer_data[1]
            primer_sequence2=primer_sequence.replace("U","T")
            primer_reference=primer_data[2]
            primer_length=primer_data[3]
            primer_start=primer_data[4]
            primer_end=primer_data[5]
            primer_end2=primer_end.replace("\n","")
            if primer_reference=="in silico primer reverse":
                primer_start=alignment_length-int(primer_start)
                primer_end2=alignment_length-int(primer_end2)
            Primers_dict[primer_id]=[primer_id, primer_sequence2, primer_reference, primer_length, primer_start, primer_end2]

    def create_blastn_query(self,path, seq_name, sequence):
        root=os.getcwd()
        os.chdir(root+'/'+path)
        file_to_save="query.fasta"
        query_fasta=open(file_to_save,'w')
        query_fasta.write(">"+seq_name+"\n")
        query_fasta.write(sequence)
        os.chdir(root)

    def create_blast_db(self,path,fasta_file):
        root=os.getcwd()
        os.chdir(root+'/'+path)
        command_line="makeblastdb -in " +fasta_file+ " -dbtype nucl -logfile log1.txt"
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print (blast_results, err)
        status = child_process.wait()
        print (status)
        os.chdir(root)

    def create_fasta_ref_from_alig(self,path,fasta_file):
        root=os.getcwd()
        os.chdir(root+'/'+path)
        ref_fasta=open("REFblastdbalig.fasta", "w")
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_name=record.id
            sequence=record.seq
            if seq_name == "NC_001543":
                ref_fasta.write(">"+seq_name+"\n")
                ref_fasta.write(str(sequence.upper())+"\n")
                #print(record.id)
                time.sleep(3)

    def create_blast_db_alig(self,path,fasta_file):
        root=os.getcwd()
        os.chdir(root+'/'+path)
        command_line="makeblastdb -in " +"REFblastdbalig.fasta"+ " -dbtype nucl -logfile log1.txt"
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print (blast_results, err)
        status = child_process.wait()
        print (status)
        os.chdir(root)

    def run_blastn_ref(self,path,db_file,query_file,primer_id):
        root=os.getcwd()
        os.chdir(root+'/'+path)
        command_line="blastn -gapopen 2 -word_size 18 -query "+ query_file +" -db "+ db_file +" -task blastn-short -out RESULTS/RHDV_"+ primer_id+".out -logfile RESULTS/log2_"+primer_id+".txt"
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print (blast_results, err)
        status = child_process.wait()
        print (status)
        os.chdir(root)

    def run_blastn_alig(self,path,db_file,query_file,primer_id):
        root=os.getcwd()
        os.chdir(root+'/'+path)
        command_line="blastn -gapopen 2 -word_size 18 -query "+ query_file +" -db "+ db_file +" -task blastn-short -out RESULTS2/RHDV_"+ primer_id+".out -logfile RESULTS2/log2_"+primer_id+".txt"
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print (blast_results, err)
        status = child_process.wait()
        print (status)
        os.chdir(root)

    def read_blastn_output_ref(self, input_file, primer_id):
        blastn_results=open(input_file,"r")
        Blastn_result_list=[]
        for line in blastn_results:
            Blastn_result_list.append(line)
            #print (line)
        #for data in Blastn_result_list_FC:
        position_primer=Blastn_result_list[31]
        position_alignment=Blastn_result_list[33]
        #print (position_primer)
        primer_blast_data=position_primer.split("  ")
        ref_blast_data=position_alignment.split("  ")
        #print (primer_blast_data)
        #print (len (primer_blast_data))
        if len(primer_blast_data)==5:
            del primer_blast_data[2]
        elif len(primer_blast_data)==3:
            primer_blast_data=primer_blast_data.split(" ")
        elif len(primer_blast_data)==4:
            pass
        if len(ref_blast_data)==5:
            del ref_blast_data[2]
        elif len(ref_blast_data)==3:
            ref_blast_data=ref_blast_data.split(" ")
        elif len(ref_blast_data)==4:
            pass
            ref_primer_start=primer_blast_data[1]
            ref_primer_start2=ref_primer_start.replace(" ","")
            ref_primer_end=primer_blast_data[3]
            ref_primer_end2=ref_primer_end.replace("\n","")
            ref_primer_end3=ref_primer_end2.replace(" ","")
            ref_primer_seq=primer_blast_data[2]
            ref_primer_seq2=ref_primer_seq.replace(" ","")
            ref_start=ref_blast_data[1]
            ref_start2=ref_start.replace(" ","")
            ref_end=ref_blast_data[3]
            ref_end2=ref_end.replace("\n","")
            ref_end3=ref_end2.replace(" ","")
            ref_seq=ref_blast_data[2]
            ref_seq2=ref_seq.replace(" ","")
            #print(ref_primer_start, ref_primer_end2, ref_primer_seq)
            #print (primer_id)
            #if primer_id in Primers_dict:
                #print (primer_id)
                #Primers_dict[primer_id].append(str(Blast_dict[primer_id][1]))
                #print (Primers_dict[primer_id])
            Primers_dict[primer_id].append(ref_primer_start2)
            Primers_dict[primer_id].append(ref_primer_end3)
            Primers_dict[primer_id].append(ref_primer_seq2)
            Primers_dict[primer_id].append(ref_start2)
            Primers_dict[primer_id].append(ref_end3)
            Primers_dict[primer_id].append(ref_seq2)
        #print (Primers_dict[primer_id])
        #print (Primers_dict[str(1)])
        Blastn_result_list=[]

    def read_blastn_output_alig(self, input_file, primer_id):
        blastn_results_alig=open(input_file,"r")
        Blastn_result_list_alig=[]
        for line in blastn_results_alig:
            Blastn_result_list_alig.append(line)
            #print (line)
        #for data in Blastn_result_list_FC:
        position_primer_alig=Blastn_result_list_alig[31]
        position_alignment_alig=Blastn_result_list_alig[33]
        #print (position_alignment_alig)
        primer_blast_data_alig=position_primer_alig.split("  ")
        alignment_blast_data_alig=position_alignment_alig.split("  ")
        #print (len (primer_blast_data_alig))
        #print (str(alignment_blast_data_alig))
        if len(primer_blast_data_alig)==5:
            del primer_blast_data_alig[2]
        elif len(primer_blast_data_alig)==3:
            primer_blast_data_alig=position_primer_alig.split(" ")
        elif len(primer_blast_data_alig)==4:
            pass
        if len (alignment_blast_data_alig)==5:
            del alignment_blast_data_alig[2]
        elif len(alignment_blast_data_alig)==3:
            alignment_blast_data_alig=position_alignment_alig.split(" ")
        elif len(alignment_blast_data_alig)==4:
            pass
            alig_primer_start=primer_blast_data_alig[1]
            alig_primer_start2=alig_primer_start.replace(" ","")
            alig_primer_end=primer_blast_data_alig[3]
            alig_primer_end2=alig_primer_end.replace("\n","")
            alig_primer_end3=alig_primer_end2.replace(" ","")
            alig_primer_seq=primer_blast_data_alig[2]
            alig_primer_seq2=alig_primer_seq.replace(" ","")
            alignment_start=alignment_blast_data_alig[1]
            alignment_start2=alignment_start.replace(" ","")
            alignment_end=alignment_blast_data_alig[3]
            alignment_end2=alignment_end.replace("\n","")
            alignment_end3=alignment_end2.replace(" ","")
            alignment_seq=alignment_blast_data_alig[2]
            alignment_seq2=alignment_seq.replace(" ","")
#if primer_id in Primers_dict:
    #print (Primers_dict[primer_id])
            Primers_dict[primer_id].append(alig_primer_start2)
            Primers_dict[primer_id].append(alig_primer_end3)
            Primers_dict[primer_id].append(alig_primer_seq2)
            Primers_dict[primer_id].append(alignment_start2)
            Primers_dict[primer_id].append(alignment_end3)
            Primers_dict[primer_id].append(alignment_seq2)
            #print (Primers_dict[primer_id])
        Blastn_result_list_alig=[]

    def search_primer_position_reference(self, alignment_file):
        #root=os.getcwd()
        #os.chdir(root+'/'+path)
        #Conservation_scores.create_blast_db("CONSERVATION_SCORES\Virus\RHDV",alignment_file)
        for primers in Primers_dict:
            primer_data=Primers_dict[primers]
            primer_id=primer_data[0]
            primer_sequence=primer_data[1]
            primer_reference=primer_data[2]
            primer_length=primer_data[3]
            #print(primer_sequence)
            RHDV_primer_lit.create_blastn_query(self,"CONSERVATION_SCORES\Virus\RHDV",primer_id,primer_sequence)
            RHDV_primer_lit.run_blastn_ref(self,"CONSERVATION_SCORES\Virus\RHDV","ReferenceRHDV.fasta","query.fasta", primer_id) 
            #RHDV_primer_lit.read_blastn_output_ref(self,"CONSERVATION_SCORES\Virus\RHDV\RESULTS\RHDV_"+ primer_id+".out",primer_id)
        
    def search_primer_position_alignment(self, alignment_file):
        #root=os.getcwd()
        #os.chdir(root+'/'+path)
        #Conservation_scores.create_fasta_ref_from_alig("CONSERVATION_SCORES\Virus\RHDV",alignment_file)
        #Conservation_scores.create_blast_db_alig("CONSERVATION_SCORES\Virus\RHDV",alignment_file)
        #percent=0
        for primers in Primers_dict:
            primer_data=Primers_dict[primers]
            primer_id=primer_data[0]
            primer_sequence=primer_data[1]
            primer_reference=primer_data[2]
            primer_length=primer_data[3]
            #percent=percent+1
            #total_primers=len(Primers_dict)
            #percent_total=round((percent/total_primers)*100,2)
            #print("calculating:"+str(percent_total))
            #if int(primer_id) > 1315:
            #print(primer_id)
            #RHDV_primer_lit.create_blastn_query(self,"CONSERVATION_SCORES\Virus\RHDV",primer_id,primer_sequence)
            #RHDV_primer_lit.run_blastn_alig(self,"CONSERVATION_SCORES\Virus\RHDV","REFblastdbalig.fasta","query.fasta", primer_id)
            RHDV_primer_lit.read_blastn_output_alig(self,"CONSERVATION_SCORES\Virus\RHDV\RESULTS2\/RHDV_"+ primer_id+".out",primer_id)

    def mt_and_gc_calcultation(self):
        for primers in Primers_dict:
            primer_data=Primers_dict[primers]
            primer_id=primer_data[0]
            primer_sequence=primer_data[1]
            primer_length=primer_data[3]
            n_A = primer_sequence.count('A')
            n_T = primer_sequence.count('T')
            n_G = primer_sequence.count('G')
            n_C = primer_sequence.count('C')
            #print(n_A,n_C,n_G,n_T)
            Total_n=n_A+n_T+n_G+n_C
            if Total_n==0:
                pass
            else:
                GC_content = (n_G + n_C)/(n_G + n_C + n_A + n_T) 
                GC_content = GC_content * 100
                GC_content=round(GC_content,2)
                if int(primer_length) <13:
                    Melting = 2*(n_A + n_T) + 4*(n_G + n_C)
                if int(primer_length) >=13:
                    Melting = 64.9 + 41*(n_G + n_C - 16.4)/(n_A + n_T + n_G + n_C)
                    Melting=round(Melting,2)
                if primer_id in Primers_dict:
                    #print(Primers_dict[primer_id])
                    Primers_dict[primer_id].append(GC_content)
                    Primers_dict[primer_id].append(Melting)
            #print(Primers_dict[str(1)])
            #print(GC_content, Melting)

    def calculate_PIS_window_alignment(self,file,start,end,primer_id):
        alignment_data = AlignIO.read(open(file), "fasta")
        alignment=alignment_data[:,start:end]
        #print (alignment)
        summary_align = AlignInfo.SummaryInfo(alignment)
        #print (summary_align)
        alignment_len=alignment.get_alignment_length()
        alignment_number_seqs=len(alignment)
        dict_column_freqs={}
        Number_identical_sites=0
        #Calculates the freqs per column for Indels and nucleotides
        for col in range(alignment_len):
            #print ("Retrieving column:"+str(col))
            alignment_columns=summary_align.get_column(col)
            #print (alignment_columns)
            countA=alignment_columns.count("A")/alignment_number_seqs
            countT=alignment_columns.count("T")/alignment_number_seqs
            countG=alignment_columns.count("G")/alignment_number_seqs
            countC=alignment_columns.count("C")/alignment_number_seqs
            countN=alignment_columns.count("N")/alignment_number_seqs
            countDEL=alignment_columns.count("-")/alignment_number_seqs
            dict_column_freqs[col]=(countA,countT,countG,countC,countN,countDEL)
            #Ref_nucleotide=alignment_columns[0]
            #print (alignment_columns[1])
        #Calculates percentage of identical sites
        for col in dict_column_freqs:
            data_col=dict_column_freqs[col]
            #print (data_col)
            for value in data_col:
                if value==1:
                    Number_identical_sites+=1
                else:
                    pass
        #print (Number_identical_sites)
        if alignment_len >0:
            PIS=(Number_identical_sites/alignment_len)*100
        if alignment_len ==0:
            PIS=0
        #print (dict_column_freqs[0])
        print (PIS)
        #print("Alignment length %i" % alignment.get_alignment_length())
        #if primer_id in Primers_dict:
            #Primers_dict[primer_id].append(PIS)
        PIS=round(PIS,2)
        return PIS
        #for record in alignment:
        #    print(record.seq + " " + record.id)
            #print(record.id+"\n")
    
    def calculate_PPI_window_alignment_fast(self,file,start,end):
        alignment_data = AlignIO.read(open(file), "fasta")
        alignment=alignment_data[:,start:end]
        #print (alignment)
        summary_align = AlignInfo.SummaryInfo(alignment)
        alignment_len=alignment.get_alignment_length()
        pairwise_differences_total_fast=0
        pairwise_equal_fast=0
        pairwise_dif_fast=0
        alignment_number_seqs=len(alignment)
        #Calculates percentage of pairwise identity
        for col in range(alignment_len):
            data_col=summary_align.get_column(col)
            list_nuc=list(data_col)
            array_nuc=np.asarray(list_nuc)
            #print (array_nuc)
            mask = array_nuc[:,None] > array_nuc
            array_ppi_fast = mask[~np.eye(array_nuc.size,dtype=bool)]
            #print (array_ppi_fast)
            pairwise_dif_fast+=np.count_nonzero(array_ppi_fast)
            #pairwise_dif_fast+=array_ppi_fast.sum()
            pairwise_differences_total_fast+=len(array_ppi_fast)/2
            #print (pairwise_equal_fast)
            #print (array_ppi_fast)
            #print (pairwise_equal_fast)
            Ref_nucleotide=data_col[0]
        pairwise_equal_fast+=pairwise_differences_total_fast-pairwise_dif_fast
        if pairwise_differences_total_fast >0:
            PPI=(pairwise_equal_fast/pairwise_differences_total_fast)*100
        elif pairwise_differences_total_fast ==0:
            PPI=0
        PPI=round(PPI,2)
        return PPI

    def run_RHDV_PIS_PPI_and_SCORE(self):
        for primers in Primers_dict:
            primer_data=Primers_dict[primers]
            primer_id=primer_data[0]
            if len(primer_data)==14:
                print ("calculating conservation scores for:"+str(primer_data)+"\n")
                start=int(primer_data[9])
                end=int(primer_data[10])
                if start >end:
                    start=int(primer_data[10])
                    end=int(primer_data[9])
                    primer_type="Reverse"
                elif start <end:
                    start=int(primer_data[9])
                    end=int(primer_data[10])
                    primer_type="Forward"
                #PIS=RHDV_primer_lit.calculate_PIS_window_alignment(self,"Cracas_all_species_alignment_04072023.fasta",start,end,primer_id)
                PPI=RHDV_primer_lit.calculate_PPI_window_alignment_fast(self,"RHDV_all.fasta",start,end)
                if primer_type=="Forward":
                    PPI3=RHDV_primer_lit.calculate_PPI_window_3end_alignment_fast(self,"RHDV_all.fasta",start,end,primer_type)
                elif primer_type=="Reverse":
                    PPI3=RHDV_primer_lit.calculate_PPI_window_3end_alignment_fast(self,"RHDV_all.fasta",start,end,primer_type)
            elif len(primer_data) <13:
                PPI=0
                PPI3=0
                #PIS=0
                primer_type="not found"
                #del Primers_dict[primers]
                pass
            Score=((PPI+PPI3)/2)
            Score=round(Score,2)
            if primer_id in Primers_dict:
                Primers_dict[primer_id].append(PPI)
                Primers_dict[primer_id].append(primer_type)
                Primers_dict[primer_id].append(PPI3)
                #Primers_dict[primer_id].append(PIS)
                Primers_dict[primer_id].append(Score)
        #print(Primers_dict[primer_id])

    def calculate_PPI_window_3end_alignment_fast(self,file,start,end,primer_type):
        #calculate pairwise identity for the last 3 end primer nucleotides
        alignment_data = AlignIO.read(open(file), "fasta")
        if primer_type=="Forward":
            alignment=alignment_data[:,end-3:end]
            #print (alignment)
        elif primer_type=="Reverse":
            alignment=alignment_data[:,start:end]
            alignment_rev_comp=MultipleSeqAlignment([])
            for record in alignment:
                record2seq=record.seq.reverse_complement()
                record2_name=record.id
                seq=SeqRecord(record2seq, id=record2_name)
                alignment_rev_comp.append(seq)
            alignment_5end=alignment_rev_comp[:,0:3]
            #print (alignment_5end)
            #print(alignment_rev_comp)
            #print ("original seq:"+record+"rev comp:"+record2seq)
        #print("\ncalculation completed\n")
        summary_align = AlignInfo.SummaryInfo(alignment)
        alignment_len=alignment.get_alignment_length()
        pairwise_differences_total_fast=0
        pairwise_equal_fast=0
        pairwise_dif_fast=0
        alignment_number_seqs=len(alignment)
        #Calculates percentage of pairwise identity
        for col in range(alignment_len):
            data_col=summary_align.get_column(col)
            list_nuc=list(data_col)
            array_nuc=np.asarray(list_nuc)
            #print (array_nuc)
            mask = array_nuc[:,None] > array_nuc
            array_ppi_fast = mask[~np.eye(array_nuc.size,dtype=bool)]
            #print (array_ppi_fast)
            pairwise_dif_fast+=np.count_nonzero(array_ppi_fast)
            #pairwise_dif_fast+=array_ppi_fast.sum()
            pairwise_differences_total_fast+=len(array_ppi_fast)/2
            #print (pairwise_equal_fast)
            #print (array_ppi_fast)
            #print (pairwise_equal_fast)
            Ref_nucleotide=data_col[0]
        pairwise_equal_fast+=pairwise_differences_total_fast-pairwise_dif_fast
        if pairwise_differences_total_fast >0:
            PPI3=(pairwise_equal_fast/pairwise_differences_total_fast)*100
        elif pairwise_differences_total_fast ==0:
            PPI3=0
        #print ("\nNumber of comparisons "+str(pairwise_differences_total_fast))
        #print (str(PPI)+"- Number of total comparisons:"+str(pairwise_differences_total_fast))
        PPI3=round(PPI3,2)
        return PPI3

    def export_rhdv_primers_db(self,output_file):
            file_name_validated=open("CONSERVATION_SCORES\CSVs" + output_file+"_validated.csv","w")
            file_name_validated.write("Primer ID;Primer Sequence;Reference;Length;PrimerStart;PrimerEnd;BlastPrimerStart;BlastPrimerEnd;BlastPrimerSeq;AligStart;AligEnd;AligSeq;Gc Content;Melting;PPI;Primer Type;PPI3;Conservation Score\n")
            file_name_not_validated=open("CONSERVATION_SCORES\CSVs" + output_file+"_not_validated.csv","w")
            file_name_not_validated.write("Primer ID;Primer Sequence;Reference;Length;PrimerStart;PrimerEnd;BlastPrimerStart;BlastPrimerEnd;BlastPrimerSeq;AligStart;AligEnd;AligSeq;Gc Content;Melting;PPI;Primer Type;PPI3;Conservation Score\n")
            for primer_data in Primers_dict:
                final_data=str(Primers_dict[primer_data])
                final_data2=final_data.replace(",",";")
                final_data3=final_data2.replace("[","")
                final_data4=final_data3.replace("]","")
                final_data5=final_data4.replace("'","")
                final_data6=final_data5.replace(" ","")
                if len (Primers_dict[primer_data])==18:
                    file_name_validated.write(final_data6+"\n")
                if len (Primers_dict[primer_data])<18:
                    file_name_not_validated.write(final_data6+"\n")
                #print(Primers_dict[primer_data])
            file_name_validated.close()
            file_name_not_validated.close()

    def calculate_primer_combinations(self):
        for primer in Primers_dict2:
            primer_data=Primers_dict2[primer]
            primer1_name=primer_data[0]
            primer1_seq=primer_data[1]
            primer1_ref=primer_data[2]
            primer1_start=primer_data[4]
            primer1_end=int(primer_data[5])
            primer1_gc=primer_data[6]
            primer1_mt=primer_data[7]
            primer1_type=primer_data[8]
            primer1_scores=float(primer_data[9])
            #print(primer1_type)
            if primer1_type=="Forward":
                for primer in Primers_dict2:
                    primer_data=Primers_dict2[primer]
                    primer2_name=primer_data[0]
                    primer2_seq=primer_data[1]
                    primer2_ref=primer_data[2]
                    primer2_start=int(primer_data[4])
                    primer2_end=primer_data[5]
                    primer2_gc=primer_data[6]
                    primer2_mt=primer_data[7]
                    primer2_type=primer_data[8]
                    primer2_scores=float(primer_data[9])
                    #print(primer2_type)
                    if primer2_type=="Reverse":
                        if primer1_end<primer2_start:
                            #length=primer2_end-primer1_start
                            #print((length))
                            #PCR product with 70 to 1000 nucleotides of length
                            if 70<=primer2_end-primer1_start<=1000:
                                #print(primer1_name,primer2_name)
                                if (primer1_name,primer2_name) not in Primers_dict_comb:
                                    primer_comb_scores=(primer1_scores+primer2_scores)/2
                                    primer_comb_scores=round(primer_comb_scores,2)
                                    amplicon_length=primer2_end-primer1_start
                                    gc_content=(float(primer1_gc)+float(primer2_gc))/2
                                    gc_content=round(gc_content,2)
                                    melting=(float(primer1_mt)+float(primer2_mt))/2
                                    melting=round(melting,2)
                                    primer1_self_fold=primer3.calc_hairpin(primer1_seq)
                                    primer2_self_fold=primer3.calc_hairpin(primer2_seq)
                                    primer1_homo_fold=primer3.calc_homodimer(primer1_seq)
                                    primer2_homo_fold=primer3.calc_homodimer(primer2_seq)
                                    primer_dimer_fold=primer3.calc_heterodimer(primer1_seq,primer2_seq)
                                    primer1_self_fold_list=str(primer1_self_fold).split(",")
                                    primer2_self_fold_list=str(primer2_self_fold).split(",") 
                                    primer1_homo_fold_list=str(primer1_homo_fold).split(",")
                                    primer2_homo_fold_list=str(primer2_homo_fold).split(",")
                                    primer_dimer_fold_list=str(primer_dimer_fold).split(",")
                                    Primers_dict_comb[primer1_name,primer2_name]=[primer1_name,primer2_name,primer1_seq,primer2_seq,primer1_gc,primer2_gc,gc_content,primer1_mt,
                                                                                  primer2_mt,melting,primer1_self_fold_list[2],primer2_self_fold_list[2],primer1_homo_fold_list[2],
                                                                                  primer2_homo_fold_list[2],primer_dimer_fold_list[2],primer1_start,primer2_end,
                                                                                  amplicon_length,primer1_scores,primer2_scores,primer_comb_scores,primer1_ref,primer2_ref]
                                    print(primer1_name)
                                    #print(primer1_homo_fold,primer_dimer_fold)
                                    #print(Primers_dict_comb)

    def calculate_fold(self):
        for primers in Primers_dict_comb:
            primer_data=Primers_dict_comb[primers]
            primer1_name=primer_data[0]
            primer2_name=primer_data[1]
            primer_f_self_a=primer_data[10]
            primer_f_self=primer_f_self_a.replace("dg=","")
            primer_r_self_a=primer_data[11]
            primer_r_self=primer_r_self_a.replace("dg=","")
            primer_f_homo_a=primer_data[12]
            primer_f_homo=primer_f_homo_a.replace("dg=","")
            primer_r_homo_a=primer_data[13]
            primer_r_homo=primer_r_homo_a.replace("dg=","")
            dimer_a=primer_data[14]
            dimer=dimer_a.replace("dg=","")
            if 0.0>=float(primer_f_self)>=(-2000):
                Perc_self_f=(((0.05)*float(primer_f_self))+100)
            else: Perc_self_f=0
            if 0.0>=float(primer_r_self)>=(-2000):
                Perc_self_r=(((0.05)*float(primer_r_self))+100)
            else: Perc_self_r=0
            if 0.0>=float(primer_f_homo)>=(-5000):
                Perc_homo_f=(((0.02)*float(primer_f_homo))+100)
            else: Perc_homo_f=0
            if 0.0>=float(primer_r_homo)>=(-5000):
                Perc_homo_r=(((0.02)*float(primer_r_homo))+100)
            else: Perc_homo_r=0
            if 0.0>=float(dimer)>=(-5000):
                Dimer_perc=(((0.02)*float(dimer))+100)
            else: Dimer_perc=0
            Perc_self_m=((Perc_self_f+Perc_self_r)/2)
            Perc_homo_m=((Perc_homo_f+Perc_homo_r)/2)
            Perc_fold=((Perc_self_m+Perc_homo_m+Dimer_perc)/3)
            Perc_fold=round(Perc_fold,2)
            if (primer1_name,primer2_name) in Primers_dict_comb:
                #Primers_dict_comb[primer1_name,primer2_name].append(Perc_self_f)
                #Primers_dict_comb[primer1_name,primer2_name].append(Perc_self_r)
                #Primers_dict_comb[primer1_name,primer2_name].append(Perc_homo_f)
                #Primers_dict_comb[primer1_name,primer2_name].append(Perc_homo_r)
                #Primers_dict_comb[primer1_name,primer2_name].append(Dimer_perc)
                Primers_dict_comb[primer1_name,primer2_name].append(Perc_fold)
        
    def regions(self):
        for primers in Primers_dict_comb:
            primer_data=Primers_dict_comb[primers]
            primer1_name=primer_data[0]
            primer2_name=primer_data[1]
            primer_f_start=primer_data[15]
            primer_r_end=primer_data[16]
            if 438>=float(primer_f_start)>=0:
                primer_f_region="p16"
            elif 1110>=float(primer_f_start)>=439:
                primer_f_region="p23-p26"
            elif 2163>=float(primer_f_start)>=1111:
                primer_f_region="2C"
            elif 2988>=float(primer_f_start)>=2164:
                primer_f_region="p29"
            elif 3333>=float(primer_f_start)>=2999:
                primer_f_region="VPg"
            elif 3762>=float(primer_f_start)>=3334:
                primer_f_region="3C"
            elif 5310>=float(primer_f_start)>=3763:
                primer_f_region="RdRp"
            elif 7024>=float(primer_f_start)>=5311:
                primer_f_region="VP60"
            elif 7378>=float(primer_f_start)>=7041:
                primer_f_region="VP10"
            elif 7042>=float(primer_f_start)>=7025:
                primer_f_region="VP60/VP10"
            if 438>=float(primer_r_end)>=0:
                primer_r_region="p16"
            elif 1110>=float(primer_r_end)>=439:
                primer_r_region="p23-p26"
            elif 2163>=float(primer_r_end)>=1111:
                primer_r_region="2C"
            elif 2988>=float(primer_r_end)>=2164:
                primer_r_region="p29"
            elif 3333>=float(primer_r_end)>=2999:
                primer_r_region="VPg"
            elif 3762>=float(primer_r_end)>=3334:
                primer_r_region="3C"
            elif 5310>=float(primer_r_end)>=3763:
                primer_r_region="RdRp"
            elif 7024>=float(primer_r_end)>=5311:
                primer_r_region="VP60"
            elif 7378>=float(primer_r_end)>=7042:
                primer_r_region="VP10"
            elif 7042>=float(primer_r_end)>=7025:
                primer_r_region="VP60/VP10"
            if primer_f_region==primer_r_region:
                region=primer_f_region
            else:
                region=primer_f_region + " - " + primer_r_region
            if (primer1_name,primer2_name) in Primers_dict_comb:
                Primers_dict_comb[primer1_name,primer2_name].append(primer_f_region)
                Primers_dict_comb[primer1_name,primer2_name].append(primer_r_region)
                Primers_dict_comb[primer1_name,primer2_name].append(region)

    def export_primer_combinations(self,output_file):
        file_name=open("CONSERVATION_SCORES\CSVs" + output_file,"w")
        file_name.write("ID Forward;Forward Ref;ID Reverse;Reverse Ref;Forward Sequence;Reverse Sequence;Forward GC;Reverse GC;Comb GC;Forward Tm;Reverse Tm;Comb Tm;Forward Homo;Reverse Homo;Forward Self;Reverse Self;Dimer;Fold;Forward Start;Reverse End;Reverse Region;Amplicon Length;Region;Primer Combination Score\n")
        for primer in Primers_dict_comb:
            Primer_Forward_ID=Primers_dict_comb[primer][0]
            Primer_Forward_Ref=Primers_dict_comb[primer][21]
            Primer_Reverse_ID=Primers_dict_comb[primer][1]
            Primer_Reverse_Ref=Primers_dict_comb[primer][22]
            Primer_Forward_Seq=Primers_dict_comb[primer][2]
            Primer_Reverse_Seq=Primers_dict_comb[primer][3]
            Primer_Forward_GC=Primers_dict_comb[primer][4]
            Primer_Reverse_GC=Primers_dict_comb[primer][5]
            GC_comb=Primers_dict_comb[primer][6]
            Primer_Forward_Tm=Primers_dict_comb[primer][7]
            Primer_Reverse_Tm=Primers_dict_comb[primer][8]
            Melting_temp=Primers_dict_comb[primer][9]
            Forward_homo_fold=Primers_dict_comb[primer][10]
            Forward_homo_fold2=Forward_homo_fold.replace("dg=","")
            Reverse_homo_fold=Primers_dict_comb[primer][11]
            Reverse_homo_fold2=Reverse_homo_fold.replace("dg=","")
            Forward_self_fold=Primers_dict_comb[primer][12]
            Forward_self_fold2=Forward_self_fold.replace("dg=","")
            Reverse_self_fold=Primers_dict_comb[primer][13]
            Reverse_self_fold2=Reverse_self_fold.replace("dg=","")
            Dimer=Primers_dict_comb[primer][14]
            Dimer2=Dimer.replace("dg=","")
            Fold=Primers_dict_comb[primer][23]
            Primer_Forward_Start=Primers_dict_comb[primer][15]
            #Primer_Forward_Region=Primers_dict_comb[primer][24]
            Primer_Reverse_End=Primers_dict_comb[primer][16] 
            #Primer_Reverse_Region=Primers_dict_comb[primer][25]           
            Amplicon_length=Primers_dict_comb[primer][17]
            Region=Primers_dict_comb[primer][26]
            Primer_Combination_Score=Primers_dict_comb[primer][20]
            final_data=str(Primer_Forward_ID)+";"+str(Primer_Forward_Ref)+";"+str(Primer_Reverse_ID)+";"+str(Primer_Reverse_Ref)+";"+str(Primer_Forward_Seq)+";"+str(Primer_Reverse_Seq)+";"+str(Primer_Forward_GC)+";"+str(Primer_Reverse_GC)+";"+str(GC_comb)+";"+str(Primer_Forward_Tm)+";"+str(Primer_Reverse_Tm)+";"+str(Melting_temp)+";"+str(Forward_homo_fold2)+";"+str(Reverse_homo_fold2)+";"+str(Forward_self_fold2)+";"+str(Reverse_self_fold2)+";"+str(Dimer2)+";"+str(Fold)+";"+str(Primer_Forward_Start)+";"+str(Primer_Reverse_End)+";"+str(Amplicon_length)+";"+str(Region)+";"+str(Primer_Combination_Score)
            file_name.write(final_data+"\n")
        #print(Primers_dict_comb)
        file_name.close()

    def generate_in_silico_primers(self,input_file,name_analysis,range_start,range_end,step):
        '''input_file is the name of the file with reference genome
        name_analysis is the name of the analysed species
        window_len_range is the length of the primers range: 18-30 nuc
        step is the leap between the generated primers: 1 nuc'''
        reference_data = SeqIO.parse(open(input_file), "fasta")
        for reference_record in reference_data:
            reference_record_id=reference_record.id
            if reference_record_id=="Balanus_trigonus_-_JQ035523.1":
                reference_record_name=reference_record_id
                reference_record_seq_forward=reference_record.seq
                reference_record_seq_reverse=reference_record.seq.reverse_complement()
                reference_len=len(reference_record_seq_forward)
        in_silico_primer_id=0
        for primer_length in range(range_start,range_end+1,1):  
            for window_position in range(0,reference_len-primer_length+1,step):
                #print (window_position)
                start_pos=window_position
                end_pos=start_pos+primer_length
                in_silico_seq_forward=str(reference_record_seq_forward[start_pos:end_pos])
                primer_length_final_forward=len(in_silico_seq_forward)
                in_silico_seq_reverse=str(reference_record_seq_reverse[start_pos:end_pos])
                primer_length_final_reverse=len(in_silico_seq_reverse)
                if (start_pos,end_pos,primer_length_final_forward,"forward") not in Primers_dict_in_silico:
                    in_silico_primer_id+=1
                    Primers_dict_in_silico[start_pos,end_pos,primer_length_final_forward,"forward"]=[in_silico_primer_id,in_silico_seq_forward,"in silico primer forward",primer_length_final_forward,start_pos,end_pos]
                if (start_pos,end_pos,primer_length_final_reverse,"reverse") not in Primers_dict_in_silico:
                    in_silico_primer_id+=1
                    Primers_dict_in_silico[start_pos,end_pos,primer_length_final_reverse,"reverse"]=[in_silico_primer_id,in_silico_seq_reverse,"in silico primer reverse",primer_length_final_reverse,start_pos,end_pos]
                elif (start_pos,end_pos,primer_length_final_forward,"forward") in Primers_dict_in_silico:
                    print("in silico primer forward already exists")
                elif (start_pos,end_pos,primer_length_final_reverse,"reverse") in Primers_dict_in_silico:
                    print("in silico primer reverse already exists")
        for primer in Primers_dict_in_silico:
            print(Primers_dict_in_silico[primer])

    def export_in_silico_primers(self,output_file):
        file_name=open("CONSERVATION_SCORES\CSVs" + output_file,"w")
        file_name.write("ID;Sequence;Type;Length;Start;End\n")
        for primer in Primers_dict_in_silico:
            Primers_in_silico_ID=Primers_dict_in_silico[primer][0]
            Primers_in_silico_seq=Primers_dict_in_silico[primer][1]
            Primers_in_silico_type=Primers_dict_in_silico[primer][2]
            Primers_in_silico_length=Primers_dict_in_silico[primer][3]
            Primers_in_silico_start=Primers_dict_in_silico[primer][4]
            Primers_in_silico_end=Primers_dict_in_silico[primer][5]
            final_data=str(Primers_in_silico_ID)+";"+str(Primers_in_silico_seq)+";"+str(Primers_in_silico_type)+";"+str(Primers_in_silico_length)+";"+str(Primers_in_silico_start)+";"+str(Primers_in_silico_end)
            file_name.write(final_data+"\n")
            #print(Primers_dict_in_silico[primer])
        file_name.close()

Conservation_scores=RHDV_primer_lit()
#Conservation_scores.parse_primers_RHDV("finalprimersrhdv03072023.csv")
Conservation_scores.parse_primers_final("Final_AROLIT.csv")
#Conservation_scores.parse_primers_in_silico("CSVsInSilicoPrimers1704.csv",7437)
#Conservation_scores.parse_primers_RHDV_test("Test.csv")
#Conservation_scores.create_blast_db("CONSERVATION_SCORES\Virus\RHDV","ReferenceRHDV.fasta")
#Conservation_scores.create_blast_db_alig("CONSERVATION_SCORES\Virus\RHDV",alignment_file)
#Conservation_scores.run_blastn("CONSERVATION_SCORES\Virus\RHDV","RHDV_all.fasta")
#Conservation_scores.search_primer_position_reference("ReferenceRHDV.fasta")
#Conservation_scores.search_primer_position_alignment("RHDV_all.fasta")
#Conservation_scores.mt_and_gc_calcultation()
#Conservation_scores.run_RHDV_PIS_PPI_and_SCORE()
#Conservation_scores.export_rhdv_primers_db("PrimersInSilico150723Gap2WS18")
Conservation_scores.calculate_primer_combinations()
Conservation_scores.calculate_fold()
Conservation_scores.regions()
Conservation_scores.export_primer_combinations("PrimerCombinationsAROLit_SelfHomoDimer.csv")
#Conservation_scores.generate_in_silico_primers("CONSERVATION_SCORES\Virus\RHDV\Balanus_trigonus_-_JQ035523.1.fasta","RHDV",18,30,1)
#Conservation_scores.export_in_silico_primers("InSilicoPrimersCracas.csv")