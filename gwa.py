import numpy as np
import ftplib
import gzip
import iupred2a_lib
from itertools import groupby
from Bio import SeqIO
from Bio.Seq import Seq
import  matplotlib.pyplot as plt
import pandas as pd
ftp = ftplib.FTP("ftp.uniprot.org")
ftp.login()
ftp_file = "pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
ftp.retrbinary("RETR " + ftp_file , open("human_proteome.fasta", 'wb').write)
#ftp.retrbinary("RETR " + ftp_file ,open("human_proteome.fasta", 'wb').write)
#df = pd.DataFrame(columns=[ 'Fraction Prolines'])
write_to_list=[]
write_id_list=[]
with gzip.open("human_proteome.fasta", mode="rt") as file_handler:
#with open("test.fasta", "rt") as file_handler:     
     for record in SeqIO.parse(file_handler, "fasta"):   
         list_record = list(record.seq)
         
         list_record.append('G') 
         #print(list_record)
         #print("SEQUENCE UNDER DISORDER PREDICTION ALGORITHM NOW:", record.seq)    
         iupred_scores = iupred2a_lib.iupred(record.seq, "long")
         #print (len(record.seq), type(iupred_scores))
         #print((iupred_scores)) 
         #for x in iupred_scores:
             #print((x))
         
         #plt.figure()
         #plt.xlabel('Residue Number')
         #plt.ylabel('Disorder Score')
         #plt.plot(iupred_scores[0] )  
         #plt.show(block=False)
         #plt.pause(3)
         #plt.close() 
         
         iupred_score_list=list(iupred_scores[0][:])
         iupred_score_list.append(0)
         temporary_list=[]
         final_list=[]
         #print(len(mylist))
         for x in range(0, len(list_record)):
             #print(type(iupred_score_list[x]))      
             if iupred_score_list[x] > 0.5 : 
       
                temporary_list.append(list_record[x])
       #final_list.append(temporary_list)
       #print(x)
             if iupred_score_list[x] <= 0.5  :
                #print(len(temporary_list))   
                #if len(temporary_list) != 0:
                #print(temporary_list)
                final_list.append(list(temporary_list))
                #print(final_list)
                temporary_list.clear()
         for elem in final_list :
             if (len(elem)) >= 20:
                 print( "your IDR search gives the following sequence region", (''.join(str(items) for items in elem) ))       
                 write_to_list.append(''.join(str(items) for items in elem))
                 write_id_list.append(record.id)
print("done with record IDR sorting")
#print(write_to_list)
#with open("analisis_proline_in_seq.txt", "w") as writefile:
proline_analysis_list=[]

for item in write_to_list   : 
    idr_seq_record=Seq(item)
    proline_analysis_list.append(idr_seq_record.count('P')/len(idr_seq_record))
       
    #df = df.append(pd.DataFrame([proline_analisys_temp_list], columns=['Fraction Prolines']),     ignore_index=True)
    #df = pd.concat(proline_analisys_temp_list)
    
    #proline_analisys_temp_list.clear()
#writefile.close()
print(proline_analysis_list)    
#with open('proline.txt', 'w') as f:
#    for item in proline_analysis_list:
#        f.write("%s\n" % item)

#with open('idr_sequences_human_proteome.txt', 'w') as write_file:
#    for item in write_to_list:
#        write_file.write("%s\n" % item)
with open('idr_sequences_ids_human_proteome.txt', 'w') as write_file:
    for item in write_id_list:
        write_file.write("%s\n" % item)
