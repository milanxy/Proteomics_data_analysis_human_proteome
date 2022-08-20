import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns
import plotly.express as px
inner_length=
final_count=[]
state=[]
count_residues=0
fraction_fp =0
fraction_pf =0
fraction_wp =0
fraction_pw =0
fraction_py=0
fraction_yp=0
fraction_proline1=0
fraction_proline2=0
fraction_proline3=0
aa_list =['A', 'C', 'D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
calculated_abundance={'A':[], 'C':[], 'D':[],'E':[],'F':[],'G':[],'H':[],'I':[],'K':[],'L':[],'M':[],'N':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'V':[],'W':[],'Y':[]}
data = {'Fraction_Proline1':[], 'Fraction_Proline2':[], 'Fraction_Proline3':[], 'Fraction_PhePro':[],'Fraction_ProPhe':[], 'Fraction_TrpPro':[],'Fraction_ProTrp':[], 'Fraction_ProTyr':[], 'Fraction_TyrPro':[],'seq_length':[],'total_prolines':[] }
fraction_aromatic={'tot_frac_pf':[], 'tot_frac_pw':[], 'tot_frac_py':[], 'tot_frac_all_aromatic':[], 'tot_pro':[], 'seq_length':[]}
#aromatic_data = pd.DataFrame(columns=["Fraction_PheProline", "Fraction_ProlinePhe", "Fraction_TrpProline", "Fraction_ProlineTrp"])
with open("../../../idr_sequences_human_proteome.txt", "r") as readfile:
   for record in readfile:
     final_count.clear()
     state.clear()
     fraction_fp =0
     fraction_pf =0
     fraction_wp =0
     fraction_pw =0
     fraction_py=0
     fraction_yp=0
     fraction_proline1=0
     fraction_proline2=0
     fraction_proline3=0
     if len(record) > inner_length and len(record) < inner_length+100:
       #print(len(record))
       idr_seq_record=Seq(record)
       for item in aa_list:
           calculated_abundance[item].append(100*idr_seq_record.count(item)/len(record))
       for item in record:          
           count_residues =count_residues+1
      
       record =record+'G'
       #record=("".join((record, 'G')))
       old_item= record[0]
       count_seq_proline_adjacent=0
      
#############################################################################
       for item in record:
           if item == 'P' and old_item !='P':
              count_seq_proline_adjacent =1
              old_item = item
           elif item == 'P' and old_item =='P':
              count_seq_proline_adjacent = count_seq_proline_adjacent +1
              old_item = item
           elif item != 'P' and old_item =='P':
                final_count.append(count_seq_proline_adjacent)
                count_seq_proline_adjacent=0
                old_item = item
           elif item != 'P' and old_item !='P':
                old_item = item
       for item in final_count:
           if item == 1:
              fraction_proline1 = fraction_proline1+ 1
           elif item == 2 :
              fraction_proline2 = fraction_proline2+ 1
           elif item == 3 :
              fraction_proline3 = fraction_proline3+ 1
#########################################################################       



       for index, item in enumerate(record):
           if index == 0:
              item_next =record[index+1]
              if item == 'P'  and item_next != 'P':
                 state.append(1)
              else: 
                 state.append(0)
           elif index == len(record)-1:
               old_item = record[index-1]
               if old_item != 'P' and item =='P' :
                 state.append(1)
               else:
                 state.append(0)
           elif index !=0 and index != len(record)-1:
           
              old_item = record[index-1]
              item_next =record[index+1] 
              if item == 'P'  and item_next != 'P' and old_item != 'P':
                 state.append(1)
              else:
                 state.append(0)
           
       
       for index, item in enumerate(record): 
          if index == 0:
            if state[index] == 1 and record[index+1] == 'F': 
              fraction_pf = fraction_pf  +1
            
          elif index == len(record)-1:  
            if state[index] == 1 and record[index-1] == 'F':
               fraction_fp = fraction_fp  +1
          elif index != 0 and index != len(record)-1:   
            if state[index] == 1 and record[index+1] == 'F' and record[index-1] != 'F': 
              fraction_pf = fraction_pf  +1
            elif state[index] == 1 and record[index-1] == 'F' and record[index+1] != 'F':
              fraction_fp = fraction_fp  +1
            elif state[index] == 1 and record[index-1] == 'F' and record[index+1] == 'F':
                 fraction_fp = fraction_fp  +1
                 fraction_pf =fraction_pf   +1
       for index, item in enumerate(record):
          if index == 0:
            if state[index] == 1 and record[index+1] == 'W':
              fraction_pw = fraction_pw  +1

          elif index == len(record)-1:
            if state[index] == 1 and record[index-1] == 'W':
               fraction_wp = fraction_wp  +1
          elif index != 0 and index != len(record)-1:
            if state[index] == 1 and record[index+1] == 'W' and record[index-1] != 'W':
              fraction_pw = fraction_pw  +1
            elif state[index] == 1 and record[index-1] == 'W' and record[index+1] != 'W':
              fraction_wp = fraction_wp  +1
            elif state[index] == 1 and record[index-1] == 'W' and record[index+1] == 'W':
                 fraction_wp = fraction_wp  +1
                 fraction_pw = fraction_pw  +1
       for index, item in enumerate(record):
          if index == 0:
            if state[index] == 1 and record[index+1] == 'Y':
              fraction_py = fraction_py  +1
       
            if index == len(record)-1:
               if state[index] == 1 and record[index-1] == 'Y':
                  fraction_yp = fraction_yp  +1
          elif index != 0 and index != len(record)-1:
            if state[index] == 1 and record[index+1] == 'Y' and record[index-1] != 'Y':
              fraction_py = fraction_py  +1
            elif state[index] == 1 and record[index-1] == 'Y'  and record[index+1] != 'Y':
              fraction_yp = fraction_yp  +1

            elif state[index] == 1 and record[index-1] == 'Y' and record[index+1] == 'Y':
                 fraction_yp = fraction_yp  +1
                 fraction_py = fraction_py  +1
       idr_seq=Seq(record)
       #fraction_fp = idr_seq.count_overlap("FP")
       #fraction_pf = idr_seq.count_overlap("PF")
       #fraction_wp = idr_seq.count_overlap("WP")
       #fraction_pw = idr_seq.count_overlap("PW")
       #fraction_py = idr_seq.count_overlap("PY")
       #fraction_yp = idr_seq.count_overlap("YP")
       #if idr_seq.count_overlap("P") != 0:
###########################################################default################################################
       data['Fraction_Proline1'].append((fraction_proline1/len(record))*100)
       data['Fraction_Proline2'].append((fraction_proline2*2/len(record))*100)
       data['Fraction_Proline3'].append((fraction_proline3*3/len(record))*100)
       data['Fraction_PhePro'].append((fraction_fp/(len(record)-1))*100)
       data['Fraction_ProPhe'].append((fraction_pf/(len(record)-1))*100)
       data['Fraction_TrpPro'].append((fraction_wp/(len(record)-1))*100)
       data['Fraction_ProTrp'].append((fraction_pw/(len(record)-1))*100)
       data['Fraction_TyrPro'].append((fraction_yp/(len(record)-1))*100)
       data['Fraction_ProTyr'].append((fraction_py/(len(record)-1))*100)
       data['seq_length'].append(len(record)-1)
       data[ 'total_prolines'].append(100*(idr_seq.count_overlap("P")/(len(record)-1)))
################################################################################################################
       if fraction_proline1 != 0:
          fraction_aromatic['tot_frac_pf'].append(100*float(fraction_pf+fraction_fp)/fraction_proline1)
          fraction_aromatic['tot_frac_pw'].append(100*float(fraction_pw+fraction_wp)/fraction_proline1)
          fraction_aromatic['tot_frac_py'].append(100*float(fraction_py+fraction_yp)/fraction_proline1)
          fraction_aromatic['tot_frac_all_aromatic'].append(100*(float(fraction_pf+fraction_fp+fraction_pw+fraction_wp+fraction_py+fraction_yp)/fraction_proline1))
          fraction_aromatic['tot_pro'].append(100*float(fraction_proline1)/(len(record)-1))
          fraction_aromatic['seq_length'].append(len(record)-1)
       elif fraction_proline1 == 0:
          fraction_aromatic['tot_frac_pf'].append(0)
          fraction_aromatic['tot_frac_pw'].append(0)
          fraction_aromatic['tot_frac_py'].append(0)
          fraction_aromatic['tot_frac_all_aromatic'].append(0)
          fraction_aromatic['tot_pro'].append(0)
          fraction_aromatic['seq_length'].append(len(record)-1)
df = pd.DataFrame(data)
print(df.shape[0])

df_frac_aromatic = pd.DataFrame(fraction_aromatic)
       #idr_seq=Seq(record)
       #fraction_fp=idr_seq.count_overlap("FP")/len(record)
       #fraction_pf=idr_seq.count_overlap("PF")/len(record)
       #fraction_wp=idr_seq.count_overlap("WP")/len(record)
       #fraction_pw=idr_seq.count_overlap("PW")/len(record)
#df['Fraction_Proline1'] = df['Fraction_Proline1'].replace({0:np.nan})
#df['Fraction_Proline2'] = df['Fraction_Proline2'].replace({0:np.nan})
#df['Fraction_Proline3'] = df['Fraction_Proline3'].replace({0:np.nan})
df['Fraction_PhePro'] = df['Fraction_PhePro'].replace({0:np.nan})
df['Fraction_ProPhe'] = df['Fraction_ProPhe'].replace({0:np.nan})
df['Fraction_TrpPro'] = df['Fraction_TrpPro'].replace({0:np.nan})
df['Fraction_ProTrp'] = df['Fraction_ProTrp'].replace({0:np.nan})
df['Fraction_TyrPro'] = df['Fraction_TyrPro'].replace({0:np.nan})
df['Fraction_ProTyr'] = df['Fraction_ProTyr'].replace({0:np.nan})
#print (count_residues)  
df_frac_aromatic['tot_frac_pf']=df_frac_aromatic['tot_frac_pf'].replace({0:np.nan})
df_frac_aromatic['tot_frac_pw']=df_frac_aromatic['tot_frac_pw'].replace({0:np.nan}) 
df_frac_aromatic['tot_frac_py']=df_frac_aromatic['tot_frac_py'].replace({0:np.nan})
print(df["Fraction_Proline1"].mean(), df["Fraction_Proline2"].mean(), df["Fraction_Proline3"].mean(), df["Fraction_PhePro"].mean(),df["Fraction_ProPhe"].mean(), df["Fraction_TrpPro"].mean()  , df["Fraction_ProTrp"].mean(), df["Fraction_TyrPro"].mean(), df["Fraction_ProTyr"].mean()  )


dfca=pd.DataFrame(calculated_abundance)

for item in aa_list:
    #dfca[item] = dfca[item].replace({0:np.nan})
    print(item, dfca[item].mean(), dfca[item].std()/2)
#for item in aa_list:
#    print("Abundace of", item, "in IDR regions (%)","=", 100*(dataframe_calculated_abundance[item].mean))
effective_length = int(inner_length) + 50
#proline_abundance = 100*(dataframe_calculated_abundance['P'].sum())/count_residues
#f = open("sample.txt", "a")

#f.write(f"{effective_length} {proline_abundance}\n")

#f.close()

#f = open("proline_single.txt", "a")
#f.write(f"{effective_length} {0.25*(df.loc[3].values[1]+df.loc[4].values[1]+df.loc[5].values[1]+df.loc[6].values[1])}\n")
#f.close()
 ########################new files commented out now as they may get overwritten############################
#f = open("proline_abundance.txt", "a")
#f.write(f"{effective_length} {dfca['P'].mean()} {0.5* dfca['P'].std()}\n")
#f.close()
f = open("proline_single.txt", "a")
f.write(f"{effective_length} {df['Fraction_Proline1'].mean()} {0.5* df['Fraction_Proline1'].std()}\n")
f.close()

f = open("proline_single_median.txt", "a")
f.write(f"{effective_length} {df['Fraction_Proline1'].median()} {0.5* df['Fraction_Proline1'].std()}\n")
f.close()

f = open("proline_double.txt", "a")
f.write(f"{effective_length} {df['Fraction_Proline2'].mean()} {0.5* df['Fraction_Proline2'].std()}\n")
f.close()


f = open("proline_double_median.txt", "a")
f.write(f"{effective_length} {df['Fraction_Proline2'].median()} {0.5* df['Fraction_Proline2'].std()}\n")
f.close()

f = open("proline_triple.txt", "a")
f.write(f"{effective_length} {df['Fraction_Proline3'].mean()} {0.5* df['Fraction_Proline3'].std()}\n")
f.close()

f = open("proline_triple_median.txt", "a")
f.write(f"{effective_length} {df['Fraction_Proline3'].median()} {0.5* df['Fraction_Proline3'].std()}\n")
f.close()

f = open("proline_total.txt", "a")
f.write(f"{effective_length} {df['total_prolines'].mean()} {0.5* df['total_prolines'].std()}\n")
f.close()

f = open("proline_total_median.txt", "a")
f.write(f"{effective_length} {df['total_prolines'].median()} {0.5* df['total_prolines'].std()}\n")
f.close()

f = open("seq_with_zero_prolines_with_length.txt", "a")
f.write(f"{effective_length} {(df['total_prolines']==0).sum()}\n")
f.close()

f = open("seq_with_zero_isolatedprolines_with_length.txt", "a")
f.write(f"{effective_length} {(df['Fraction_Proline1']==0).sum()}\n")
f.close()




#f = open("proline_double.txt", "a")
#f.write(f"{effective_length} {df['Fraction_Proline2'].mean()} {0.5* df['Fraction_Proline2'].std()}\n")
#f.close()
#f = open("proline_triple.txt", "a")
#f.write(f"{effective_length} {df['Fraction_Proline3'].mean()} {0.5* df['Fraction_Proline3'].std()}\n")
#f.close()

#f = open("proline_phe.txt", "a")
#f.write(f"{effective_length} {0.5*(df['Fraction_PhePro'].mean()+df['Fraction_ProPhe'].mean())} {0.25* df['Fraction_PhePro'].std()+0.25* df['Fraction_ProPhe'].std()}\n")
#f.close()


#f = open("proline_trp.txt", "a")
#f.write(f"{effective_length} {0.5*(df['Fraction_TrpPro'].mean()+df['Fraction_ProTrp'].mean())} {0.25* df['Fraction_TrpPro'].std()+0.25* df['Fraction_ProTrp'].std()}\n")
#f.close()

#f = open("proline_tyr.txt", "a")
#f.write(f"{effective_length} {0.5*(df['Fraction_TyrPro'].mean()+df['Fraction_ProTyr'].mean())} {0.25* df['Fraction_TyrPro'].std()+0.25* df['Fraction_ProTyr'].std()}\n")
#f.close()





#####################################################################
#f = open("proline_phe_new_norm.txt", "a")
#f.write(f"{effective_length} {0.5*df_frac_aromatic['tot_frac_pf'].mean()} {0.5* df_frac_aromatic['tot_frac_pf'].std()}\n")
#f.close()


#f = open("proline_trp_new_norm.txt", "a")
#f.write(f"{effective_length} {0.5*df_frac_aromatic['tot_frac_pw'].mean()} {0.5* df_frac_aromatic['tot_frac_pw'].std()}\n")
#f.close()


#f = open("proline_tyr_new_norm.txt", "a")
#f.write(f"{effective_length} {0.5*df_frac_aromatic['tot_frac_py'].mean()} {0.5* df_frac_aromatic['tot_frac_py'].std()}\n")
#f.close()
#df_frac_aromatic.to_csv("aromatic_data2.csv")
print(df_frac_aromatic.mean())
#df_frac_aromatic.sort_values(by =["tot_frac_pf","tot_frac_pw", "tot_frac_py"], ascending=False).to_csv("temp3.csv")




#df_frac_aromatic.sort_values(by ="tot_frac_pw", ascending=False)
#df_frac_aromatic.sort_values(by ="tot_frac_ps", ascending=False)
#print(df_frac_aromatic['tot_frac_pf'].idxmax(), df_frac_aromatic['tot_frac_pw'].idxmax(), df_frac_aromatic['tot_frac_py'].idxmax(), df_frac_aromatic['tot_pro'].idxmax(), df_frac_aromatic['seq_length'].idxmax())



#########################################plotting contours#################################################################

#plot=sns.kdeplot( data=df_frac_aromatic, x='seq_length', y='tot_frac_pf',color='r', fill=True,cmap='Reds', cbar =True)
#plot=sns.kdeplot( data=df, x='seq_length', y='Fraction_Proline1',color='r', fill=True,cmap='Reds', cbar =True)
#plot=sns.jointplot( data=df_frac_aromatic, x='seq_length', y='tot_frac_all_aromatic',kind='reg')
#plt.xlim(0, 300)
#plt.ylim( 0, 30)
#plt.savefig("proline1.png",bbox_inches = "tight" )
#plt.show()
#df_frac_aromatic.plot.scatter(x = 'seq_length', y = 'tot_frac_py', s = 100)
#dz = px.data.df_frac_aromatic()


#fig=px.density_contour(dz, x="seq_length", y="tot_frac_py")
#fig.show()
#fig = px.density_contour(df_frac_aromatic, x="seq_length", y="tot_frac_pf")
#fig.update_traces(contours_coloring="fill", contours_showlabels = True)
#fig.show()

