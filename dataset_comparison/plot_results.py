from audit import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from textwrap import shorten
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
POLICIES=['PSI_original','PSI_NetSolP','NESG_unique_E45_S45_vs0','NESG_NetSolP']
TITLES=['PSI original (12,216)','NetSolP PSI (11,226)','NESG strict scores (3,762)','NetSolP NESG (1,323)']
def main():
 figures=ROOT/'figures';figures.mkdir(exist_ok=True);d=read(OUT/'family_clan_prevalence.tsv.gz')
 for level in ['family','clan']:
  fig,axs=plt.subplots(2,2,figsize=(17,13));fig.subplots_adjust(left=.22,right=.98,wspace=.95,bottom=.07,top=.93,hspace=.33)
  for ax,pol,title in zip(axs.flat,POLICIES,TITLES):
   g=d[(d.policy==pol)&(d.level==level)&(d.combined_rank<=20)];ids=g.sort_values(['combined_rank','accession']).accession.drop_duplicates().tolist();y=np.arange(len(ids));names=g.drop_duplicates('accession').set_index('accession')
   for lab,delta,color,label in [(1,-.18,'#8270ad','Soluble'),(0,.18,'#b2bac8','Insoluble')]:
    v=g[g.label==lab].set_index('accession').percent.reindex(ids,fill_value=0);ax.barh(y+delta,v,.35,color=color,label=label)
   ax.set_yticks(y,[shorten(str(names.loc[x,'description'] or names.loc[x,'name']),width=43,placeholder='…')+' ['+x+']' for x in ids],fontsize=7);ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('Proteins with this assignment (% of class)');ax.legend(fontsize=8)
  fig.suptitle('Actual Pfam '+('families' if level=='family' else 'clans/superfamilies')+' — top 20 by combined protein prevalence',fontsize=14)
  fig.text(.5,.01,'A protein counts once per category, but can carry several categories. These are prevalence bars, not a partition summing to 100%.',ha='center',fontsize=9)
  for ext in ['png','svg']:fig.savefig(figures/f'top20_{level}.{ext}',dpi=200)
  plt.close(fig)
 f=read(OUT/'family_support.tsv.gz');fig,ax=plt.subplots(figsize=(7,4.5))
 for pol,title in zip(POLICIES,TITLES):
  z=f[(f.policy==pol)&(f['filter']=='all_GA')].proteins.sort_values(ascending=False).to_numpy();ax.loglog(np.arange(1,len(z)+1),z,label=title)
 ax.set_xlabel('Family rank');ax.set_ylabel('Protein records per family');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(figures/'family_long_tail.png',dpi=200);fig.savefig(figures/'family_long_tail.svg');plt.close(fig)
 e=pd.read_csv(OUT/'NESG_E_by_S.tsv',sep='\t').set_index('E');fig,ax=plt.subplots(figsize=(6,6));im=ax.imshow(np.log1p(e.to_numpy()),cmap='Purples');ax.set_xticks(range(6),e.columns);ax.set_yticks(range(len(e)),e.index);ax.set_xlabel('Solubility score S (as distributed)');ax.set_ylabel('Expression score E');ax.set_title('NESG archive: all 9,703 raw records')
 for i in range(len(e)):
  for j in range(6):ax.text(j,i,str(e.iloc[i,j]),ha='center',va='center',fontsize=8,color='white' if e.iloc[i,j]>150 else '#17273c')
 fig.text(.5,.01,'E = 0 is not interpreted as measured insolubility; fractional E is unresolved.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.03,1,1]);fig.savefig(figures/'NESG_score_contingency.png',dpi=200);plt.close(fig)
if __name__=='__main__':main()
