"""Reanalyse supplied RNA-seq counts using ENA run annotations, not inferred SRR order.

The count quantification itself has not been regenerated from FASTQ. Results
are a corrective analysis conditional on these supplied counts.
"""
import os
os.environ['OMP_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import json, sys, re
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pydeseq2

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/derived'

def main():
    meta=pd.read_csv(ROOT/'references/PRJNA800740_ENA_run_metadata.tsv',sep='\t').set_index('run_accession')
    meta['genotype']=np.where(meta.scientific_name.eq('Solanum pennellii'),'Sp','M82')
    meta['treatment']=meta.sample_title.str.extract('(control|drought)',expand=False)
    meta['day']=meta.sample_title.str.extract('day (\d+)',expand=False).astype(int)
    meta['replicate']=meta.experiment_title.str.extract('replicate (\d+)',expand=False).astype(int)
    meta['group']=meta.genotype+'_'+meta.treatment+'_D'+meta.day.astype(str)
    c=pd.read_csv(ROOT/'data/source/PRJNA800740_Solyc_expr_counts.csv',index_col=0)
    c.index=c.index.str.replace(r'\.\d+$','',regex=True);c=c.groupby(level=0).sum().T
    assert np.all(np.isfinite(c)) and np.all(c>=0) and np.all(c==c.round())
    c=c.astype(int);meta=meta.loc[c.index]
    assert meta.groupby(['genotype','treatment','day']).size().eq(3).all()
    meta.to_csv(OUT/'leaf_metadata_corrected.tsv',sep='\t')
    old=pd.read_csv(ROOT/'data/source/PRJNA800740_M82_metadata.tsv',sep='\t').set_index('sample_id')
    audit=old.join(meta,rsuffix='_ENA')
    audit['species_match']=audit.genotype.eq('M82')
    audit['treatment_match']=audit.treatment.eq(audit.treatment_ENA)
    audit['timepoint_match']=audit.timepoint.eq('D'+audit.day.astype(str))
    audit['replicate_match']=audit.replicate.eq(audit.replicate_ENA)
    audit['all_four_match']=audit[['species_match','treatment_match','timepoint_match','replicate_match']].all(axis=1)
    audit.to_csv(ROOT/'audit/leaf_metadata_conflicts.tsv',sep='\t')
    keep=c.ge(10).sum(axis=0)>=3
    filtered=c.loc[:,keep]
    design=pd.get_dummies(meta.group,dtype=float)
    design.to_csv(OUT/'leaf_design_matrix.tsv',sep='\t')
    dds=DeseqDataSet(counts=filtered,metadata=meta,design=design,refit_cooks=True,n_cpus=4,quiet=False)
    print('Fit',filtered.shape,'rank',np.linalg.matrix_rank(design),flush=True)
    dds.deseq2()
    norm=pd.DataFrame(dds.layers['normed_counts'],index=c.index,columns=filtered.columns)
    norm.to_csv(OUT/'leaf_deseq2_normalized_counts.tsv',sep='\t')
    pd.Series(dds.obs['size_factors'],index=c.index,name='size_factor').to_csv(OUT/'leaf_size_factors.tsv',sep='\t')
    log=np.log2(norm+1)
    topvar=log.var().nlargest(min(2000,len(log.columns))).index
    pca=PCA(n_components=2);xy=pca.fit_transform(log[topvar])
    pd.DataFrame(xy,index=log.index,columns=['PC1','PC2']).join(meta).to_csv(OUT/'leaf_pca_corrected.tsv',sep='\t')
    # No cross-platform comparison of absolute expression levels.
    gmt={}
    for line in (ROOT/'data/source/pathway_genesets_tomato.gmt').read_text().splitlines():
        a=line.split('\t');gmt[a[0]]=a[2:]
    z=(log-log.mean())/log.std(ddof=1).replace(0,np.nan)
    pw=[]
    for name,genes in gmt.items():
        covered=sorted(set(genes)&set(z.columns))
        if not covered:continue
        s=z[covered].mean(axis=1)
        for sample,value in s.items():
            pw.append(dict(sample=sample,pathway=name,score=float(value),n_genes=len(covered),**meta.loc[sample,['genotype','treatment','day','replicate']].to_dict()))
    pw=pd.DataFrame(pw);pw.to_csv(OUT/'leaf_pathway_scores_corrected.tsv',sep='\t',index=False)
    # Effects averaged equally over the three timepoints; panel data retain all timepoints.
    cols=list(design.columns);contrasts={}
    def vec(weights):return np.array([weights.get(x,0.) for x in cols])
    contrasts['ddE_average']=vec({f'{g}_{t}_D{d}':s/3 for g,t,s in [('Sp','drought',1),('Sp','control',-1),('M82','drought',-1),('M82','control',1)] for d in [5,8,11]})
    for genotype in ['Sp','M82']:
        contrasts[f'dE_{genotype}_average']=vec({f'{genotype}_{t}_D{d}':s/3 for t,s in [('drought',1),('control',-1)] for d in [5,8,11]})
    for day in [5,8,11]:
        contrasts[f'ddE_D{day}']=vec({f'{g}_{t}_D{day}':s for g,t,s in [('Sp','drought',1),('Sp','control',-1),('M82','drought',-1),('M82','control',1)]})
    results={}
    for name,v in contrasts.items():
        print('Contrast',name,flush=True)
        ds=DeseqStats(dds,contrast=v,independent_filter=False,cooks_filter=True,n_cpus=4,quiet=True)
        ds.summary();r=ds.results_df.copy();r.index.name='Solyc_ID'
        r.to_csv(OUT/f'leaf_corrected_{name}.tsv',sep='\t');results[name]=r
        print(name,'FDR<.05',int(r.padj.lt(.05).sum()),'RPE',r.loc['Solyc03g115820'].to_dict(),flush=True)
    root=pd.read_csv(ROOT/'data/source/ddE_interaction_GSE39894.tsv',sep='\t')
    root.Solyc_ID=root.Solyc_ID.str.replace(r'\.\d+$','',regex=True)
    cross=root.merge(results['ddE_average'].reset_index(),on='Solyc_ID',validate='one_to_one',suffixes=('_root','_leaf'))
    cross['same_sign']=np.sign(cross.ddE)==np.sign(cross.log2FoldChange)
    cross['candidate_core_corrected']=cross.interaction_FDR.lt(.05)&cross.padj.lt(.05)&cross.same_sign
    cross['was_original_core']=cross.Solyc_ID.isin(pd.read_csv(ROOT/'data/source/systemic_core_genes.tsv',sep='\t').Solyc_ID)
    cross.to_csv(OUT/'root_leaf_corrected_comparison.tsv',sep='\t',index=False)
    fact=dict(pydeseq2_version=pydeseq2.__version__,input_genes=c.shape[1],analysed_genes=filtered.shape[1],samples=len(c),
        formula='12-cell means: genotype x treatment x day; 3 samples per cell',n_parameters=len(cols),
        filter='count >=10 in at least 3 of 36 samples',normalization='DESeq2 median-of-ratios',
        pca_variance_ratio=pca.explained_variance_ratio_.tolist(),
        original_metadata_fully_correct_rows=int(audit.all_four_match.sum()),original_metadata_rows=len(audit),
        corrected_candidate_core=int(cross.candidate_core_corrected.sum()),
        original_core_retained=int((cross.candidate_core_corrected&cross.was_original_core).sum()),
        caveat='Corrective analysis of supplied count matrix. FASTQ mapping, genotype-specific mapping bias and orthology assignments have not been revalidated.')
    (ROOT/'audit/leaf_reanalysis.json').write_text(json.dumps(fact,indent=2),encoding='utf-8')
    print(json.dumps(fact,indent=2),flush=True)

if __name__=='__main__':main()
