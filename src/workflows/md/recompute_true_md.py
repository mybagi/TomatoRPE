"""Audit real 100 ns project trajectories and reproduce molecular observables.

Phase frames extracts actual 0/50/100 ns coordinates from existing whole/PBC-
prepared trajectories, already fitted on RPE CA residues 48-279. No ligand fit,
fabricated coordinates, interpolation, or movement towards the binding pocket.
Phase metrics recomputes core CA RMSD and core-CA-fit ligand HEAVY-atom RMSD
from the prepared trajectories. Periodic distances use original raw TPR/XTC
coordinates with their original box. Never apply this box to rotated frames.
"""
from pathlib import Path
import argparse, hashlib, json, warnings, os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from pymbar import timeseries

WS=Path(os.environ.get("TOMATO_MD_WORKSPACE", ".")).resolve()
ROOT=WS/'outputs/fig1_structures_md_real_20261005'
OUT=ROOT/'md'; PREP=WS/'outputs/fig3-6_redrawn_comprehensive_20260917/data/md_inputs'
RAW=WS/'work/source_md_20260917'
OLD=WS/'outputs/fig6_focused_20260919/data/derived'
SYSTEMS=[('CNP0542861.1','Indole alkaloid'),('CNP0511659.1','Indolmycin'),
 ('CNP0184442.0','Scopoletin'),('CNP0196064.0','IBA'),('CNP0211069.0','Bergapten'),
 ('substrate','Substrate'),('apo','Apo control')]
TIMES=[0.,50.,100.]
warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()

def loaded(sid):
 src=PREP/sid; mp=pd.read_csv(src/'atom_mapping.tsv',sep='\t')
 u=mda.Universe(str(src/'complex.pdb'),str(src/'trajectory.xtc'))
 assert len(u.atoms)==len(mp)
 assert np.array_equal(u.atoms.names,mp['name'].to_numpy())
 assert np.array_equal(u.atoms.resids,mp.residue_number.to_numpy())
 assert np.array_equal(u.atoms.resnames,mp.residue_name.to_numpy())
 p=np.flatnonzero(mp.chain.eq('A')); l=np.flatnonzero(mp.chain.eq('L'))
 ca=p[(mp.loc[p,'name']=='CA').to_numpy()]
 core=ca[mp.loc[ca,'residue_number'].between(48,279).to_numpy()]
 lh=l[mp.loc[l,'element'].ne('H').to_numpy()]
 assert len(ca)==282 and len(core)==232
 assert (sid=='apo')==(len(l)==0)
 return u,mp,p,l,ca,core,lh

def frames():
 rows=[]; systems=[]
 for sid,name in SYSTEMS:
  src=PREP/sid; u,mp,p,l,ca,core,lh=loaded(sid)
  dest=OUT/'frames'/sid;dest.mkdir(parents=True,exist_ok=True)
  meta=json.loads((src/'preparation.json').read_text(encoding='utf-8'))
  raw=mda.Universe(str(RAW/f'complex_{sid}/md.tpr'),str(RAW/f'complex_{sid}/md.xtc'))
  rt=[float(raw.trajectory[0].time)/1000,float(raw.trajectory[-1].time)/1000]
  assert len(raw.trajectory)==len(u.trajectory)==1001 and np.allclose(rt,[0,100])
  assert np.allclose([u.trajectory[0].time/1000,u.trajectory[-1].time/1000],rt)
  for t in TIMES:
   i=int(round(t*1000/u.trajectory.dt));ts=u.trajectory[i]
   assert abs(ts.time/1000-t)<1e-8
   path=dest/f'{int(t):03d}ns.pdb';u.atoms.write(str(path),bonds='all')
   check=mda.Universe(str(path));err=float(np.max(np.abs(check.atoms.positions-ts.positions)))
   assert err<.0006
   old=WS/f'outputs/fig3-6_refined_20260918/structures/frames/{sid}/{int(t):03d}ns.pdb'
   olderr=float(np.max(np.abs(mda.Universe(str(old)).atoms.positions-check.atoms.positions))) if old.exists() else None
   if olderr is not None:assert olderr<.0011
   rows.append(dict(system=sid,name=name,time_ns=t,frame_index_zero_based=i,pdb=str(path),
    sha256=digest(path),prepared_source_topology=str(src/'complex.pdb'),
    prepared_source_trajectory=str(src/'trajectory.xtc'),raw_topology=str(RAW/f'complex_{sid}/md.tpr'),
    raw_trajectory=str(RAW/f'complex_{sid}/md.xtc'),n_atoms=len(u.atoms),n_core_ca=len(core),
    n_ligand_heavy=len(lh),max_pdb_rounding_error_A=err,existing_frame_max_difference_A=olderr,
    coordinate_convention='whole molecules; nearest periodic ligand image; rigid protein core CA48-279 fit to AlphaFold RPE reference; no ligand self-fit'))
  systems.append(dict(system=sid,name=name,raw_frames=len(raw.trajectory),raw_start_ns=rt[0],raw_end_ns=rt[1],
   dt_ns=raw.trajectory.dt/1000,raw_atoms=len(raw.atoms),solute_atoms=len(u.atoms),protein_residues=282,
   ligand_atoms=len(l),ligand_heavy_atoms=len(lh),original_preparation=meta))
  print('FRAMES_READY',sid,flush=True)
 pd.DataFrame(rows).to_csv(OUT/'frame_manifest.tsv',sep='\t',index=False)
 (OUT/'system_inputs.json').write_text(json.dumps(systems,indent=2),encoding='utf-8')
 print('ALL_FRAMES_READY',str(OUT/'frame_manifest.tsv'),flush=True)

def fit(x,ref):
 cx=x.mean(0);cr=ref.mean(0);a=x-cx;b=ref-cr
 uu,ss,vt=np.linalg.svd(a.T@b);d=np.eye(3);d[-1,-1]=np.linalg.det(uu@vt)
 rot=uu@d@vt
 return rot,cx,cr,float(np.sqrt(np.mean(np.sum((a@rot-b)**2,axis=1)))/10)

def frame_contacts():
 records=[];checks=[]
 pocket=pd.read_csv(OLD/'fixed_pocket_definition.tsv',sep='\t').resid.to_numpy(int)
 for sid,name in SYSTEMS[:-1]:
  mp=pd.read_csv(PREP/sid/'atom_mapping.tsv',sep='\t')
  ph=np.flatnonzero(mp.chain.eq('A') & mp.element.ne('H'))
  lh=np.flatnonzero(mp.chain.eq('L') & mp.element.ne('H'))
  metrics=pd.read_csv(OUT/f'md_metrics_{sid}.tsv',sep='\t')
  for t in TIMES:
   source=OUT/f'frames/{sid}/{int(t):03d}ns.pdb';u=mda.Universe(str(source))
   d=distance_array(u.atoms.positions[ph],u.atoms.positions[lh])
   row=metrics[metrics.time_ns.eq(t)].iloc[0]
   minimum=float(d.min()/10);expected=float(row.protein_ligand_heavy_minimum_image_distance_nm)
   checks.append(dict(system=sid,time_ns=t,render_coordinate_heavy_min_nm=minimum,
     raw_box_minimum_image_heavy_min_nm=expected,absolute_difference_nm=abs(minimum-expected),
     definition='Direct geometry in whole/core-fitted closest-image render frame; not a continuous diffusion trajectory'))
   for rid in np.unique(mp.loc[ph,'residue_number']):
    local=np.flatnonzero(mp.loc[ph,'residue_number'].to_numpy()==rid)
    ai,bi=np.unravel_index(d[local].argmin(),d[local].shape)
    pi=ph[local[ai]];li=lh[bi];dist=float(d[local[ai],bi])
    if dist>5:continue
    records.append(dict(system=sid,time_ns=t,source=str(source),resid=int(rid),
      resname=mp.iloc[pi].residue_name,protein_atom=mp.iloc[pi]['name'],ligand_atom=mp.iloc[li]['name'],
      protein_index=int(pi),ligand_index=int(li),distance_A=dist,fixed_pocket=bool(rid in pocket),
      contact_definition='Minimum residue-ligand heavy-atom distance in this actual PDB; proximity only, not hydrogen-bond assignment'))
 pd.DataFrame(records).to_csv(OUT/'frame_interface_contacts.tsv',sep='\t',index=False)
 pd.DataFrame(checks).to_csv(OUT/'frame_render_geometry_crosscheck.tsv',sep='\t',index=False)
 auditpath=ROOT/'audit/md_audit.json';report=json.loads(auditpath.read_text(encoding='utf-8'))
 report['render_frame_heavy_distance_crosscheck']=checks
 report['render_contact_table']=str(OUT/'frame_interface_contacts.tsv')
 auditpath.write_text(json.dumps(report,indent=2),encoding='utf-8')
 print('CONTACTS_READY',len(records),'render geometry max difference nm',max(x['absolute_difference_nm'] for x in checks),flush=True)

def metrics():
 OUT.mkdir(parents=True,exist_ok=True);(ROOT/'audit').mkdir(parents=True,exist_ok=True)
 inputs=json.loads((OUT/'system_inputs.json').read_text(encoding='utf-8'))
 fixed=pd.read_csv(OLD/'fixed_pocket_definition.tsv',sep='\t');pocket=fixed.resid.to_numpy(int)
 summaries=[];audits=[];sampling=[]
 for sid,name in SYSTEMS:
  src=PREP/sid;u,mp,p,l,ca,core,lh=loaded(sid)
  u.trajectory[0];ref=u.atoms.positions.astype(float).copy();rows=[]
  previous=pd.read_csv(OLD/f'plot_md_{sid}.tsv',sep='\t')
  for i,ts in enumerate(u.trajectory):
   x=u.atoms.positions.astype(float);rot,cx,cr,core_rmsd=fit(x[core],ref[core]);_,_,_,full_rmsd=fit(x[ca],ref[ca])
   row=dict(system=sid,frame_index_zero_based=i,time_ns=float(ts.time)/1000,
    core_ca_rmsd_nm=core_rmsd,full_protein_ca_rmsd_nm=full_rmsd)
   if len(lh):
    aligned=(x[lh]-cx)@rot+cr
    row['ligand_heavy_core_ca_fit_rmsd_nm']=float(np.sqrt(np.mean(np.sum((aligned-ref[lh])**2,axis=1)))/10)
   rows.append(row)
  df=pd.DataFrame(rows);assert len(df)==1001 and np.allclose(df.time_ns,np.arange(1001)*.1)
  full_error=float(np.max(np.abs(df.full_protein_ca_rmsd_nm-previous.full_ca_rmsd_nm)))
  assert full_error<1e-6
  rawtop=RAW/f'complex_{sid}/md.tpr';rawtraj=RAW/f'complex_{sid}/md.xtc'
  meta=json.loads((src/'preparation.json').read_text(encoding='utf-8'))
  tophash=digest(rawtop);trajhash=digest(rawtraj)
  assert tophash==meta['source_topology_sha256'] and trajhash==meta['source_trajectory_sha256']
  raw=mda.Universe(str(rawtop),str(rawtraj));orig=mp.source_atom.to_numpy(int)-1
  assert np.array_equal(raw.atoms[orig].names,mp['name']) and np.array_equal(raw.atoms[orig].resnames,mp.residue_name)
  raw_residue_counts=pd.Series(raw.atoms.resnames).value_counts().to_dict()
  catalytic_metal_names=[str(n) for n in raw_residue_counts if str(n).upper() in ['ZN','ZN2','MN','MG','FE','CO','NI','CU']]
  audit=dict(system=sid,name=name,source_topology_sha256=tophash,source_trajectory_sha256=trajhash,
   original_hash_match=True,prepared_topology_sha256=digest(src/'complex.pdb'),prepared_trajectory_sha256=digest(src/'trajectory.xtc'),
   raw_frames=len(raw.trajectory),raw_start_ns=raw.trajectory[0].time/1000,raw_end_ns=raw.trajectory[-1].time/1000,
   original_atom_mapping_verified=True,full_ca_vs_previous_max_abs_difference_nm=full_error,
   raw_topology_residue_atom_counts=raw_residue_counts,catalytic_metal_resnames_found=catalytic_metal_names)
  if len(lh):
   ph=p[mp.loc[p,'element'].ne('H').to_numpy()];rd=mp.loc[ph,'residue_number'].to_numpy(int)-1
   rp=raw.atoms[orig[ph]];rl=raw.atoms[orig[lh]];ispocket=np.isin(rd+1,pocket)
   distance_rows=[];charge=float(raw.atoms[orig[l]].charges.sum())
   resdist=np.empty((len(raw.trajectory),282),np.float32)
   for i,ts in enumerate(raw.trajectory):
    d=distance_array(rp.positions,rl.positions,box=ts.dimensions)/10
    rr=np.full(282,np.inf);np.minimum.at(rr,rd,d.min(1));resdist[i]=rr
    distance_rows.append(dict(protein_ligand_heavy_minimum_image_distance_nm=float(d.min()),
      fixed_pocket_heavy_minimum_image_distance_nm=float(d[ispocket].min()),
      protein_contact_residues_below_0_4nm=int((rr<.4).sum()),
      fixed_pocket_contact_residues_below_0_4nm=int((rr[pocket-1]<.4).sum())))
    if i%250==0: print('PERIODIC_DISTANCE_RECOMPUTED',sid,i,flush=True)
   ds=pd.DataFrame(distance_rows);df=pd.concat([df,ds],axis=1)
   oldfixed=pd.read_csv(OLD/f'fixed_pocket_{sid}.tsv',sep='\t')
   pockerr=float(np.max(np.abs(df.fixed_pocket_heavy_minimum_image_distance_nm-oldfixed.pocket_min_nm)))
   minerr=float(np.max(np.abs(df.protein_ligand_heavy_minimum_image_distance_nm-previous.mindist_heavy_nm)))
   assert pockerr<1e-6 and minerr<1e-6
   assert np.array_equal(df.fixed_pocket_contact_residues_below_0_4nm,oldfixed.pocket_contact_residues_4A)
   np.savez_compressed(OUT/f'residue_minimum_image_distances_{sid}.npz',time_ns=df.time_ns.to_numpy(),residue=np.arange(1,283),minimum_distance_nm=resdist)
   audit.update(ligand_heavy_atoms=len(lh),ligand_atoms=len(l),ligand_total_charge_e=charge,
    raw_minimum_image_distance_vs_previous_max_abs_difference_nm=minerr,
    fixed_pocket_vs_previous_max_abs_difference_nm=pockerr,fixed_pocket_contact_counts_equal=True,
    ligand_rmsd_change='New observable: heavy atoms and core CA fit. Prior ligand_corefit_rmsd_nm included ligand H and used core backbone fitting; do not relabel it as heavy-atom RMSD.')
  df.to_csv(OUT/f'md_metrics_{sid}.tsv',sep='\t',index=False,float_format='%.10g')
  late=df.time_ns>=50
  for col in ['core_ca_rmsd_nm']+(['ligand_heavy_core_ca_fit_rmsd_nm'] if len(lh) else []):
   xx=df.loc[late,col].to_numpy(float);tt=df.loc[late,'time_ns'].to_numpy(float)
   sd=float(xx.std(ddof=1));constant=bool(np.ptp(xx)<1e-12)
   gg=float(timeseries.statistical_inefficiency(xx,fast=False,mintime=3)) if not constant else None
   trend=float(np.polyfit(tt-tt[0],xx,1)[0]*(tt[-1]-tt[0])/sd) if not constant else None
   shift=float((xx[len(xx)//2:].mean()-xx[:len(xx)//2].mean())/sd) if not constant else None
   sampling.append(dict(system=sid,observable=col,start_ns=50,end_ns=100,n_frames=len(xx),
    statistical_inefficiency_g=gg,conditional_neff=float(len(xx)/gg) if gg else None,
    correlation_time_ns=(gg-1)*.1/2 if gg else None,linear_drift_in_sd=trend,half_mean_shift_in_sd=shift,
    descriptive_drift_flag=bool(abs(trend)>1 or abs(shift)>.5) if trend is not None else False,
    constant=constant,independent_trajectories=1,interpretation='Conditional observable ESS assuming stationarity; not a convergence certificate, not biological n'))
  s=dict(system=sid,name=name,n_frames=len(df),start_ns=df.time_ns.iloc[0],end_ns=df.time_ns.iloc[-1],
    core_ca_rmsd_mean_50_100_ns_nm=float(df.loc[late,'core_ca_rmsd_nm'].mean()),core_ca_rmsd_end_nm=float(df.core_ca_rmsd_nm.iloc[-1]))
  if len(lh):s.update(ligand_heavy_rmsd_end_nm=float(df.ligand_heavy_core_ca_fit_rmsd_nm.iloc[-1]),
   pocket_distance_end_nm=float(df.fixed_pocket_heavy_minimum_image_distance_nm.iloc[-1]),
   pocket_contact_fraction_0_100_ns_below_0_5nm=float((df.fixed_pocket_heavy_minimum_image_distance_nm<.5).mean()))
  summaries.append(s);audits.append(audit);print('METRICS_READY',sid,flush=True)
 pd.DataFrame(summaries).to_csv(OUT/'md_summary.tsv',sep='\t',index=False)
 pd.DataFrame(sampling).to_csv(OUT/'sampling_neff_50_100ns.tsv',sep='\t',index=False)
 fixed.to_csv(OUT/'fixed_pocket_definition.tsv',sep='\t',index=False)
 report=dict(software={'MDAnalysis':mda.__version__,'numpy':np.__version__,'pandas':pd.__version__},
  systems=audits,frames_per_system=1001,spacing_ns=.1,times_ns=TIMES,independent_trajectories_per_system=1,
  observable_definitions={'core_ca_rmsd_nm':'232 Cα atoms, RPE residues48–279; unweighted rigid fit to the same systems actual first trajectory frame; Å /10 to nm',
   'full_protein_ca_rmsd_nm':'282 Cα atoms; independent unweighted full-CA fit to same-system first trajectory frame',
   'ligand_heavy_core_ca_fit_rmsd_nm':'Indexed ligand non-H atoms; rigid transform from protein core Cα; same-system first-frame reference; no ligand self-alignment',
   'distances':'Minimum-image heavy-atom pair distances recalculated from original raw TPR/XTC with original triclinic box, nm',
   'contacts':'Number of residues whose minimum ligand-heavy distance is <0.4nm; geometry not hydrogen-bond assignment'},
  rendering_geometry='Existing whole-molecule/no-jump processing; closest periodic ligand image; protein core CA alignment to common AlphaFold RPE. This is not a continuous ligand diffusion path. Does not pull escaping ligand into pocket.',
  conditional_sampling_diagnostics=sampling,
  neff_method='PyMBAR 4.0.3 timeseries.statistical_inefficiency fast=False,mintime=3; N_eff=N/g for fixed50–100ns window501frames; constant observable g/ESS undefined; descriptive drift flags only, no convergence assertion, no CI/P value',
  limitations=['0/50/100ns snapshots are descriptive frames, not independent samples','Each system has one trajectory; timepoints are not biological replicates','Conditional effective sample sizes depend on observable and stationarity; equal temporal spacing does not establish independence','Protein is a computational model; no catalytic metal present in project systems','Topology file/system ID linkage does not replace independent compound identity validation'],
  completed='Fresh per-frame RMSD and raw periodic distance/contact recomputation; original hashes and previous audited metrics crosschecked')
 (ROOT/'audit/md_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print('MD_AUDIT_READY',str(ROOT/'audit/md_audit.json'),flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--phase',choices=['frames','metrics','contacts','all'],default='all');args=a.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 if args.phase in ['frames','all']:frames()
 if args.phase in ['metrics','all']:metrics()
 if args.phase in ['contacts','all']:frame_contacts()
