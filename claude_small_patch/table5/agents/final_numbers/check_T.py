import json, statistics as st
R = json.load(open('/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/final_numbers/myparse.json'))
meshes = [(20,20),(40,40),(80,80),(80,160),(160,160)]
out = {}
scat = []
for nx, nz in meshes:
    for s in ('fs','fspmat','monofix'):
        tags = [f'T_{s}_rep{k}_np16_nx{nx}_nz{nz}' for k in (1,2,3)]
        pcs = [R[t]['off_pc'][1] for t in tags]; npc = [R[t]['off_pc'][0] for t in tags]
        ks = [R[t]['off_ksp'][1] for t in tags]
        its = [R[t]['its'] for t in tags]
        med = st.median(pcs)
        out[(nx,nz,s)] = dict(its=its, npc=npc, pcs=pcs, med=med, mn=min(pcs), mx=max(pcs), kmed=st.median(ks),
                              per=med/its[0][-1], sc=(max(pcs)-min(pcs))/med, patch=[R[t]['patch'] for t in tags])
        scat.append(((max(pcs)-min(pcs))/med, nx, nz, s))
        o = out[(nx,nz,s)]
        print(f"({nx},{nz}) {s:8s} its {['/'.join(map(str,i)) for i in its]} nPC {npc} PCApply {pcs} med {med:.4f} [{min(pcs):.4f},{max(pcs):.4f}] "
              f"s/it {o['per']:.4f} KSP med {o['kmed']:.4f} scatter {100*o['sc']:.2f}% patch {o['patch']}")
print()
for nx, nz in meshes:
    f, m, p = out[(nx,nz,'fs')], out[(nx,nz,'monofix')], out[(nx,nz,'fspmat')]
    print(f"({nx},{nz}) per-app mono/fs-1 = {100*(m['per']/f['per']-1):.2f}%  cum mono/fs-1 = {100*(m['med']/f['med']-1):.2f}%  "
          f"1-fs/mono = {100*(1-f['med']/m['med']):.2f}%  fspmat/fs-1 = {100*(p['med']/f['med']-1):.2f}%  "
          f"ranges overlap fs/mono: {not (m['mn']>f['mx'] or f['mn']>m['mx'])}  fs/fspmat: {not (p['mn']>f['mx'] or f['mn']>p['mx'])}")
scat.sort()
print('scatter sorted:', [f"{100*a:.2f}% ({b},{c},{d})" for a,b,c,d in scat])
print('median scatter', 100*st.median([a for a,_,_,_ in scat]), 'max', 100*max(a for a,_,_,_ in scat))
print('scaling 160/80,160 fs', out[(160,160,'fs')]['med']/out[(80,160,'fs')]['med'], 'mono', out[(160,160,'monofix')]['med']/out[(80,160,'monofix')]['med'])
for nz in (20,40,80,160):
    print('patch ratio', nz, (30*nz+8)/(38*nz+8), 1-(30*nz+8)/(38*nz+8))
